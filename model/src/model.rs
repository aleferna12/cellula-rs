use crate::cell::Cell;
use crate::cellular_automata::CellularAutomata;
use crate::chem_environment::ChemEnvironment;
use crate::clonal_adhesion::ClonalAdhesion;
use crate::constants::{BoundaryType, NeighbourhoodType};
use crate::ecology::disperser::{Disperser, SelectiveDispersion};
use crate::ecology::transporter::{Transporter, WipeOut};
use crate::genetics::grn::Grn;
use crate::io::io_manager::IoManager;
use crate::io::movie_maker::MovieMaker;
use crate::io::parameters::Parameters;
use crate::pond::Pond;
use cellulars_lib::adhesion::StaticAdhesion;
use cellulars_lib::cell_container::CellContainer;
use cellulars_lib::environment::Environment;
use cellulars_lib::evolution::selector::WeightedOrderedSelection;
use cellulars_lib::lattice_entity::LatticeEntity;
use cellulars_lib::positional::boundaries::Boundaries;
use cellulars_lib::positional::rect::Rect;
use rand::distr::{Distribution, Uniform};
use rand::{RngCore, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use std::error::Error;

pub struct Model {
    pub ponds: Vec<Pond>,
    pub io: IoManager,
    pub rng: Xoshiro256StarStar,
    pub dispersion_period: u32,
    time_steps: u32
}

impl Model {
    pub fn initialise_from_parameters(mut parameters: Parameters) -> Result<Model, Box<dyn Error>> {
        log::info!("Initialising model");
        // TOML doesnt support large u64s so we use a u32 seed
        let seed = parameters.general.seed.unwrap_or(Xoshiro256StarStar::from_os_rng().next_u32() as u64);
        parameters.general.seed = seed.into();
        let mut rng = Xoshiro256StarStar::seed_from_u64(seed);

        let movie_maker = if parameters.io.movie.show {
            match MovieMaker::new(
                parameters.io.movie.width,
                parameters.io.movie.height,
                parameters.io.movie.frame_period
            ) {
                Ok(mm) => {
                    log::info!("Creating window for real-time movie display");
                    Some(mm)
                },
                Err(e) => {
                    log::warn!("Failed to initialise movie maker with error `{e}`");
                    None
                }
            }
        } else { None };
        let io = IoManager::builder()
            .outdir(parameters.io.outdir.clone().into())
            .image_format(parameters.io.image_format.clone())
            .image_period(parameters.io.image_period)
            .cell_period(parameters.io.cell_period)
            .genome_period(parameters.io.genome_period)
            .lattice_period(parameters.io.lattice_period)
            .plots(parameters.io.plot.clone())
            .maybe_movie_maker(movie_maker)
            .build();

        log::info!("Creating output directories and copy of parameter file");
        io.create_directories(parameters.io.replace_outdir)?;
        io.create_parameters_file(&parameters)?;

        let mut env = ChemEnvironment::new(
            Environment::new(
                CellContainer::default(),
                NeighbourhoodType::new(parameters.pond.neigh_r),
                Boundaries::new(BoundaryType::new(Rect::new(
                    (0., 0.).into(),
                    (parameters.pond.width as f32, parameters.pond.height as f32).into(),
                ))).expect("failed to create boundaries during initialisation, lattice size is too big")
            ).expect("failed to create environment during initialisation, lattice size is too big")
        );
        if parameters.pond.enclose {
            env.make_border(true, true, true, true);
        }

        let ca = CellularAutomata::builder()
            .boltz_t(parameters.ca.boltz_t)
            .size_lambda(parameters.ca.size_lambda)
            .chemotaxis_mu(parameters.ca.chemotaxis_mu)
            .enable_migration(parameters.cell.migrate)
            .adhesion(
                ClonalAdhesion::new(
                    parameters.cell.max_cells + LatticeEntity::first_cell_spin(),
                    parameters.ca.adhesion.clone_energy,
                    StaticAdhesion {
                        cell_energy: parameters.ca.adhesion.cell_energy,
                        medium_energy: parameters.ca.adhesion.medium_energy,
                        solid_energy: parameters.ca.adhesion.solid_energy,
                    }
                )
            )
            .build();

        let empty_pond = Pond::builder()
            .env(env)
            .ca(ca)
            .rng(rng.clone())
            .update_period(parameters.cell.update_period)
            .cell_target_area(parameters.cell.target_area)
            .cell_search_scaler(parameters.cell.search_radius)
            .division_enabled(parameters.cell.divide)
            .max_cells(parameters.cell.max_cells)
            .build();

        let mut ponds = vec![];
        for pond_i in 0..parameters.pond.n_ponds {
            log::info!("Making pond #{pond_i}");

            let mut pond = empty_pond.clone();
            for _ in 0..parameters.cell.starting_cells {
                let cell = Cell::new_empty(
                    parameters.cell.target_area,
                    parameters.cell.div_area,
                    Grn::new(
                        [1. / pond.env.height() as f32],
                        parameters.cell.genome.n_regulatory,
                        parameters.cell.genome.mutation_rate,
                        parameters.cell.genome.mutation_std,
                        || Uniform::new(-1., 1.).unwrap().sample(&mut rng)
                    )
                );
                pond.spawn_cell_random(
                    cell,
                    parameters.cell.starting_area
                );
            }
            log::info!(
                "Created {} out of the {} cells requested",
                pond.env.cells.n_valid(),
                parameters.cell.starting_cells
            );
            ponds.push(pond);
        }

        Ok(Self {
            ponds,
            io,
            rng,
            dispersion_period: parameters.general.dispersion_period,
            time_steps: parameters.general.time_steps
        })
    }

    pub fn run_for(&mut self, time_steps: u32) {
        for time_step in 0..=time_steps {
            let saved = self.io.try_io(
                time_step,
                &self.ponds
            );
            if let Err(e) = saved {
                log::warn!("Failed to save image at time step {time_step} with error `{e}`")
            }
            for pond in &mut self.ponds {
                pond.step();
            }

            if time_step > 0 && time_step % self.dispersion_period == 0 {
                let dispersed = SelectiveDispersion { 
                    selector: WeightedOrderedSelection{
                        rng: &mut self.rng 
                    } 
                }.disperse(&self.ponds);
                for event in dispersed {
                    let [from, to] = self.ponds
                        .get_disjoint_mut([event.from, event.to])
                        .expect("dispersion event `from` and `to` are the same");
                    WipeOut.transport(
                        from,
                        to,
                        event.spins
                    );
                }
            }
        }
    }

    pub fn run(&mut self) {
        self.run_for(self.time_steps);
    }
}

#[cfg(test)]
mod tests {
    use rand::{Rng, SeedableRng};
    use rand_xoshiro::Xoshiro256StarStar;

    #[test]
    fn test_seed() {
        let mut rng = Xoshiro256StarStar::seed_from_u64(1241254152);
        let s = (0..50)
            .map(|_| rng.random_range(0..9).to_string())
            .collect::<Vec<_>>()
            .join("");
        let res = "15515320360704325727185856564110164830043067488704";
        assert_eq!(res, s);
    }
}