use crate::adhesion::{ClonalAdhesion, StaticAdhesion};
use crate::cell::Cell;
use crate::cell_container::CellContainer;
use crate::cellular_automata::Ca;
use crate::constants::{BoundaryType, NeighbourhoodType};
use crate::ecology::disperser::{Disperser, SelectiveDispersion};
use crate::ecology::selector::WeightedOrderedSelection;
use crate::ecology::transporter::{Transporter, WipeOut};
use crate::environment::{Environment, LatticeEntity};
use crate::genetics::grn::Grn;
use crate::io::io_manager::IoManager;
use crate::io::movie_maker::MovieMaker;
use crate::io::parameters::Parameters;
use crate::pond::Pond;
use crate::positional::rect::Rect;
use crate::space::Space;
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

        let io = IoManager::new(
            &parameters.io.outdir,
            parameters.io.image_format.clone(),
            parameters.io.image_period,
            parameters.io.cell_period,
            parameters.io.genome_period,
            parameters.io.lattice_period,
            parameters.io.plot.clone(),
            if parameters.io.movie.show {
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
            } else {
                None
            }
        );

        log::info!("Creating output directories and copy of parameter file");
        io.create_directories(parameters.io.replace_outdir)?;
        io.create_parameters_file(&parameters)?;

        let mut ponds = vec![];
        // TODO: if we make everything clonable then that helps here and also in model_bench
        //  (takes less time to reinitialise everything which means more samples)
        for pond_i in 0..parameters.general.n_ponds {
            log::info!("Making pond #{pond_i}");
            let mut env = Environment::new(
                parameters.pond.cell_search_radius,
                parameters.pond.max_cells,
                CellContainer::new(
                    parameters.cell.target_area,
                    parameters.cell.divide,
                    parameters.cell.migrate,
                ),
                Space::new(
                    BoundaryType::new(Rect::new(
                        (0., 0.).into(),
                        (parameters.pond.width as f32, parameters.pond.height as f32).into(),
                    ))
                )?,
                NeighbourhoodType::new(parameters.pond.neigh_r)
            );

            if parameters.pond.enclose {
                env.make_border();
            }

            let mut pop_n = 0;
            for _ in 0..parameters.pond.starting_cells {
                let cell = Cell::new_empty(
                    parameters.cell.target_area,
                    parameters.cell.div_area,
                    Grn::new(
                        [1. / env.height() as f32],
                        parameters.cell.n_regulatory_genes,
                        parameters.cell.mutation_rate,
                        parameters.cell.mutation_std,
                        || Uniform::new(-1., 1.).unwrap().sample(&mut rng)
                    )
                );
                let spawned = env.spawn_cell_random(
                    parameters.cell.starting_area,
                    cell,
                    &mut rng
                );
                if spawned.is_some() {
                    pop_n += 1;
                }
            }
            log::info!("Created {pop_n} out of the {} cells requested", parameters.pond.starting_cells);

            let ca= Ca::new(
                parameters.cellular_automata.boltz_t,
                parameters.cellular_automata.size_lambda,
                parameters.cellular_automata.chemotaxis_mu,
                ClonalAdhesion::new(
                    parameters.pond.max_cells + LatticeEntity::first_cell_spin(),
                    StaticAdhesion {
                        cell_energy: parameters.cellular_automata.adhesion.cell_energy,
                        medium_energy: parameters.cellular_automata.adhesion.medium_energy,
                        solid_energy: parameters.cellular_automata.adhesion.solid_energy,
                    }
                )
            );
            ponds.push(Pond::new(
                env,
                ca,
                rng.clone(),
                parameters.cell.update_period,
            ));
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
                let dispersed = SelectiveDispersion{ 
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