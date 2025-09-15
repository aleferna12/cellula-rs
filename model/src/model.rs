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
use anyhow::Context;
use cellulars_lib::adhesion::StaticAdhesion;
use cellulars_lib::environment::Environment;
use cellulars_lib::evolution::selector::WeightedOrderedSelection;
use cellulars_lib::lattice_entity::LatticeEntity;
use cellulars_lib::positional::boundaries::Boundaries;
use cellulars_lib::positional::rect::Rect;
use cellulars_lib::symmetric_table::SymmetricTable;
use rand::distr::{Distribution, Uniform};
use rand::{RngCore, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use std::path::Path;

pub struct Model {
    pub ponds: Vec<Pond>,
    pub io: IoManager,
    pub rng: Xoshiro256StarStar,
    pub dispersion_period: u32,
    pub info_period: u32,
    time_steps: u32
}

impl Model {
    pub fn initialise_from_parameters(
        parameters: Parameters
    ) -> anyhow::Result<Self> {
        log::info!("Initialising model");

        let seed = Self::determine_seed(parameters.general.seed);
        let mut rng = Xoshiro256StarStar::seed_from_u64(seed);
        Ok(Self {
            ponds: Self::make_new_ponds(
                &parameters,
                Self::make_ca(&parameters, None),
                &mut rng
            )?,
            io: Self::setup_io(&parameters, seed)?,
            rng,
            dispersion_period: parameters.general.dispersion_period,
            info_period: parameters.io.info_period,
            time_steps: parameters.general.time_steps
        })
    }

    pub fn initialise_from_backup(
        parameters: Parameters,
        sim_path: impl AsRef<Path>,
        time_step: u32
    ) -> anyhow::Result<Self> {
        let sim_path = sim_path.as_ref();
        log::info!("Resuming simulation at {}", sim_path.display());
        log::info!("Starting from time step {time_step}");

        let seed = Self::determine_seed(parameters.general.seed);
        let mut rng = Xoshiro256StarStar::seed_from_u64(seed);
        Ok(Self {
            ponds: Self::read_backup_ponds(
                &parameters,
                &mut rng,
                sim_path,
                time_step
            )?,
            io: Self::setup_io(&parameters, seed)?,
            rng,
            dispersion_period: parameters.general.dispersion_period,
            info_period: parameters.io.info_period,
            time_steps: parameters.general.time_steps
        })
    }

    fn setup_io(parameters: &Parameters, new_seed: u64) -> anyhow::Result<IoManager> {
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
            .cells_period(parameters.io.data.cells_period)
            .genomes_period(parameters.io.data.genomes_period)
            .clones_period(parameters.io.data.clones_period)
            .lattices_period(parameters.io.data.lattices_period)
            .plots(parameters.io.plot.clone())
            .maybe_movie_maker(movie_maker)
            .build();

        log::info!("Creating output directories and copy of parameter file");
        if parameters.io.replace_outdir {
            log::info!("Cleaning contents of '{}'", io.outdir.display());
        }
        io.create_directories(parameters.io.replace_outdir, parameters.pond.n_ponds)?;
        let mut params_new_seed = parameters.clone();
        params_new_seed.general.seed = new_seed.into();
        io.create_parameters_file(&params_new_seed)?;
        Ok(io)
    }

    fn make_ca(parameters: &Parameters, clones: Option<SymmetricTable<bool>>) -> CellularAutomata<ClonalAdhesion> {
        CellularAutomata::builder()
            .boltz_t(parameters.ca.boltz_t)
            .size_lambda(parameters.ca.size_lambda)
            .chemotaxis_mu(parameters.ca.chemotaxis_mu)
            .enable_migration(parameters.cell.migrate)
            .adhesion(
                ClonalAdhesion::new(
                    parameters.ca.adhesion.clone_energy,
                    StaticAdhesion {
                        cell_energy: parameters.ca.adhesion.cell_energy,
                        medium_energy: parameters.ca.adhesion.medium_energy,
                        solid_energy: parameters.ca.adhesion.solid_energy,
                    },
                    clones.unwrap_or(SymmetricTable::new(
                        (parameters.cell.max_cells + LatticeEntity::first_cell_spin()) as usize)
                    )
                )
            )
            .build()
    }

    fn determine_seed(seed_param: Option<u64>) -> u64 {
        // TOML doesnt support large u64s so we use a u32 seed
        seed_param.unwrap_or(Xoshiro256StarStar::from_os_rng().next_u32() as u64)
    }

    fn make_empty_pond(
        parameters: &Parameters,
        env: ChemEnvironment,
        ca: CellularAutomata<ClonalAdhesion>,
        rng: &mut Xoshiro256StarStar
    ) -> Pond {
        Pond::builder()
            .env(env)
            .ca(ca)
            .rng(Xoshiro256StarStar::seed_from_u64(rng.next_u64()))
            .update_period(parameters.cell.update_period)
            .cell_target_area(parameters.cell.target_area)
            .cell_search_scaler(parameters.cell.search_radius)
            .division_enabled(parameters.cell.divide)
            .build()
    }

    fn make_new_ponds(
        parameters: &Parameters,
        ca: CellularAutomata<ClonalAdhesion>,
        rng: &mut Xoshiro256StarStar
    ) -> anyhow::Result<Vec<Pond>> {
        let mut env = ChemEnvironment::new(
            Environment::new_empty(
                NeighbourhoodType::new(parameters.pond.neigh_r),
                Boundaries::new(BoundaryType::new(Rect::new(
                    (0., 0.).into(),
                    (parameters.pond.width as f32, parameters.pond.height as f32).into(),
                ))).context("lattice size is too big")?
            ).context("lattice size is too big")?,
            parameters.cell.max_cells
        );
        if parameters.pond.enclose {
            env.make_border(true, true, true, true);
        }

        let mut ponds = vec![];
        for pond_i in 0..parameters.pond.n_ponds {
            log::info!("Making pond #{pond_i}");

            let mut pond = Self::make_empty_pond(parameters, env.clone(), ca.clone(), rng);
            for _ in 0..parameters.cell.starting_cells {
                let cell = Cell::new_empty(
                    parameters.cell.target_area,
                    parameters.cell.div_area,
                    Grn::from_sampler(
                        [1. / pond.env.height() as f32],
                        parameters.cell.genome.n_regulatory,
                        parameters.cell.genome.mutation_rate,
                        parameters.cell.genome.mutation_std,
                        || Uniform::new(-1., 1.).unwrap().sample(rng)
                    )
                );
                pond.env.spawn_cell_random(
                    cell,
                    parameters.cell.starting_area,
                    &mut pond.rng
                );
            }
            log::info!(
                "Created {} out of the {} cells requested",
                pond.env.cells.n_valid(),
                parameters.cell.starting_cells
            );
            ponds.push(pond);
        }

        Ok(ponds)
    }

    fn read_backup_ponds(
        parameters: &Parameters,
        rng: &mut Xoshiro256StarStar,
        sim_path: impl AsRef<Path>,
        time_step: u32
    ) -> anyhow::Result<Vec<Pond>> {
        let sim_path = sim_path.as_ref();
        let mut ponds = vec![];
        for pond_i in 0..parameters.pond.n_ponds {
            log::info!("Reading pond #{pond_i}");
            let cells = IoManager::read_cells(
                IoManager::resolve_cells_path(sim_path, time_step, pond_i),
                IoManager::resolve_genomes_path(sim_path, time_step, pond_i),
            )?;

            let rect = Rect::new(
                (0., 0.).into(),
                (parameters.pond.width as f32, parameters.pond.height as f32).into(),
            );
            let lattice = IoManager::read_lattice(
                IoManager::resolve_lattice_path(sim_path, time_step, pond_i),
                rect.clone().try_into()?,
            )?;

            let mut env = ChemEnvironment::new(
                Environment::new(
                    cells,
                    lattice,
                    NeighbourhoodType::new(parameters.pond.neigh_r),
                    Boundaries::new(BoundaryType::new(rect))?,
                ),
                parameters.cell.max_cells
            );
            let env_ptr: *mut _ = &mut env;
            for pos in env.cell_lattice.iter_positions() {
                // We do this to avoid two lattices in memory
                unsafe { (*env_ptr).update_edges(pos); };
            }

            let mut pond = Self::make_empty_pond(
                parameters, 
                env,
                Self::make_ca(
                    parameters,
                    Some(IoManager::read_clones(IoManager::resolve_clones_path(
                        sim_path,
                        time_step,
                        pond_i
                    ))?)
                ), 
                rng
            );
            pond.time_step = time_step;
            ponds.push(pond);
        }
        Ok(ponds)
    }

    pub fn run_for(&mut self, time_steps: u32) {
        for time_step in self.ponds[0].time_step..=time_steps {
            if self.ponds[0].time_step % self.info_period == 0 {
                self.log_info();
            }
            
            let saved = self.io.write_if_time(
                time_step,
                &self.ponds
            );
            if let Err(e) = saved {
                log::warn!("Failed to save data at time step {time_step} with error `{e}`")
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
    
    pub fn goodbye(&self) {
        log::info!("Finished after {} time steps", self.time_steps);
    }

    fn log_info(&self) {
        log::info!("Time step {}:", self.ponds[0].time_step);
        for (i, pond) in self.ponds.iter().enumerate() {
            let valid = pond.env.cells.n_valid();
            log::info!("\tPond #{i} - {valid} cells");
        }
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