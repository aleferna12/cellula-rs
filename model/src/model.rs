//! Contains logic for creating and running the master [Model] struct.

use crate::cell::{Cell, CellType};
use crate::constants::{BoundaryType, NeighbourhoodType};
use crate::io::io_manager::IoManager;
#[cfg(feature = "movie")]
use crate::io::movie_maker::MovieMaker;
use crate::io::parameters::Parameters;
use crate::my_environment::MyEnvironment;
use crate::my_pond::MyPond;
use crate::my_potts::MyPotts;
use anyhow::Context;
use cellulars_lib::static_adhesion::StaticAdhesion;
use cellulars_lib::environment::Environment;
use cellulars_lib::pond::Pond;
use cellulars_lib::positional::boundaries::{Boundaries, RectConversionError};
use cellulars_lib::positional::rect::Rect;
use cellulars_lib::traits::step::Step;
use rand::{Rng, RngCore, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use std::path::Path;

/// This is the master struct that runs the simulation in a [MyPond] and manages IO through an [IoManager].
pub struct Model {
    /// Pond containing all cells and the model Potts algorithm.
    pub pond: MyPond,
    /// Instance responsible for managing IO for the model.
    pub io: IoManager,
    /// Unique random number generator of this model.
    pub rng: Xoshiro256StarStar,
    /// Period with which information is logged.
    pub info_period: u32,
    time_steps: u32
}

impl Model {
    /// Initialises a brand-new model from some `parameters`.
    pub fn new_from_parameters(
        parameters: Parameters
    ) -> anyhow::Result<Self> {
        log::info!("Initialising model");

        let seed = Self::determine_seed(parameters.general.seed);
        let mut rng = Xoshiro256StarStar::seed_from_u64(seed);
        Ok(Self {
            pond: Self::make_new_pond(
                &parameters,
                Self::make_potts(&parameters),
                &mut rng
            )?,
            io: Self::setup_io(&parameters, seed)?,
            rng,
            info_period: parameters.io.info_period,
            time_steps: parameters.general.time_steps
        })
    }

    /// Initialises the model from a previous state.
    ///
    /// `sim_path` should point to the main folder of a simulation, while `time_step` specifies which files from this
    /// folder will be reloaded.
    pub fn new_from_backup(
        parameters: Parameters,
        sim_path: impl AsRef<Path>,
        time_step: u32
    ) -> anyhow::Result<Self> {
        let sim_path = sim_path.as_ref();
        log::info!("Resuming simulation at {}", sim_path.display());
        log::info!("Starting from time step {time_step}");

        let seed = Self::determine_seed(parameters.general.seed);
        let mut rng = Xoshiro256StarStar::seed_from_u64(seed);
        let pond = Self::read_backup_pond(
            &parameters,
            &mut rng,
            sim_path,
            time_step
        )?;
        Ok(Self {
            io: Self::setup_io(&parameters, seed)?,
            info_period: parameters.io.info_period,
            time_steps: parameters.general.time_steps,
            pond,
            rng,
        })
    }

    fn setup_io(parameters: &Parameters, new_seed: u64) -> anyhow::Result<IoManager> {
        #[cfg(feature = "movie")]
        let movie_maker = if let Some(movie_params) = &parameters.io.movie {
            if movie_params.show {
                match MovieMaker::new(
                    movie_params.width,
                    movie_params.height,
                    movie_params.frame_period
                ) {
                    Ok(mm) => {
                        log::info!("Creating window for real-time movie display");
                        Some(mm)
                    },
                    Err(e) => {
                        log::warn!("Failed to initialise movie window with error `{e}`");
                        None
                    }
                }
            } else {
                None
            }
        } else {
            log::info!("Not displaying movie since movie parameters were omitted");
            None
        };
        #[cfg(not(feature = "movie"))]
        if parameters.io.movie.is_some() {
            log::info!("Not displaying movie since feature flag `movie` was not set");
        }

        let io_builder = IoManager::builder()
            .outdir(parameters.io.outdir.clone().into())
            .image_format(parameters.io.image_format.clone())
            .image_period(parameters.io.image_period)
            .cells_period(parameters.io.data.cells_period)
            .lattice_period(parameters.io.data.lattice_period)
            .plots(parameters.io.plot.clone().try_into()?);
        #[cfg(feature = "movie")]
        let io = io_builder.maybe_movie_maker(movie_maker).build();
        #[cfg(not(feature = "movie"))]
        let io = io_builder.build();

        log::info!("Creating output directories and copy of parameter file");
        if parameters.io.replace_outdir {
            log::info!("Cleaning contents of '{}'", io.outdir.display());
        }
        io.create_directories(parameters.io.replace_outdir)?;
        let mut params_new_seed = parameters.clone();
        params_new_seed.general.seed = new_seed.into();
        io.create_parameters_file(&params_new_seed)?;
        Ok(io)
    }

    fn make_potts(parameters: &Parameters) -> MyPotts {
        MyPotts::builder()
            .boltz_t(parameters.potts.boltz_t)
            .size_lambda(parameters.potts.size_lambda)
            .chemotaxis_mu(parameters.potts.chemotaxis_mu)
            .enable_migration(parameters.cell.migrate)
            .adhesion(
                StaticAdhesion {
                    cell_energy: parameters.potts.adhesion.cell_energy,
                    medium_energy: parameters.potts.adhesion.medium_energy,
                    solid_energy: parameters.potts.adhesion.solid_energy,
                }
            )
            .build()
    }

    fn determine_seed(seed_param: Option<u64>) -> u64 {
        // TOML doesnt support large u64s so we use a u32 seed
        seed_param.unwrap_or(Xoshiro256StarStar::from_os_rng().next_u32() as u64)
    }

    fn make_empty_pond(
        parameters: &Parameters,
        env: MyEnvironment,
        ca: MyPotts,
        rng: &mut Xoshiro256StarStar,
        time_step: u32,
    ) -> MyPond {
        MyPond::new(
            Pond::new(
                env,
                ca,
                Xoshiro256StarStar::seed_from_u64(rng.next_u64()),
                time_step
            ),
            parameters.cell.update_period,
            parameters.cell.divide
        )
    }

    fn make_new_pond(
        parameters: &Parameters,
        ca: MyPotts,
        rng: &mut Xoshiro256StarStar
    ) -> anyhow::Result<MyPond> {
        let env = MyEnvironment::new(
            Environment::new_empty(
                NeighbourhoodType::new(parameters.pond.neigh_r),
                Boundaries::new(BoundaryType::new(Rect::new(
                    (0., 0.).into(),
                    (parameters.pond.width as f32, parameters.pond.height as f32).into(),
                ))).context("lattice size is too big")?
            ).context("lattice size is too big")?,
            parameters.cell.max_cells,
            parameters.cell.search_radius
        );

        log::info!("Making pond");
        let mut pond = Self::make_empty_pond(parameters, env.clone(), ca.clone(), rng, 0);
        for _ in 0..parameters.cell.starting_cells {
            let cell = Cell::new_empty(
                parameters.cell.target_area,
                parameters.cell.div_area,
                if rng.random_bool(0.5) { CellType::Migrating } else { CellType::Dividing }
            );
            pond.pond.env.spawn_cell_random(
                cell,
                parameters.cell.starting_area,
                &mut pond.pond.rng
            );
        }
        log::info!(
                "Created {} out of the {} cells requested",
                pond.pond.env.env.cells.n_valid(),
                parameters.cell.starting_cells
            );

        if parameters.pond.enclose {
            pond.pond.env.make_border(true, true, true, true);
        }
        Ok(pond)
    }

    fn read_backup_pond(
        parameters: &Parameters,
        rng: &mut Xoshiro256StarStar,
        sim_path: impl AsRef<Path>,
        time_step: u32
    ) -> anyhow::Result<MyPond> {
        let sim_path = sim_path.as_ref();

        log::info!("Reading pond");
        let cells = IoManager::read_cells(
            IoManager::resolve_cells_path(sim_path, time_step),
        )?;

        let rect = Rect::new(
            (0., 0.).into(),
            (parameters.pond.width as f32, parameters.pond.height as f32).into(),
        );
        let lattice = IoManager::read_lattice(
            IoManager::resolve_lattice_path(sim_path, time_step),
            rect.to_usize().ok_or(RectConversionError {})?,
        )?;

        let mut env = MyEnvironment::new(
            Environment::new(
                cells,
                lattice,
                NeighbourhoodType::new(parameters.pond.neigh_r),
                Boundaries::new(BoundaryType::new(rect))?,
            ),
            parameters.cell.max_cells,
            parameters.cell.search_radius
        );
        for pos in env.env.cell_lattice.iter_positions() {
            env.env.update_edges(pos);
        }

        let pond = Self::make_empty_pond(
            parameters,
            env,
            Self::make_potts(parameters),
            rng,
            time_step,
        );
        Ok(pond)
    }

    /// Runs the model for the number of time-steps specified when creating the model.
    pub fn run(&mut self) {
        self.run_for(self.time_steps);
    }

    /// Logs some information at end of the simulation.
    pub fn goodbye(&self) {
        log::info!("Finished after {} time steps", self.time_steps);
    }

    fn log_info(&self) {
        log::info!("Time step {}:", self.pond.time_step());
        let valid = self.pond.pond.env.env.cells.n_valid();
        log::info!("\t{valid} cells");
    }
}

impl Step for Model {
    fn step(&mut self) {
        if self.pond.time_step().is_multiple_of(self.info_period) {
            self.log_info();
        }

        let saved = self.io.write_if_time(
            self.pond.time_step(),
            &self.pond.pond.env
        );
        if let Err(e) = saved {
            log::warn!("Failed to save data at time step {} with error `{e}`", self.pond.time_step())
        }
        self.pond.step();
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