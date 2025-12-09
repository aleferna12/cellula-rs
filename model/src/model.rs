use crate::cell::Cell;
use crate::constants::{BoundaryType, NeighbourhoodType};
use crate::evolution::grn::Grn;
use crate::io::io_manager::IoManager;
use crate::io::movie_maker::MovieMaker;
use crate::io::parameters::Parameters;
use crate::my_environment::MyEnvironment;
use crate::my_potts::MyPotts;
use crate::pond::Pond;
use anyhow::Context;
use cellulars_lib::adhesion::StaticAdhesion;
use cellulars_lib::environment::Environment;
use cellulars_lib::positional::boundaries::Boundaries;
use cellulars_lib::positional::rect::Rect;
use cellulars_lib::step::Step;
use rand::distr::{Distribution, Uniform};
use rand::{RngCore, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use std::path::Path;

pub struct Model {
    pub pond: Pond,
    pub io: IoManager,
    pub rng: Xoshiro256StarStar,
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
            .lattices_period(parameters.io.data.lattice_period)
            .plots(parameters.io.plot.clone().try_into()?)
            .maybe_movie_maker(movie_maker)
            .build();

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
            .perimeter_lambda(parameters.potts.perimeter_lambda)
            .act_lambda(parameters.potts.act_lambda)
            .chemotaxis_mu(parameters.potts.chemotaxis_mu)
            .enable_migration(parameters.cell.migrate)
            .adhesion(
                StaticAdhesion {
                    cell_energy: parameters.potts.adhesion.cell_energy,
                    medium_energy: parameters.potts.adhesion.medium_energy,
                    solid_energy: parameters.potts.adhesion.solid_energy,
                }
            )
            .chemotaxis_min(parameters.potts.chemotaxis_min)
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
        rng: &mut Xoshiro256StarStar
    ) -> Pond {
        Pond::builder()
            .env(env)
            .potts(ca)
            .rng(Xoshiro256StarStar::seed_from_u64(rng.next_u64()))
            .update_period(parameters.cell.update_period)
            .cell_target_area(parameters.cell.target_area)
            .enable_division(parameters.cell.divide)
            .enable_cell_updates(parameters.cell.update)
            .season_duration(parameters.pond.season_duration)
            .half_fitness(parameters.pond.half_fitness)
            .reproduction_steps(parameters.pond.reproduction_steps)
            .build()
    }

    fn make_new_pond(
        parameters: &Parameters,
        ca: MyPotts,
        rng: &mut Xoshiro256StarStar
    ) -> anyhow::Result<Pond> {
        let env = MyEnvironment::builder()
            .env(Environment::new_empty(
                NeighbourhoodType::new(parameters.pond.neigh_r),
                Boundaries::new(BoundaryType::new(Rect::new(
                    (0., 0.).into(),
                    (parameters.pond.width as f32, parameters.pond.height as f32).into(),
                ))).context("lattice size is too big")?
            ).context("lattice size is too big")?)
            .max_cells(parameters.cell.max_cells)
            .act_max(parameters.potts.act_max)
            .cell_search_scaler(parameters.cell.search_radius)
            .build();

        log::info!("Making pond");
        let mut pond = Self::make_empty_pond(parameters, env.clone(), ca.clone(), rng);
        // Gradient has to be initialised before cells are added so that chem center pos are tracked
        pond.env.make_next_chem_gradient(&mut pond.rng);
        for _ in 0..parameters.cell.starting_cells {
            let cell = Cell::new_empty(
                parameters.cell.target_area,
                parameters.cell.target_perimeter,
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

        if parameters.pond.enclose {
            pond.env.make_border(true, true, true, true);
        }
        Ok(pond)
    }

    fn read_backup_pond(
        parameters: &Parameters,
        rng: &mut Xoshiro256StarStar,
        sim_path: impl AsRef<Path>,
        time_step: u32
    ) -> anyhow::Result<Pond> {
        let sim_path = sim_path.as_ref();

        log::info!("Reading pond");
        let cells = IoManager::read_cells(
            IoManager::resolve_cells_path(sim_path, time_step),
            IoManager::resolve_genomes_path(sim_path, time_step),
        )?;

        let rect = Rect::new(
            (0., 0.).into(),
            (parameters.pond.width as f32, parameters.pond.height as f32).into(),
        );
        let rect_usize: Rect<usize> = rect.clone().try_into()?;
        let cell_lattice = IoManager::read_cell_lattice(
            IoManager::resolve_cell_lattice_path(sim_path, time_step),
            rect_usize.clone(),
        )?;
        let chem_lattice = IoManager::read_lattice_u32(
            IoManager::resolve_chem_lattice_path(sim_path, time_step),
            rect_usize.clone(),
        )?;
        let act_lattice = IoManager::read_lattice_u32(
            IoManager::resolve_act_lattice_path(sim_path, time_step),
            rect_usize,
        )?;

        let mut env = MyEnvironment::new_from_backup()
            .env(Environment::new(
                cells,
                cell_lattice,
                NeighbourhoodType::new(parameters.pond.neigh_r),
                Boundaries::new(BoundaryType::new(rect))?,
            ))
            .chem_lattice(chem_lattice)
            .act_lattice(act_lattice)
            .max_cells(parameters.cell.max_cells)
            .act_max(parameters.potts.act_max)
            .cell_search_scaler(parameters.cell.search_radius)
            .call();
        for pos in env.cell_lattice.iter_positions() {
            // We do this to avoid two lattices in memory
            env.update_edges(pos);
        }

        let mut pond = Self::make_empty_pond(
            parameters,
            env,
            Self::make_potts(parameters),
            rng
        );
        pond.time_step = time_step;
        Ok(pond)
    }

    pub fn run(&mut self) {
        self.run_for(self.time_steps);
    }

    pub fn goodbye(&self) {
        log::info!("Finished after {} time steps", self.time_steps);
    }

    fn log_info(&self) {
        log::info!("Time step {}:", self.pond.time_step);
        let valid = self.pond.env.cells.n_valid();
        log::info!("\t{valid} cells");
    }
}

impl Step for Model {
    fn step(&mut self) {
        if self.pond.time_step % self.info_period == 0 {
            self.log_info();
        }

        let saved = self.io.write_if_time(
            self.pond.time_step,
            &self.pond.env
        );
        if let Err(e) = saved {
            log::warn!("Failed to save data at time step {} with error `{e}`", self.pond.time_step)
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