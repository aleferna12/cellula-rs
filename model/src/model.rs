use crate::cell::Cell;
use crate::chem_environment::ChemEnvironment;
use crate::contact_potts::ContactPotts;
use crate::constants::{BoundaryType, NeighbourhoodType};
use crate::evolution::grn::Grn;
use crate::io::io_manager::IoManager;
use crate::io::movie_maker::MovieMaker;
use crate::io::parameters::Parameters;
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
    pub dispersion_period: u32,
    pub info_period: u32,
    time_steps: u32,
    time_step: u32,
}

impl Model {
    pub fn initialise_from_parameters(
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
            dispersion_period: parameters.general.dispersion_period,
            info_period: parameters.io.info_period,
            time_steps: parameters.general.time_steps,
            time_step: 0
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
        let pond = Self::read_backup_pond(
            &parameters,
            &mut rng,
            sim_path,
            time_step
        )?;
        Ok(Self {
            io: Self::setup_io(&parameters, seed)?,
            dispersion_period: parameters.general.dispersion_period,
            info_period: parameters.io.info_period,
            time_steps: parameters.general.time_steps,
            time_step: pond.time_step,
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
            .clones_period(parameters.io.data.clones_period)
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

    fn make_potts(parameters: &Parameters) -> ContactPotts {
        ContactPotts::builder()
            .boltz_t(parameters.potts.boltz_t)
            .size_lambda(parameters.potts.size_lambda)
            .act_max(parameters.potts.act_max)
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
            .build()
    }

    fn determine_seed(seed_param: Option<u64>) -> u64 {
        // TOML doesnt support large u64s so we use a u32 seed
        seed_param.unwrap_or(Xoshiro256StarStar::from_os_rng().next_u32() as u64)
    }

    fn make_empty_pond(
        parameters: &Parameters,
        env: ChemEnvironment,
        ca: ContactPotts,
        rng: &mut Xoshiro256StarStar
    ) -> Pond {
        Pond::builder()
            .env(env)
            .potts(ca)
            .rng(Xoshiro256StarStar::seed_from_u64(rng.next_u64()))
            .update_period(parameters.cell.update_period)
            .cell_target_area(parameters.cell.target_area)
            .cell_search_scaler(parameters.cell.search_radius)
            .division_enabled(parameters.cell.divide)
            .build()
    }

    fn make_new_pond(
        parameters: &Parameters,
        ca: ContactPotts,
        rng: &mut Xoshiro256StarStar
    ) -> anyhow::Result<Pond> {
        let mut env = ChemEnvironment::new(
            Environment::new_empty(
                NeighbourhoodType::new(parameters.pond.neigh_r),
                Boundaries::new(BoundaryType::new(Rect::new(
                    (0., 0.).into(),
                    (parameters.pond.width as f32, parameters.pond.height as f32).into(),
                ))).context("lattice size is too big")?
            ).context("lattice size is too big")?,
            parameters.cell.max_cells,
            parameters.potts.act_max
        );
        if parameters.pond.enclose {
            env.make_border(true, true, true, true);
        }

        log::info!("Making pond");
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
        let lattice = IoManager::read_lattice(
            IoManager::resolve_lattice_path(sim_path, time_step),
            rect.clone().try_into()?,
        )?;

        let mut env = ChemEnvironment::new(
            Environment::new(
                cells,
                lattice,
                NeighbourhoodType::new(parameters.pond.neigh_r),
                Boundaries::new(BoundaryType::new(rect))?,
            ),
            parameters.cell.max_cells,
            parameters.potts.act_max
        );
        env.clones_table = IoManager::read_clones(IoManager::resolve_clones_path(
            sim_path,
            time_step
        ))?;
        let env_ptr: *mut _ = &mut env;
        for pos in env.cell_lattice.iter_positions() {
            // We do this to avoid two lattices in memory
            unsafe { (*env_ptr).update_edges(pos); };
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
            self.time_step,
            &self.pond.env
        );
        if let Err(e) = saved {
            log::warn!("Failed to save data at time step {} with error `{e}`", self.time_step)
        }
        self.pond.step();
        self.time_step = self.pond.time_step;
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