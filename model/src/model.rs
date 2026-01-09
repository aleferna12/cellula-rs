use crate::bit_adhesion::BitAdhesion;
use crate::cell::Cell;
use crate::constants::{BoundaryType, NeighbourhoodType};
use crate::evolution::bit_genome::BitGenome;
use crate::io::io_manager::IoManager;
use crate::io::movie_maker::MovieMaker;
use crate::io::parameters::Parameters;
use crate::my_environment::MyEnvironment;
use crate::my_potts::MyPotts;
use crate::pond::Pond;
use anyhow::{anyhow, Context};
use cellulars_lib::adhesion::StaticAdhesion;
use cellulars_lib::basic_cell::{Alive, Cellular};
use cellulars_lib::constants::CellIndex;
use cellulars_lib::environment::Environment;
use cellulars_lib::habitable::Habitable;
use cellulars_lib::positional::boundaries::Boundaries;
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::positional::rect::Rect;
use cellulars_lib::step::Step;
use polars::polars_utils::itertools::Itertools;
use rand::{Rng, RngCore, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use std::collections::HashMap;
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
        parameters: Parameters,
        maybe_template_path: Option<String>,
    ) -> anyhow::Result<Self> {
        log::info!("Initialising model");

        let seed = Self::determine_seed(parameters.general.seed);
        let mut rng = Xoshiro256StarStar::seed_from_u64(seed);
        Ok(Self {
            pond: Self::make_new_pond(
                &parameters,
                &mut rng,
                maybe_template_path
            )?,
            io: Self::setup_io(&parameters, seed)?,
            rng,
            info_period: parameters.io.info_period,
            time_steps: parameters.general.time_steps
        })
    }

    /// Makes a new model from a layout file.
    ///
    /// Layout specifications are documented in the CLI.
    pub fn new_from_layout(
        parameters: Parameters,
        layout_path: impl AsRef<Path>,
        maybe_templates_path: Option<String>
    ) -> anyhow::Result<Self> {
        let layout_path = layout_path.as_ref();
        log::info!("Initializing model with layout \"{}\"", layout_path.display());

        let seed = Self::determine_seed(parameters.general.seed);
        let mut rng = Xoshiro256StarStar::seed_from_u64(seed);
        let pond = Self::read_layout_pond(&parameters, layout_path, &mut rng, maybe_templates_path)?;
        Ok(Self {
            pond,
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
            .adhesion(BitAdhesion { 
                static_adhesion: StaticAdhesion {
                    cell_energy: parameters.potts.adhesion.cell_energy,
                    medium_energy: parameters.potts.adhesion.medium_energy,
                    solid_energy: parameters.potts.adhesion.solid_energy,
                },
                gene_energy: parameters.potts.adhesion.gene_energy,
            })
            .chemotaxis_min(parameters.potts.chemotaxis_min)
            .build()
    }

    fn determine_seed(seed_param: Option<u64>) -> u64 {
        // TOML doesnt support large u64s so we use a u32 seed
        seed_param.unwrap_or(Xoshiro256StarStar::from_os_rng().next_u32() as u64)
    }

    fn make_env(parameters: &Parameters) -> anyhow::Result<MyEnvironment> {
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
        Ok(env)
    }

    fn make_empty_pond(parameters: &Parameters, rng: &mut Xoshiro256StarStar) -> anyhow::Result<Pond> {
        let pond = Pond::builder()
            .env(Self::make_env(parameters)?)
            .potts(Self::make_potts(parameters))
            .rng(Xoshiro256StarStar::seed_from_u64(rng.next_u64()))
            .time_step(0)
            .enable_division(parameters.cell.divide)
            .season_duration(parameters.pond.season_duration)
            .reproduction_steps(parameters.pond.reproduction_steps)
            .half_fitness(parameters.pond.half_fitness)
            .cell_target_area(parameters.cell.target_area)
            .build();
        Ok(pond)
    }

    fn templates_path_to_box(
        maybe_templates_path: Option<String>,
        mut_rate: f32,
        genome_len: u8
    ) -> anyhow::Result<Option<Box<[Cell]>>> {
        maybe_templates_path.map(|path| {
            // This is required to obtain a clonable iterator that we can cycle over
            // TODO!: why do cells not save mut rate and genome len?
            let templates_cells = IoManager::read_cells(path, mut_rate, genome_len)?
                .into_iter()
                .map(|rel_cell| rel_cell.cell)
                .collect::<Box<[_]>>();
            anyhow::Ok(templates_cells)
        }).transpose()
    }

    fn read_layout_pond(
        parameters: &Parameters,
        layout_path: impl AsRef<Path>,
        rng: &mut Xoshiro256StarStar,
        maybe_templates_path: Option<String>
    ) -> anyhow::Result<Pond> {
        let layout_path = layout_path.as_ref();

        let layout = IoManager::read_layout(
            layout_path,
            parameters.pond.width,
            parameters.pond.height
        )?;

        // Using floor bc thats what we use in spawn_cell_random
        let cell_side = parameters.cell.starting_area.isqrt() as usize;
        let mut solid_positions = vec![];
        // luma values -> (grid_indexes -> positions)
        let mut luma_cell_positions = HashMap::new();
        for j in 0..parameters.pond.height {
            for i in 0..parameters.pond.width {
                let luma = layout[(i as u32, j as u32)].0[0];
                if luma == 255 {
                    continue;
                }

                let pos = Pos::new(i, j);
                if luma == 0 {
                    solid_positions.push(pos);
                    continue;
                }

                let grid_index = Pos::new(
                    i / cell_side,
                    j / cell_side
                ).col_major(parameters.pond.height) as CellIndex;
                let cell_positions = luma_cell_positions
                    .entry(luma)
                    .or_insert(HashMap::new());
                let positions = cell_positions
                    .entry(grid_index)
                    .or_insert_with(Vec::new);
                positions.push(pos);
            }
        }
        log::info!("Found {} different groups in the layout", luma_cell_positions.len());

        let mut sorted_luma = luma_cell_positions.keys().copied().collect_vec();
        sorted_luma.sort();

        let mut pond = Self::make_empty_pond(parameters, rng)?;
        // Gradient has to be initialised before cells are added so that chem center pos are tracked
        pond.env.make_next_chem_gradient(&mut pond.rng);
        let maybe_templates_box = Self::templates_path_to_box(
            maybe_templates_path,
            parameters.cell.genome.mutation_rate,
            parameters.cell.genome.length
        )?;
        for (group_index, luma) in sorted_luma.into_iter().enumerate() {
            let cell_positions = luma_cell_positions
                .remove(&luma)
                .expect("missing luma key");
            for positions in cell_positions.values() {
                if !pond.env.can_add_cell() {
                    continue;
                }
                let cell = match &maybe_templates_box {
                    None => Self::empty_cell_from_parameters(parameters, rng)?,
                    Some(templates_box) => templates_box
                        .get(group_index)
                        .ok_or(anyhow::anyhow!("there were more groups in the layout than in the template"))?
                        .clone()
                };
                pond.env.spawn_cell(cell.birth(), positions.iter().copied());
            }
        }
        pond.env.spawn_solid(solid_positions.into_iter());
        if parameters.pond.enclose {
            pond.env.make_border(true, true, true, true);
        }
        Ok(pond)
    }

    fn empty_cell_from_parameters(
        parameters: &Parameters,
        rng: &mut impl Rng
    ) -> anyhow::Result<Cell> {
        Ok(Cell::new_empty(
            parameters.cell.target_area,
            parameters.cell.target_perimeter,
            BitGenome::new_random(
                parameters.cell.genome.mutation_rate,
                parameters.cell.genome.length,
                rng,
            ).ok_or(anyhow!("invalid `parameters.cell.genome.length`"))?
        ))
    }

    fn make_new_pond(
        parameters: &Parameters,
        rng: &mut Xoshiro256StarStar,
        maybe_templates_path: Option<String>
    ) -> anyhow::Result<Pond> {
        log::info!("Making pond");
        let mut pond = Self::make_empty_pond(parameters, rng)?;
        // Gradient has to be initialised before cells are added so that chem center pos are tracked
        pond.env.make_next_chem_gradient(&mut pond.rng);

        // Obtains an iterator over cell templates if a templates_path is present
        let maybe_templates_box = Self::templates_path_to_box(
            maybe_templates_path,
            parameters.cell.genome.mutation_rate,
            parameters.cell.genome.length
        )?;
        let mut maybe_templates_it = maybe_templates_box.map(|templates_box| templates_box.into_iter().cycle());
        let mut spawn_attempts = 0;
        while pond.env.cells.n_valid() < parameters.cell.starting_cells && pond.env.can_add_cell() {
            let cell = match &mut maybe_templates_it {
                None => Self::empty_cell_from_parameters(parameters, rng)?,
                Some(templates_it) => templates_it
                    .next()
                    .ok_or(anyhow::anyhow!("failed to obtain cell from template iterator"))?
            };
            let cell_area = if cell.area() == 0 { parameters.cell.starting_area } else { cell.area() };
            pond.env.spawn_cell_random(
                cell.birth(),
                cell_area,
                &mut pond.rng
            );
            spawn_attempts += 1;

            if spawn_attempts == parameters.cell.starting_cells * 2 {
                log::warn!("Parameters have led to high cell density and difficulties placing cells in the simulation");
                log::warn!("Consider decreasing `cell.starting_cells` or increasing the pond area");
            } else if spawn_attempts > parameters.cell.starting_cells * 20 {
                log::error!(
                    "Only {} cells were initialized out of {} cells requested",
                    pond.env.cells.n_valid(),
                    parameters.cell.starting_cells);
                break;
            }
        }
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
            parameters.cell.genome.mutation_rate,
            parameters.cell.genome.length
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

        let mut pond = Pond::builder()
            .env(env)
            .potts(Self::make_potts(parameters))
            .rng(Xoshiro256StarStar::seed_from_u64(rng.next_u64()))
            .time_step(time_step)
            .half_fitness(parameters.pond.half_fitness)
            .season_duration(parameters.pond.season_duration)
            .cell_target_area(parameters.cell.target_area)
            .enable_division(parameters.cell.divide)
            .reproduction_steps(parameters.pond.reproduction_steps)
            .build();
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
        if self.pond.time_step.is_multiple_of(self.info_period) {
            self.log_info();
        }

        let saved = self.io.write_if_time(
            self.pond.time_step,
            &mut self.pond.env
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