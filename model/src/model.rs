//! Contains logic for creating and running the master [`Model`] struct.

use crate::constants::{BoundaryType, KinectNeighbourhoodType, NeighbourhoodType};
use crate::io::io_manager::IoManager;
use crate::io::kinect_listener::KinectListener;
#[cfg(feature = "movie")]
use crate::io::movie_maker::MovieMaker;
use crate::io::parameters::Parameters;
use crate::my_cell::{CellType, MyCell};
use crate::my_environment::MyEnvironment;
use crate::my_pond::MyPond;
use crate::pairwise_adhesion::PairwiseAdhesion;
use crate::physics::Physics;
use crate::potts::Potts;
use cellulars_lib::base::environment::Environment;
use cellulars_lib::base::pond::Pond;
use cellulars_lib::constants::FloatType;
use cellulars_lib::perimeter_constraint::PerimeterConstraint;
use cellulars_lib::positional::boundaries::Boundaries;
use cellulars_lib::positional::pos::CastCoords;
use cellulars_lib::positional::rect::Rect;
use cellulars_lib::prelude::{Alive, CellIndex, Cellular, Habitable, Pos, Spin};
use cellulars_lib::traits::cellular::EmptyCell;
use cellulars_lib::traits::step::Step;
use image::imageops::{FilterType, crop_imm, flip_vertical_in_place, overlay, resize};
use image::{ImageReader, Rgba, RgbaImage};
use minifb::Key;
use polars::polars_utils::itertools::Itertools;
use rand::{Rng, RngCore, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::path::Path;
use std::thread::sleep;
use std::time::{Duration, Instant};

/// This is the master struct that runs the simulation in a [`MyPond`] and manages IO through an [`IoManager`].
pub struct Model {
    /// Pond containing all cells and the model Potts algorithm.
    pub my_pond: MyPond,
    pub physics: Physics,
    /// Instance responsible for managing IO for the model.
    pub io: IoManager,
    /// Unique random number generator of this model.
    pub rng: Xoshiro256StarStar,
    /// Period with which information is logged.
    pub info_period: u32,
    time_steps: u32,
    display_mode: DisplayMode,
    last_mode_change: Instant,
}

impl Model {
    /// Initialises a brand-new model from some `parameters`.
    pub fn new_from_parameters(
        parameters: Parameters,
        maybe_templates_path: Option<String>,
    ) -> anyhow::Result<Self> {
        log::info!("Initializing model");

        let seed = Self::determine_seed(parameters.general.seed);
        let mut rng = Xoshiro256StarStar::seed_from_u64(seed);
        Ok(Self {
            my_pond: Self::make_new_pond(
                &parameters,
                &mut rng,
                maybe_templates_path
            )?,
            io: Self::setup_io(&parameters, seed)?,
            physics: Physics::new(20., 5, &mut rng),
            rng,
            info_period: parameters.io.info_period,
            time_steps: parameters.general.time_steps,
            display_mode: DisplayMode::Macro,
            last_mode_change: Instant::now(),
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
            my_pond: pond,
            io: Self::setup_io(&parameters, seed)?,
            physics: Physics::new(20., 5, &mut rng),
            rng,
            info_period: parameters.io.info_period,
            time_steps: parameters.general.time_steps,
            display_mode: DisplayMode::Macro,
            last_mode_change: Instant::now(),
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
        log::info!("Resuming simulation at \"{}\"", sim_path.display());
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
            physics: Physics::new(20., 5, &mut rng),
            io: Self::setup_io(&parameters, seed)?,
            info_period: parameters.io.info_period,
            time_steps: parameters.general.time_steps,
            my_pond: pond,
            rng,
            display_mode: DisplayMode::Macro,
            last_mode_change: Instant::now(),
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

        let kinect_listener = match &parameters.io.kinect {
            None => {
                log::info!("Not displaying movie since movie parameters were omitted");
                None
            }
            Some(params) => {
                if !params.allow {
                    None
                } else {
                    KinectListener::new(
                        params.min_depth,
                        params.max_depth,
                        params.frame_period,
                        KinectNeighbourhoodType::new(1)
                    )
                }
            }
        };

        let io_builder = IoManager::builder()
            .outdir(parameters.io.outdir.clone().into())
            .image_format(parameters.io.image_format.clone())
            .image_period(parameters.io.image_period)
            .cells_period(parameters.io.data.cells_period)
            .lattice_period(parameters.io.data.lattice_period)
            .plots(parameters.io.plot.clone().try_into()?)
            .micro_bg(RgbaImage::new(512, 424))
            .macro_bg(RgbaImage::new(512, 424))
            .kinect_img(RgbaImage::new(512, 424))
            .maybe_kinect_listener(kinect_listener);
        #[cfg(feature = "movie")]
        let mut io = io_builder.maybe_movie_maker(movie_maker).build();
        #[cfg(not(feature = "movie"))]
        let mut io = io_builder.build();
        io.load_bgs();

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

    fn make_potts(parameters: &Parameters, rng: &mut impl Rng) -> Potts {
        let mut adh = PairwiseAdhesion::new(
            parameters.potts.adhesion.medium_energy,
            parameters.potts.adhesion.solid_energy,
            parameters.cell.max_cells
        );
        adh.randomize_cell_energies(
            parameters.potts.adhesion.cell_energy * 0.25,
            parameters.potts.adhesion.cell_energy * 16.,
            rng
        );
        Potts::builder()
            .boltz_t(parameters.potts.boltz_t)
            .size_lambda(parameters.potts.size_lambda)
            .chemotaxis_mu(parameters.potts.chemotaxis_mu)
            .persistence_mu(parameters.potts.persistence_mu)
            .enable_migration(parameters.cell.migrate)
            .perimeter(PerimeterConstraint { lambda: parameters.potts.perimeter_lambda })
            .adhesion(adh)
            .build()
    }

    fn make_env(parameters: &Parameters) -> MyEnvironment {
        MyEnvironment::new(
            Environment::new_empty(
                NeighbourhoodType::new(parameters.pond.neigh_r),
                Boundaries::new(BoundaryType::new(Rect::new(
                    (0., 0.).into(),
                    (parameters.pond.width as FloatType, parameters.pond.height as FloatType).into(),
                )))
            ),
            parameters.cell.max_cells,
            parameters.cell.search_radius,
            Self::read_target_image(&parameters.io.target_png),
            Pos::new(
                parameters.pond.target_x,
                parameters.pond.target_y
            )
        )
    }

    fn read_target_image(path: &str) -> RgbaImage {
        let mut img = ImageReader::open(path).unwrap_or_else(|_| panic!("failed to open {}",
            path)).decode().expect("failed to decode target image").into_rgba8();
        flip_vertical_in_place(&mut img);
        img
    }

    fn determine_seed(seed_param: Option<u64>) -> u64 {
        // TOML doesnt support large u64s so we use a u32 seed
        seed_param.unwrap_or(Xoshiro256StarStar::from_os_rng().next_u32() as u64)
    }

    fn make_empty_pond(parameters: &Parameters, rng: &mut Xoshiro256StarStar) -> MyPond {
        MyPond::new(
            Pond::new(
                Self::make_env(parameters),
                Self::make_potts(parameters, rng),
                Xoshiro256StarStar::seed_from_u64(rng.next_u64()),
                0
            ),
            parameters.cell.update_period,
            parameters.cell.divide,
            parameters.pond.target_move_period
        )
    }

    fn templates_path_to_box(
        maybe_templates_path: Option<String>
    ) -> anyhow::Result<Option<Box<[MyCell]>>> {
        maybe_templates_path.map(|path| {
            // This is required to obtain a clonable iterator that we can cycle over
            let templates_cells = IoManager::read_cells(path)?
                .into_iter()
                .map(|rel_cell| rel_cell.cell)
                .collect::<Box<[_]>>();
            anyhow::Ok(templates_cells)
        }).transpose()
    }

    fn empty_cell_from_parameters(parameters: &Parameters, rng: &mut impl Rng) -> EmptyCell<MyCell> {
        let area = parameters.cell.target_area as FloatType * rng.random_range(0.25..2.);
        // Estimated from the biophys paper
        let per = 3. * 2. * PI * (area / PI).sqrt() * rng.random_range(1.0..1.25);
        let pers = parameters.cell.persistence_duration as FloatType * rng.random_range(0.5..2.);
        MyCell::new_empty(
            area as u32,
            per as u32,
            parameters.cell.div_area,
            pers as u32,
            CellType::Migrating
        )
    }

    fn make_new_pond(
        parameters: &Parameters,
        rng: &mut Xoshiro256StarStar,
        maybe_templates_path: Option<String>,
    ) -> anyhow::Result<MyPond> {
        log::info!("Making pond");
        let mut pond = Self::make_empty_pond(parameters, rng);
        pond.env_mut().update_chem_gradient();

        // Obtains an iterator over cell templates if a templates_path is present
        let maybe_templates_box = Self::templates_path_to_box(maybe_templates_path)?;
        let mut maybe_templates_it = maybe_templates_box.map(|templates_box| templates_box.into_iter().cycle());
        let mut spawn_attempts = 0;
        while pond.env().env.cells.n_non_empty() < parameters.cell.starting_cells {
            let cell = match &mut maybe_templates_it {
                None => Self::empty_cell_from_parameters(parameters, rng).into_cell(),
                Some(templates_it) => templates_it
                    .next()
                    .ok_or(anyhow::anyhow!("failed to obtain cell from template iterator"))?
            };
            let cell_area = if cell.area() == 0 {
                parameters.cell.starting_area
            } else {
                cell.area()
            };
            pond.pond.env.spawn_cell_random(
                cell.birth(),
                cell_area,
                &mut pond.pond.rng
            );
            spawn_attempts += 1;

            if spawn_attempts == parameters.cell.starting_cells * 2 {
                log::warn!("Parameters have led to high cell density and difficulties placing cells in the simulation");
                log::warn!("Consider decreasing `cell.starting_cells` or increasing the pond area");
            } else if spawn_attempts > parameters.cell.starting_cells * 20 {
                log::error!(
                    "Only {} cells were initialized out of {} cells requested",
                    pond.env().env.cells.n_non_empty(),
                    parameters.cell.starting_cells);
                break;
            }
        }
        if parameters.pond.enclose {
            pond.env_mut().make_border(true, true, true, true);
        }
        let Pond { env, rng, .. } = &mut pond.pond;
        env.init_migration(rng);
        Ok(pond)
    }

    fn read_layout_pond(
        parameters: &Parameters,
        layout_path: impl AsRef<Path>,
        rng: &mut Xoshiro256StarStar,
        maybe_templates_path: Option<String>
    ) -> anyhow::Result<MyPond> {
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

        let mut sorted_luma = luma_cell_positions.keys().copied().collect_vec();
        sorted_luma.sort();

        let mut not_spawned = 0;
        let mut pond = Self::make_empty_pond(parameters, rng);
        let maybe_templates_box = Self::templates_path_to_box(maybe_templates_path)?;
        for (group_index, luma) in sorted_luma.into_iter().enumerate() {
            let cell_positions = luma_cell_positions
                .remove(&luma)
                .expect("missing luma key");
            for positions in cell_positions.values() {
                if pond.env().env.cells.n_non_empty() >= parameters.cell.starting_cells {
                    not_spawned += 1;
                    continue;
                }
                let cell = match &maybe_templates_box {
                    None => Self::empty_cell_from_parameters(parameters, rng).into_cell(),
                    Some(templates_box) => templates_box
                        .get(group_index)
                        .ok_or(anyhow::anyhow!("there were more groups in the layout than in the template"))?
                        .clone()
                };
                pond.env_mut().spawn_cell(cell.birth(), positions.iter().copied());
            }
        }
        log::warn!("Number of cells spawned was limited by `cell.starting-cells` parameter \
                    (layout had space for {not_spawned} additional cells)");

        pond.env_mut().spawn_solid(solid_positions.into_iter());
        if parameters.pond.enclose {
            pond.env_mut().make_border(true, true, true, true);
        }
        let Pond { env, rng, .. } = &mut pond.pond;
        env.init_migration(rng);
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
            (parameters.pond.width as FloatType, parameters.pond.height as FloatType).into(),
        );
        let lattice = IoManager::read_lattice(
            IoManager::resolve_lattice_path(sim_path, time_step),
            rect.cast_coords(),
        )?;

        let mut env = MyEnvironment::new(
            Environment::new(
                cells,
                lattice,
                NeighbourhoodType::new(parameters.pond.neigh_r),
                Boundaries::new(BoundaryType::new(rect)),
            ),
            parameters.cell.max_cells,
            parameters.cell.search_radius,
            Self::read_target_image(&parameters.io.target_png),
            Pos::new(
                parameters.pond.target_x,
                parameters.pond.target_y
            )
        );
        for pos in env.env.cell_lattice.iter_positions() {
            env.env.update_edges(pos);
        }
        // The cells and the lattice are read from separate files, so the perimeters the cells were tracking
        // when the back-up was written are re-measured instead of stored
        env.restore_perimeters();
        env.init_migration(rng);

        let pond = MyPond::new(
            Pond::new(
                env,
                Self::make_potts(parameters, rng),
                Xoshiro256StarStar::seed_from_u64(rng.next_u64()),
                time_step
            ),
            parameters.cell.update_period,
            parameters.cell.divide,
            parameters.pond.target_move_period
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
        log::info!("Time step {}:", self.my_pond.time_step());
        let non_empty = self.my_pond.env().env.cells.n_non_empty();
        log::info!("\t{non_empty} cells");
    }

    fn cpm_step(&mut self, write: bool) {
        if self.my_pond.time_step().is_multiple_of(self.info_period) {
            self.log_info();
        }

        if write {
            let saved = self.io.write_if_time(
                self.my_pond.time_step(),
                self.my_pond.env()
            );
            if let Err(e) = saved {
                log::warn!("Failed to save data at time step {} with error `{e}`", self.my_pond.time_step())
            }
        }

        // Draw silhouette
        if let Some(kinect) = &mut self.io.kinect_listener
            && self.my_pond.time_step().is_multiple_of(kinect.frame_period) {
            // Clear solids
            for spin in self.my_pond.env_mut().env.cell_lattice.iter_values_mut() {
                if matches!(spin, Spin::Solid) {
                    *spin = Spin::Medium;
                }
            }

            self.io.kinect_img = kinect.draw_silhouette(self.my_pond.env_mut())
                .expect("failed to draw silhouette from kinect");
            // self.my_pond.env_mut().draw_solid_target();
        }

        self.my_pond.step();
    }

    fn physics_step(&mut self, write: bool) {
        if write && self.my_pond.time_step().is_multiple_of(100) {
            let img = self.physics_image();
            if let Some(mm) = &mut self.io.movie_maker {
                mm.update(&img).unwrap();
            }
        }

        self.physics.step(&self.my_pond.pond.env.env.cell_lattice);
    }

    fn physics_image(&mut self) -> RgbaImage {
        let mut img = self.io.macro_bg.clone();
        for ball in &self.physics.balls {
            for pos in ball.rectangle().iter_positions() {
                if pos.y >= self.my_pond.pond.env.env.height() as u32 {
                    continue;
                }
                let dist = (pos.x as FloatType - ball.center.x).hypot(pos.y as FloatType - ball.center.y);
                if dist > ball.radius {
                    continue;
                }
                let color = if dist < 17. { Rgba([128, 50, 50, 0]) } else { Rgba([0, 0, 0, 0]) };
                img.put_pixel(pos.x, pos.y, color);
            }
        }
        overlay(&mut img, &self.io.kinect_img, 0, 0);

        flip_vertical_in_place(&mut img);
        img
    }

    fn update_mode(&mut self) {
        let now = Instant::now();
        if let Some(mm) = &self.io.movie_maker
            && mm.window.is_open()
            && mm.window.is_key_released(Key::Space)
            && (now - self.last_mode_change) > Duration::from_secs(5) {
            if self.display_mode == DisplayMode::Macro {
                self.zoom_in();
                self.display_mode = DisplayMode::Micro;
            } else {
                self.display_mode = DisplayMode::Macro;
            }
            self.last_mode_change = now;
        }
    }

    fn zoom_in(&mut self) {
        let original = self.physics_image();
        let cropx = (original.width() as FloatType * 0.01) as u32;
        let cropy = (original.height() as FloatType * 0.01) as u32;
        for i in 1..50 {
            let mut img = crop_imm(
                &original,
                cropx * i,
                cropy * i,
                original.width() - 2 * cropx * i,
                original.height() - 2 * cropy * i
            ).to_image();
            img = resize(&img, 512, 424, FilterType::Gaussian);
            if let Some(mm) = &mut self.io.movie_maker {
                mm.update(&img).unwrap();
                sleep(Duration::from_millis(40));
            }
        }
    }
}

impl Step for Model {
    fn step(&mut self) {
        self.update_mode();
        // We run both steps all the time and hope the computer is fast enough lol
        // This prevents cells being erased while we interact with the ball (the lazy way to do it)
        self.physics_step(self.display_mode == DisplayMode::Macro);
        self.cpm_step(self.display_mode == DisplayMode::Micro);
    }
}

#[derive(Debug, PartialEq)]
pub enum DisplayMode {
    Micro,
    Macro
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
