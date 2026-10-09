//! Contains logic related to [`IoManager`].

// TODO!: file names should have leading zeros (account for max size of time_steps)

use crate::my_cell::{MyCell, CellType};
use crate::my_environment::MyEnvironment;
#[cfg(feature = "movie")]
use crate::io::movie_maker::MovieMaker;
use crate::io::parameters::Parameters;
use crate::io::plot::Plot;
use anyhow::{bail, Context};
use bon::Builder;
use cellulars_lib::base::cell::Cell;
use cellulars_lib::cell_container::{CellContainer, RelCell};
use cellulars_lib::constants::CellIndex;
use cellulars_lib::lattice::Lattice;
use cellulars_lib::positional::com::Com;
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::positional::rect::Rect;
use cellulars_lib::spin::Spin;
use cellulars_lib::traits::cellular::Cellular;
use cellulars_lib::traits::track_perimeter::TrackPerimeter;
use image::imageops::{flip_vertical_in_place, FilterType, flip_vertical, overlay, resize, flip_horizontal_in_place};
use image::{open, ColorType, GrayImage, ImageReader, RgbaImage};
use num_traits::NumCast;
use polars::frame::row::Row;
use polars::polars_utils::float::IsFloat;
use polars::prelude::*;
use std::collections::HashSet;
use std::f64::consts::PI;
use std::io;
use std::path::{Path, PathBuf};
use crate::io::kinect_listener::KinectListener;

static IMAGES_PATH: &str = "images";
static CELLS_PATH: &str = "cells";
static LATTICES_PATH: &str = "lattices";
static CONFIG_COPY_PATH: &str = "config.toml";
const PAD_FILE_LEN: usize = {
    let mut n = u32::MAX;
    let mut digits = 0;
    while n > 0 {
        digits += 1;
        n /= 10;
    }
    digits
};

/// Manages all io operations, including saving and loading data and displaying the simulation movie.
#[derive(Builder)]
pub struct IoManager {
    /// Path to directory where data and images of the simulation are saved.
    pub outdir: PathBuf,
    /// Image format with which to save simulation images.
    pub image_format: String,
    /// Used to update the simulation video when it's time.
    #[cfg(feature = "movie")]
    pub movie_maker: Option<MovieMaker>,
    pub kinect_listener: Option<KinectListener>,
    pub micro_bg: RgbaImage,
    pub macro_bg: RgbaImage,
    pub ball_img: RgbaImage,
    pub kinect_img: RgbaImage,
    pub eyes_img: RgbaImage,
    plots: Box<[Box<dyn Plot>]>,
    image_period: u32,
    cells_period: u32,
    lattice_period: u32
}

impl IoManager {
    /// Create the main simulation folder and all subdirectories.
    ///
    /// Fails if `replace_outdir` is `false` and the main simulation folder already exists.
    pub fn create_directories(&self, replace_outdir: bool) -> io::Result<()> {
        let outdir_exists = self.outdir.try_exists()?;
        if outdir_exists {
            if replace_outdir {
                std::fs::remove_dir_all(&self.outdir)?;
            } else {
                return Err(io::Error::new(
                    io::ErrorKind::AlreadyExists,
                    "`outdir` already exists and `replace_outdir` is `false`"
                ));
            }
        }
        std::fs::create_dir_all(&self.outdir)?;
        std::fs::create_dir(self.outdir.join(IMAGES_PATH))?;
        std::fs::create_dir(self.outdir.join(CELLS_PATH))?;
        std::fs::create_dir(self.outdir.join(LATTICES_PATH))
    }

    /// Creates a parameter file at \[`IoManager::outdir`]/config.toml\".
    pub fn create_parameters_file(&self, parameters: &Parameters) -> anyhow::Result<()> {
        let params_copy = self.outdir.join(CONFIG_COPY_PATH);
        std::fs::write(
            params_copy,
            format!(
                "{}\n{}",
                "# This is a copy of the parameters used in the simulation",
                toml::to_string(parameters)?
            )
        )?;
        Ok(())
    }

    pub fn load_imgs(&mut self) {
        self.micro_bg = flip_vertical(&open("./bg_micro.png").unwrap().into_rgba8());
        self.macro_bg = flip_vertical(&open("./bg_macro.png").unwrap().into_rgba8());
        self.ball_img = flip_vertical(&open("./ball.png").unwrap().into_rgba8());
        self.eyes_img = flip_vertical(&open("./eyes.png").unwrap().into_rgba8());
    }

    fn make_cells_from_data(celldf: DataFrame) -> anyhow::Result<CellContainer<MyCell>> {
        let last_index = celldf
            .column("index")?
            .u32()?
            .max()
            .ok_or(anyhow::anyhow!("null `index` column"))?;
        let mut cells = CellContainer::new();
        // We need this to call replace on cells later
        for _ in 0..=last_index {
            cells.push(MyCell::new_empty(0, 0, 0, 0, CellType::Migrating, 0));
        }

        for row_i in 0..celldf.height() {
            let row = celldf.get_row(row_i)?;
            let cell = Cell::new_ready(
                Self::get_col_num(&row, "area", &celldf)?,
                Self::get_col_num(&row, "target_area", &celldf)?,
                Pos::new(
                    Self::get_col_num(&row, "center_x", &celldf)?,
                    Self::get_col_num(&row, "center_y", &celldf)?,
                )
            );
            cells.replace(RelCell {
                index: Self::get_col_num(&row, "index", &celldf)?,
                cell: MyCell::builder()
                    .cell(cell)
                    .divide_area(Self::get_col_num(&row, "divide_area", &celldf)?)
                    .newborn_target_area(Self::get_col_num(&row, "newborn_target_area", &celldf)?)
                    .chem_com(Com {
                        pos: Pos::new(
                            Self::get_col_num(&row, "chem_center_x", &celldf)?,
                            Self::get_col_num(&row, "chem_center_y", &celldf)?,
                        ),
                        mass: Self::get_col_num(&row, "chem_mass", &celldf)?
                    })
                    // The perimeter of the cell is not read back, since it is re-measured from the
                    // lattice by `MyEnvironment::restore_perimeters()`
                    .target_perimeter(Self::get_col_num_or(&row, "target_perimeter", &celldf, 0)?)
                    .persistence_duration(Self::get_col_num_or(&row, "persistence_duration", &celldf, 0)?)
                    .persistence_time(Self::get_col_num_or(&row, "persistence_time", &celldf, 0)?)
                    .target_vec(Pos::new(
                        Self::get_col_num_or(&row, "target_vec_x", &celldf, 0.)?,
                        Self::get_col_num_or(&row, "target_vec_y", &celldf, 0.)?,
                    ))
                    .prev_center(Pos::new(
                        Self::get_col_num_or(&row, "prev_center_x", &celldf, 0.)?,
                        Self::get_col_num_or(&row, "prev_center_y", &celldf, 0.)?,
                    ))
                    .cell_type(Self::get_col_str(&row, "cell_type", &celldf)?.try_into()?)
                    .adh_id(0)
                    .build()
            });
        }
        Ok(cells)
    }

    fn get_col_str<'r>(row: &'r Row, col_name: &str, celldf: &DataFrame) -> anyhow::Result<&'r str> {
        let col_index = celldf
            .get_column_index(col_name)
            .ok_or(anyhow::anyhow!("missing `{col_name}`"))?;
        row.0[col_index].get_str().context("could not extract `{col_name}`")
    }

    fn get_col_num<T: NumCast + IsFloat>(row: &Row, col_name: &str, celldf: &DataFrame) -> anyhow::Result<T> {
        let col_index = celldf
            .get_column_index(col_name)
            .ok_or(anyhow::anyhow!("missing `{col_name}`"))?;
        row.0[col_index].try_extract::<T>().context("could not extract `{col_name}`")
    }

    /// Same as [`IoManager::get_col_num()`], but falls back to `default` when the column is absent.
    ///
    /// Used for the columns that were introduced after the data format was first written, so that cell
    /// files and cell templates produced by older versions of the model can still be read.
    fn get_col_num_or<T: NumCast + IsFloat>(
        row: &Row,
        col_name: &str,
        celldf: &DataFrame,
        default: T
    ) -> anyhow::Result<T> {
        let Some(col_index) = celldf.get_column_index(col_name) else {
            return Ok(default);
        };
        row.0[col_index].try_extract::<T>().context("could not extract `{col_name}`")
    }

    /// Reads a cell data file into a [`CellContainer`].
    pub fn read_cells(
        cells_path: impl AsRef<Path>
    ) -> anyhow::Result<CellContainer<MyCell>> {
        let cells_path = cells_path.as_ref();
        let file = std::fs::File::open(cells_path).context(format!("while opening {}", cells_path.display()))?;
        let celldf = ParquetReader::new(file).finish()?;
        Self::make_cells_from_data(
            celldf
        )
    }

    /// Reads a layout file at `layout_path` for a pond with dimensions
    /// `pond_width` and `pond_height` into a gray scale image.
    pub fn read_layout(
        layout_path: impl AsRef<Path>,
        pond_width: usize,
        pond_height: usize
    ) -> anyhow::Result<GrayImage> {
        let layout_path = layout_path.as_ref();
        let layout = ImageReader::open(layout_path)?
            .with_guessed_format()
            .with_context(|| format!("failed to open layout file \"{layout_path:?}\" as PNG"))?
            .decode()?;
        if !matches!(layout.color(), ColorType::L8 | ColorType::L16 | ColorType::La8 | ColorType::La16) {
            log::warn!("Layout file \"{layout_path:?}\" is not encoded in grayscale but will be converted");
        }
        Ok(layout.resize_exact(pond_width as u32, pond_height as u32, FilterType::Nearest).into_luma8())
    }

    fn pad_time_step(time_step: u32) -> String {
        format!("{time_step:0>PAD_FILE_LEN$}")
    }

    /// Given a path to the main folder of a simulation, resolve the path to the file
    /// containing the simulation parameters.
    pub fn resolve_parameters_path(sim_path: impl AsRef<Path>) -> PathBuf {
        sim_path.as_ref().join(CONFIG_COPY_PATH)
    }

    /// Given a path to the main folder of a simulation, resolve the path to the cell data file
    /// that was saved at `time_step`.
    pub fn resolve_cells_path(
        sim_path: impl AsRef<Path>,
        time_step: u32
    ) -> PathBuf {
        sim_path.as_ref()
            .join(CELLS_PATH)
            .join(format!("{}.parquet", Self::pad_time_step(time_step)))
    }

    /// Given a path to the main folder of a simulation, resolve the path to the lattice file
    /// that was saved at `time_step`.
    pub fn resolve_lattice_path(
        sim_path: impl AsRef<Path>,
        time_step: u32
    ) -> PathBuf {
        sim_path.as_ref()
            .join(LATTICES_PATH)
            .join(format!("{}.parquet", Self::pad_time_step(time_step)))
    }

    /// Reads a lattice from a backup file at `file_path`.
    pub fn read_lattice(file_path: impl AsRef<Path>, rect: Rect<usize>) -> anyhow::Result<Lattice<Spin>> {
        let file_path = file_path.as_ref();
        let file = std::fs::File::open(file_path).context(format!("while opening {}", file_path.display()))?;
        let latdf = ParquetReader::new(file).finish()?;
        if latdf.width() != rect.width()
            || latdf.height() != rect.height() {
            bail!("expected lattice dimensions do not match those in file");
        }

        let mut lattice = Lattice::new(rect);
        for (j, column) in latdf.get_columns().iter().enumerate() {
            for (i, maybe_val) in column.str()?.into_iter().enumerate() {
                match maybe_val {
                    Some(val) => {
                        let val: &str = val;
                        let spin = match val {
                            "s" => Spin::Solid,
                            "m" => Spin::Medium,
                            _ => {
                                let cell_index = val.parse::<CellIndex>().with_context(|| {
                                    format!("lattice contains invalid value {val}")
                                })?;
                                Spin::Some(cell_index)
                            },
                        };
                        lattice[(j, i).into()] = spin;
                    },
                    None => bail!("file {} contains null values", file_path.display()),
                }
            }
        }
        Ok(lattice)
    }

    /// Writes both data and simulation images (including movie frames) if its time (according to `time_step`).
    pub fn write_if_time(
        &mut self,
        time_step: u32,
        env: &MyEnvironment
    ) -> anyhow::Result<()> {
        self.write_data_if_time(time_step, env)?;
        self.write_image_if_time(time_step, env)
    }

    fn write_data_if_time(
        &self,
        time_step: u32,
        env: &MyEnvironment
    ) -> anyhow::Result<()> {
        let time_str = Self::pad_time_step(time_step);
        // We might eventually want to buffer the dataframes into an Option<Vec<DF>>
        // and write it less frequently if the volume of files become a problem
        if time_step.is_multiple_of(self.cells_period) {
            let mut celldf = env.env.cells.to_dataframe()?;
            let file_path = self.outdir
                .join(CELLS_PATH)
                .join(format!("{time_str}.parquet"));
            let file = std::fs::File::create(file_path)?;
            ParquetWriter::new(file).finish(&mut celldf)?;
        }

        if time_step.is_multiple_of(self.lattice_period) {
            let file_path = self.outdir
                .join(LATTICES_PATH)
                .join(format!("{time_str}.parquet"));
            Self::write_lattice(file_path.as_path(), &env.env.cell_lattice)?;
        }
        Ok(())
    }

    // Experimented with:
    //   - saving Medium and Solid as negative i32s
    //   - parallelization with rayon
    // and performance diff was minimal and file size became larger, keeping as is
    fn write_lattice(file_path: &Path, lattice: &Lattice<Spin>) -> PolarsResult<u64>{
        let mut cols = vec![];
        for (j, col) in lattice.as_slice().chunks_exact(lattice.height()).enumerate() {
            cols.push(Series::new(
                format!("col_{j}").into(),
                col.iter()
                    .map(|val| {
                        match val {
                            Spin::Solid => "s".into(),
                            Spin::Medium => "m".into(),
                            Spin::Some(cell_index) => cell_index.to_string()
                        }
                    })
                    .collect::<Vec<_>>(),
            ).into())
        }
        let mut latdf = DataFrame::new(cols)?;
        let file = std::fs::File::create(file_path)?;
        ParquetWriter::new(file).finish(&mut latdf)
    }

    fn write_image_if_time(
        &mut self,
        time_step: u32, 
        env: &MyEnvironment
    ) -> anyhow::Result<()> {
        // There might be a way to use LazyCell here but i got tired of fighting the borrow checker
        let mut frame = None;

        #[cfg(feature = "movie")]
        let movie_update = if let Some(mm) = &self.movie_maker {
            time_step.is_multiple_of(mm.frame_period) && mm.window_works()
        } else {
            false
        };
        #[cfg(feature = "movie")]
        if movie_update {
            frame = Some(self.make_simulation_image(env));
            let mm = self.movie_maker.as_mut().unwrap();
            let resized = image::imageops::resize(
                frame.as_ref().unwrap(),
                mm.width,
                mm.height,
                image::imageops::Nearest,
            );
            mm.update(&resized)?
        }

        if time_step.is_multiple_of(self.image_period) {
            if frame.is_none() {
                frame = Some(self.make_simulation_image(env));
            }
            frame.unwrap().save(
                &self.outdir
                    .join(IMAGES_PATH)
                    .join(format!(
                        "{}.{}",
                        Self::pad_time_step(time_step),
                        self.image_format.to_lowercase()
                    ))
            )?;
        }
        Ok(())
    }

    /// Makes a new frame of the simulation by drawing a succession of plots (see [`io::plot`](crate::io::plot)).
    pub fn make_simulation_image(
        &self, 
        env: &MyEnvironment
    ) -> RgbaImage {
        let mut image = self.micro_bg.clone();
        for plot in &self.plots {
            plot.plot(env, &mut image);
        }
        for rel_cell in env.env.cells.iter() {
            if rel_cell.cell.area() == 0 {
                continue;
            }
            let rad = (rel_cell.cell.area() as f64 / PI).sqrt();
            let ratio = rad * 0.05;
            let width = (self.eyes_img.width() as f64 * ratio) as u32;
            let height = (self.eyes_img.height() as f64 * ratio) as u32;
            let eyes = resize(
                &self.eyes_img,
                width,
                height,
                FilterType::Lanczos3,
            );
            // This looked wird but maybe some tilt could be nice
            // if rel_cell.index.is_multiple_of(2) {
            //     flip_horizontal_in_place(&mut eyes);
            // }
            let x = rel_cell.cell.center().x as i64 - (width / 2) as i64;
            let y = rel_cell.cell.center().y as i64 - (height / 2) as i64;
            overlay(
                &mut image,
                &eyes,
                x,
                y
            );
        }
        overlay(&mut image, &self.kinect_img, 0, 0);
        flip_vertical_in_place(&mut image);
        image
    }

    /// Returns the last time step in a simulation directory from which a backup can be restored.
    pub fn find_last_time_step(dir: impl AsRef<Path>) -> anyhow::Result<u32> {
        let dir = dir.as_ref();
        let paths = [CELLS_PATH, LATTICES_PATH];
        let mut intersection = HashSet::new();
        for path in paths {
            let full_path = dir.join(path);
            let file_steps = std::fs::read_dir(full_path)?
                .filter_map(|maybe_file| {
                    let file = maybe_file.ok()?;
                    let file_name = file.file_name();
                    let number_str = file_name.to_str()?.strip_suffix(".parquet")?;
                    number_str.parse::<u32>().ok()
                })
                .collect();

            if intersection.is_empty() {
                intersection = file_steps;
            } else {
                intersection = intersection.intersection(&file_steps).copied().collect();
            }
        }

        intersection
            .into_iter()
            .max()
            .ok_or(anyhow::anyhow!("directory `{dir:?}` does not contain a valid back-up"))
    }
}

trait ToDataFrame {
    fn to_dataframe(&self) -> PolarsResult<DataFrame>;
}

impl ToDataFrame for CellContainer<MyCell> {
    fn to_dataframe(&self) -> PolarsResult<DataFrame> {
        let non_empty = self.iter().filter(|rel_cell| rel_cell.cell.is_empty()).collect::<Box<_>>();
        df!(
            "index" => non_empty.iter().map(|rel_cell| rel_cell.index).collect::<Box<_>>(),
            "area" => non_empty.iter().map(|rel_cell| rel_cell.cell.area()).collect::<Box<_>>(),
            "target_area" => non_empty.iter().map(|rel_cell| rel_cell.cell.target_area()).collect::<Box<_>>(),
            "newborn_target_area" => non_empty.iter().map(|rel_cell| rel_cell.cell.newborn_target_area).collect::<Box<_>>(),
            "divide_area" => non_empty.iter().map(|rel_cell| rel_cell.cell.divide_area).collect::<Box<_>>(),
            "center_x" => non_empty.iter().map(|rel_cell| rel_cell.cell.center().x).collect::<Box<_>>(),
            "center_y" => non_empty.iter().map(|rel_cell| rel_cell.cell.center().y).collect::<Box<_>>(),
            "chem_center_x" => non_empty.iter().map(|rel_cell| rel_cell.cell.chem_center().x).collect::<Box<_>>(),
            "chem_center_y" => non_empty.iter().map(|rel_cell| rel_cell.cell.chem_center().y).collect::<Box<_>>(),
            "chem_mass" => non_empty.iter().map(|rel_cell| rel_cell.cell.chem_mass()).collect::<Box<_>>(),
            "perimeter" => non_empty.iter().map(|rel_cell| rel_cell.cell.perimeter()).collect::<Box<_>>(),
            "target_perimeter" => non_empty.iter().map(|rel_cell| rel_cell.cell.target_perimeter()).collect::<Box<_>>(),
            "persistence_duration" => non_empty.iter().map(|rel_cell| rel_cell.cell.persistence_duration).collect::<Box<_>>(),
            "persistence_time" => non_empty.iter().map(|rel_cell| rel_cell.cell.persistence_time()).collect::<Box<_>>(),
            "target_vec_x" => non_empty.iter().map(|rel_cell| rel_cell.cell.target_vec().x).collect::<Box<_>>(),
            "target_vec_y" => non_empty.iter().map(|rel_cell| rel_cell.cell.target_vec().y).collect::<Box<_>>(),
            "prev_center_x" => non_empty.iter().map(|rel_cell| rel_cell.cell.prev_center().x).collect::<Box<_>>(),
            "prev_center_y" => non_empty.iter().map(|rel_cell| rel_cell.cell.prev_center().y).collect::<Box<_>>(),
            "cell_type" => non_empty.iter().map(|rel_cell| rel_cell.cell.cell_type.to_string()).collect::<Box<[String]>>()
        )
    }
}
