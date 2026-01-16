//! Contains logic related to [`IoManager`].

use crate::cell::{Cell, CellType};
use crate::environment::Environment;
#[cfg(feature = "movie")]
use crate::io::movie_maker::MovieMaker;
use crate::io::parameters::Parameters;
use crate::io::plot::Plot;
use anyhow::{bail, Context};
use bon::Builder;
use cellulars_lib::base::base_cell::BaseCell;
use cellulars_lib::cell_container::{CellContainer, RelCell};
use cellulars_lib::constants::CellIndex;
use cellulars_lib::lattice::Lattice;
use cellulars_lib::positional::com::Com;
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::positional::rect::Rect;
use cellulars_lib::spin::Spin;
use cellulars_lib::traits::cellular::Cellular;
use image::imageops::{flip_vertical_in_place, FilterType};
use image::{ColorType, GrayImage, ImageReader, RgbaImage};
use num_traits::NumCast;
use polars::frame::row::Row;
use polars::polars_utils::float::IsFloat;
use polars::prelude::*;
use std::collections::HashSet;
use std::io;
use std::path::{Path, PathBuf};
use crate::io::kinect_listener::KinectListener;

static IMAGES_PATH: &str = "images";
static CELLS_PATH: &str = "cells";
static LATTICES_PATH: &str = "lattices";
static CONFIG_COPY_PATH: &str = "config.toml";

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

    fn make_cells_from_data(celldf: DataFrame) -> anyhow::Result<CellContainer<Cell>> {
        let last_index = celldf
            .column("index")?
            .u32()?
            .max()
            .ok_or(anyhow::anyhow!("null `index` column"))?;
        let mut cells = CellContainer::new();
        // We need this to call replace on cells later
        for _ in 0..=last_index {
            cells.push(Cell::new_empty(0, 0, CellType::Migrating));
        }

        for row_i in 0..celldf.height() {
            let row = celldf.get_row(row_i)?;
            let base_cell = BaseCell::new_ready(
                Self::get_col_num(&row, "area", &celldf)?,
                Self::get_col_num(&row, "target_area", &celldf)?,
                Pos::new(
                    Self::get_col_num(&row, "center_x", &celldf)?,
                    Self::get_col_num(&row, "center_y", &celldf)?,
                )
            );
            cells.replace(RelCell {
                index: Self::get_col_num(&row, "index", &celldf)?,
                cell: Cell::builder()
                    .base_cell(base_cell)
                    .divide_area(Self::get_col_num(&row, "divide_area", &celldf)?)
                    .newborn_target_area(Self::get_col_num(&row, "newborn_target_area", &celldf)?)
                    .chem_com(Com {
                        pos: Pos::new(
                            Self::get_col_num(&row, "chem_center_x", &celldf)?,
                            Self::get_col_num(&row, "chem_center_y", &celldf)?,
                        ),
                        mass: Self::get_col_num(&row, "chem_mass", &celldf)?
                    })
                    .cell_type(Self::get_col_str(&row, "cell_type", &celldf)?.try_into()?)
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

    /// Reads a cell data file into a [`CellContainer`].
    pub fn read_cells(
        cells_path: impl AsRef<Path>
    ) -> anyhow::Result<CellContainer<Cell>> {
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
            .join(format!("{time_step}.parquet"))
    }

    /// Given a path to the main folder of a simulation, resolve the path to the lattice file
    /// that was saved at `time_step`.
    pub fn resolve_lattice_path(
        sim_path: impl AsRef<Path>,
        time_step: u32
    ) -> PathBuf {
        sim_path.as_ref()
            .join(LATTICES_PATH)
            .join(format!("{time_step}.parquet"))
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
        env: &Environment
    ) -> anyhow::Result<()> {
        self.write_data_if_time(time_step, env)?;
        self.write_image_if_time(time_step, env)
    }

    fn write_data_if_time(
        &self,
        time_step: u32,
        env: &Environment
    ) -> anyhow::Result<()> {
        let time_str = time_step.to_string();
        // We might eventually want to buffer the dataframes into an Option<Vec<DF>>
        // and write it less frequently if the volume of files become a problem
        if time_step.is_multiple_of(self.cells_period) {
            let mut celldf = env.base_env.cells.to_dataframe()?;
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
            Self::write_lattice(file_path.as_path(), &env.base_env.cell_lattice)?;
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
        env: &Environment
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
                    .join(format!("{time_step}.{}", self.image_format.to_lowercase())
                    ))?;
        }
        Ok(())
    }

    /// Makes a new frame of the simulation by drawing a succession of plots (see [`io::plot`](crate::io::plot)).
    pub fn make_simulation_image(
        &self, 
        env: &Environment
    ) -> RgbaImage {
        let mut image = RgbaImage::new(
            env.base_env.width() as u32,
            env.base_env.height() as u32
        );
        for plot in &self.plots {
            plot.plot(env, &mut image);
        }
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

impl ToDataFrame for CellContainer<Cell> {
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
            "cell_type" => non_empty.iter().map(|rel_cell| rel_cell.cell.cell_type.to_string()).collect::<Box<[String]>>()
        )
    }
}
