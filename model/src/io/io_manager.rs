//! Contains logic related to [IoManager].

use crate::cell::{Cell, CellType};
use crate::io::movie_maker::MovieMaker;
use crate::io::parameters::Parameters;
use crate::io::plot::Plot;
use crate::my_environment::MyEnvironment;
use anyhow::{anyhow, bail, Context};
use bon::Builder;
use cellulars_lib::basic_cell::{BasicCell, Cellular, RelCell};
use cellulars_lib::cell_container::CellContainer;
use cellulars_lib::constants::CellIndex;
use cellulars_lib::lattice::Lattice;
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::positional::rect::Rect;
use cellulars_lib::spin::Spin;
use image::imageops::flip_vertical_in_place;
use image::RgbaImage;
use polars::prelude::*;
use std::collections::HashMap;
use std::io;
use std::path::{Path, PathBuf};

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
    pub movie_maker: Option<MovieMaker>,
    plots: Vec<Box<dyn Plot>>,
    image_period: u32,
    cells_period: u32,
    lattices_period: u32
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

    /// Creates a parameter file at \"[IoManager::outdir]/config.toml\".
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
        let mut cells = CellContainer::new();
        // We need this to call replace on cells later
        for _ in 0..=celldf.height() {
            cells.push(Cell::new_empty(0, 0, CellType::Migrating));
        }

        let cols: HashMap<_, _> = HashMap::from_iter(
            celldf.get_column_names()
                .into_iter()
                .enumerate()
                .map(|(i, name)| (name.as_str(), i))
        );

        for row_i in 0..celldf.height() {
            let row = celldf.get_row(row_i)?.0;
            let cell_type = row[cols["cell_type"]]
                .get_str()
                .ok_or(anyhow!("could not extract cell type string"))?
                .try_into()?;
            let basic_cell = BasicCell::new_ready(
                row[cols["area"]].try_extract::<u32>()?,
                Pos::new(
                    row[cols["center_x"]].try_extract::<f32>()?,
                    row[cols["center_y"]].try_extract::<f32>()?,
                ),
                row[cols["target_area"]].try_extract::<u32>()?
            );
            cells.replace(RelCell {
                index: row[cols["index"]].try_extract::<CellIndex>()?,
                cell: Cell::builder()
                    .basic_cell(basic_cell)
                    .divide_area(row[cols["divide_area"]].try_extract::<u32>()?)
                    .newborn_target_area(row[cols["newborn_target_area"]].try_extract::<u32>()?)
                    .chem_center(Pos::new(
                        row[cols["chem_center_x"]].try_extract::<f32>()?,
                        row[cols["chem_center_y"]].try_extract::<f32>()?, 
                    ))
                    .chem_mass(row[cols["chem_mass"]].try_extract::<u32>()?)
                    .cell_type(cell_type)
                    .build()
            });
        }
        Ok(cells)
    }

    /// Reads a cell data file into a [CellContainer].
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
        let time_str = time_step.to_string();
        // We might eventually want to buffer the dataframes into an Option<Vec<DF>>
        // and write it less frequently if the volume of files become a problem
        if time_step % self.cells_period == 0 {
            let mut celldf = env.cells.to_dataframe()?;
            let file_path = self.outdir
                .join(CELLS_PATH)
                .join(format!("{time_str}.parquet"));
            let file = std::fs::File::create(file_path)?;
            ParquetWriter::new(file).finish(&mut celldf)?;
        }

        if time_step % self.lattices_period == 0 {
            let file_path = self.outdir
                .join(LATTICES_PATH)
                .join(format!("{time_str}.parquet"));
            Self::write_lattice(file_path.as_path(), &env.cell_lattice)?;
        }
        Ok(())
    }

    // Experimented with:
    //   - saving Medium and Solid as negative i32s
    //   - parallelisation with rayon
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
        let movie_update = if let Some(mm) = &self.movie_maker {
            time_step % mm.frame_period == 0 && mm.window_works()
        } else {
            false
        };
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

        if time_step % self.image_period == 0 {
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

    /// Makes a new frame of the simulation by drawing a succession of plots (see [io::plot](crate::io::plot)).
    pub fn make_simulation_image(
        &self, 
        env: &MyEnvironment
    ) -> RgbaImage {
        let mut image = RgbaImage::new(
            env.width() as u32,
            env.height() as u32 
        );
        for plot in &self.plots {
            plot.plot(env, &mut image);
        }
        flip_vertical_in_place(&mut image);
        image
    }
}

trait ToDataFrame {
    fn to_dataframe(&self) -> PolarsResult<DataFrame>;
}

impl ToDataFrame for CellContainer<Cell> {
    fn to_dataframe(&self) -> PolarsResult<DataFrame> {
        let valid = self.iter().filter(|cell| cell.is_valid()).collect::<Vec<_>>();
        df!(
            "index" => valid.iter().map(|cell| cell.index).collect::<Vec<_>>(),
            "area" => valid.iter().map(|cell| cell.area()).collect::<Vec<_>>(),
            "target_area" => valid.iter().map(|cell| cell.target_area()).collect::<Vec<_>>(),
            "newborn_target_area" => valid.iter().map(|cell| cell.newborn_target_area).collect::<Vec<_>>(),
            "divide_area" => valid.iter().map(|cell| cell.divide_area).collect::<Vec<_>>(),
            "center_x" => valid.iter().map(|cell| cell.center().x).collect::<Vec<_>>(),
            "center_y" => valid.iter().map(|cell| cell.center().y).collect::<Vec<_>>(),
            "chem_center_x" => valid.iter().map(|cell| cell.chem_center().x).collect::<Vec<_>>(),
            "chem_center_y" => valid.iter().map(|cell| cell.chem_center().y).collect::<Vec<_>>(),
            "chem_mass" => valid.iter().map(|cell| cell.chem_mass()).collect::<Vec<_>>(),
            "cell_type" => valid.iter().map(|cell| cell.cell_type.to_string()).collect::<Vec<_>>()
        )
    }
}
