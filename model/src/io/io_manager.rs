use crate::cell::Cell;
use crate::evolution::grn::{EdgeWeight, Grn, GrnGeneType};
use crate::io::movie_maker::MovieMaker;
use crate::io::node_link::{GrnMutParams, NodeLinkData};
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
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::io;
use std::io::BufReader;
use std::path::{Path, PathBuf};

static IMAGES_PATH: &str = "images";
static CELLS_PATH: &str = "cells";
static GENOMES_PATH: &str = "genomes";
static CELL_LATTICES_PATH: &str = "lattices";
static CHEM_LATTICES_PATH: &str = "chem_lattices";
static ACT_LATTICES_PATH: &str = "act_lattices";
static CONFIG_COPY_PATH: &str = "config.toml";

#[derive(Builder)]
pub struct IoManager {
    pub outdir: PathBuf,
    pub image_format: String,
    pub movie_maker: Option<MovieMaker>,
    plots: Vec<Box<dyn Plot>>,
    image_period: u32,
    cells_period: u32,
    genomes_period: u32,
    lattices_period: u32
}

impl IoManager {
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
        std::fs::create_dir(self.outdir.join(GENOMES_PATH))?;
        std::fs::create_dir(self.outdir.join(CELL_LATTICES_PATH))?;
        std::fs::create_dir(self.outdir.join(CHEM_LATTICES_PATH))?;
        std::fs::create_dir(self.outdir.join(ACT_LATTICES_PATH))
    }

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

    fn make_cells_from_data(
        celldf: DataFrame,
        genomes: Vec<CellIndexNodeLink>,
    ) -> anyhow::Result<CellContainer<Cell>> {
        if celldf.height() > genomes.len() {
            bail!("`celldf` contains more entries than `genomes`");
        }
        if celldf.height() < genomes.len() {
            bail!("`genomes` contains more entries than `celldf`");
        }

        let mut cells = CellContainer::new();
        // We need this to call replace on cells later
        for _ in 0..=celldf.height() {
            cells.push(Cell::new_empty(0, 0, Grn::empty()));
        }

        let index_map: HashMap<_, _> = HashMap::from_iter(
            celldf.column("index")?
                .u32()?
                .into_no_null_iter()
                .enumerate()
                .map(|(i, val)| {
                    (val, i)
                })
        );

        let cols: HashMap<_, _> = HashMap::from_iter(
            celldf.get_column_names()
                .into_iter()
                .enumerate()
                .map(|(i, name)| (name.as_str(), i))
        );

        for genome in genomes {
            let cell_index = genome.index;
            let grn = Grn::<1, 1>::try_from(genome.node_link)?;
            let row_ix = index_map
                .get(&cell_index)
                .ok_or(anyhow!("cell index {cell_index} was found in `genomes` but is missing from `celldf`"))?;
            let row = celldf.get_row(*row_ix)?.0;

            cells.replace(RelCell {
                index: cell_index,
                cell: Cell {
                    basic_cell: BasicCell {
                        target_area: row[cols["target_area"]].try_extract::<u32>()?,
                        area: row[cols["area"]].try_extract::<u32>()?,
                        center: Pos::new(
                            row[cols["center_x"]].try_extract::<f32>()?,
                            row[cols["center_y"]].try_extract::<f32>()?,
                        )
                    },
                    perimeter: row[cols["perimeter"]].try_extract::<u32>()?,
                    target_perimeter: row[cols["target_perimeter"]].try_extract::<u32>()?,
                    delta_perimeter: None,
                    ancestor: Some(row[cols["ancestor"]].try_extract::<CellIndex>()?),
                    chem_center: Pos::new(
                        row[cols["chem_center_x"]].try_extract::<f32>()?,
                        row[cols["chem_center_y"]].try_extract::<f32>()?,
                    ),
                    chem_mass: row[cols["chem_mass"]].try_extract::<u32>()?,
                    genome: grn,
                }
            });
        }
        Ok(cells)
    }

    pub fn read_cells(
        cells_path: impl AsRef<Path>,
        genomes_path: impl AsRef<Path>,
    ) -> anyhow::Result<CellContainer<Cell>> {
        let cells_path = cells_path.as_ref();
        let file = std::fs::File::open(cells_path).context(format!("while opening {}", cells_path.display()))?;
        let celldf = ParquetReader::new(file).finish()?;
        let genomes = Self::read_genomes(genomes_path)?;
        Self::make_cells_from_data(
            celldf,
            genomes
        )
    }

    pub fn resolve_parameters_path(sim_path: impl AsRef<Path>) -> PathBuf {
        sim_path.as_ref().join(CONFIG_COPY_PATH)
    }

    pub fn resolve_cells_path(
        sim_path: impl AsRef<Path>,
        time_step: u32
    ) -> PathBuf {
        sim_path.as_ref()
            .join(CELLS_PATH)
            .join(format!("{time_step}.parquet"))
    }

    pub fn resolve_genomes_path(
        sim_path: impl AsRef<Path>,
        time_step: u32
    ) -> PathBuf {
        sim_path.as_ref()
            .join(GENOMES_PATH)
            .join(format!("{time_step}.json"))
    }

    fn read_genomes(file_path: impl AsRef<Path>) -> anyhow::Result<Vec<CellIndexNodeLink>> {
        let file_path = file_path.as_ref();
        let file = std::fs::File::open(file_path).context(format!("while opening {}", file_path.display()))?;
        let reader = BufReader::new(file);
        Ok(serde_json::from_reader(reader)?)
    }

    pub fn resolve_cell_lattice_path(
        sim_path: impl AsRef<Path>,
        time_step: u32
    ) -> PathBuf {
        sim_path.as_ref()
            .join(CELL_LATTICES_PATH)
            .join(format!("{time_step}.parquet"))
    }

    pub fn resolve_chem_lattice_path(
        sim_path: impl AsRef<Path>,
        time_step: u32
    ) -> PathBuf {
        sim_path.as_ref()
            .join(CHEM_LATTICES_PATH)
            .join(format!("{time_step}.parquet"))
    }

    pub fn resolve_act_lattice_path(
        sim_path: impl AsRef<Path>,
        time_step: u32
    ) -> PathBuf {
        sim_path.as_ref()
            .join(ACT_LATTICES_PATH)
            .join(format!("{time_step}.parquet"))
    }

    fn read_ladf(file_path: impl AsRef<Path>, rect: &Rect<usize>) -> anyhow::Result<DataFrame> {
        let file_path = file_path.as_ref();
        let file = std::fs::File::open(file_path).context(format!("while opening {}", file_path.display()))?;
        let latdf = ParquetReader::new(file).finish()?;
        if latdf.width() != rect.width()
            || latdf.height() != rect.height() {
            bail!("expected lattice dimensions do not match those in file");
        }
        Ok(latdf)
    }

    pub fn read_cell_lattice(file_path: impl AsRef<Path>, rect: Rect<usize>) -> anyhow::Result<Lattice<Spin>> {
        let file_path = file_path.as_ref();
        let latdf = Self::read_ladf(file_path, &rect)?;
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

    pub fn read_lattice_u32(file_path: impl AsRef<Path>, rect: Rect<usize>) -> anyhow::Result<Lattice<u32>> {
        let file_path = file_path.as_ref();
        let latdf = Self::read_ladf(file_path, &rect)?;
        let mut lattice = Lattice::new(rect);
        for (j, column) in latdf.get_columns().iter().enumerate() {
            for (i, maybe_val) in column.u32()?.into_iter().enumerate() {
                match maybe_val {
                    Some(val) => {
                        lattice[(j, i).into()] = val;
                    },
                    None => bail!("file {} contains null values", file_path.display()),
                }
            }
        }
        Ok(lattice)
    }

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

        if time_step % self.genomes_period == 0 {
            let genomes = Self::make_genome_node_links(&env.cells);
            let file_path = self.outdir
                .join(GENOMES_PATH)
                .join(format!("{time_str}.json"));
            let file = std::fs::File::create(file_path)?;
            serde_json::to_writer(file, &genomes)?;
        }

        if time_step % self.lattices_period == 0 {
            let cell_lat_file_path = self.outdir
                .join(CELL_LATTICES_PATH)
                .join(format!("{time_str}.parquet"));
            Self::write_lattice(cell_lat_file_path.as_path(), &env.cell_lattice)?;

            let chem_lat_file_path = self.outdir
                .join(CHEM_LATTICES_PATH)
                .join(format!("{time_str}.parquet"));
            Self::write_lattice_u32(chem_lat_file_path.as_path(), &env.chem_lattice)?;

            let act_lat_file_path = self.outdir
                .join(ACT_LATTICES_PATH)
                .join(format!("{time_str}.parquet"));
            Self::write_lattice_u32(act_lat_file_path.as_path(), &env.act_lattice)?;
        }
        Ok(())
    }

    fn make_genome_node_links(
        cells: &CellContainer<Cell>
    ) -> Vec<CellIndexNodeLink> {
        cells.iter()
            .filter(|cell| cell.is_valid())
            .map(|cell| CellIndexNodeLink {
                index: cell.index,
                node_link: NodeLinkData::from(cell.genome.clone())
            })
            .collect()
    }

    // Experimented with:
    //   - saving Medium and Solid as negative i32s
    //   - parallelisation with rayon
    // and performance diff was minimal and file size became larger, keeping as is
    fn write_lattice(file_path: &Path, lattice: &Lattice<Spin>) -> PolarsResult<u64>{
        let mut cols = vec![];
        for (j, col) in lattice.as_array().chunks_exact(lattice.height()).enumerate() {
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

    fn write_lattice_u32(file_path: &Path, lattice: &Lattice<u32>) -> PolarsResult<u64>{
        let mut cols = vec![];
        for (j, col) in lattice.as_array().chunks_exact(lattice.height()).enumerate() {
            cols.push(Series::new(
                format!("col_{j}").into(),
                col.to_vec(),
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
            "ancestor" => valid.iter().map(|cell| cell.ancestor).collect::<Vec<_>>(),
            "area" => valid.iter().map(|cell| cell.area()).collect::<Vec<_>>(),
            "target_area" => valid.iter().map(|cell| cell.target_area()).collect::<Vec<_>>(),
            "perimeter" => valid.iter().map(|cell| cell.perimeter).collect::<Vec<_>>(),
            "target_perimeter" => valid.iter().map(|cell| cell.target_perimeter).collect::<Vec<_>>(),
            "center_x" => valid.iter().map(|cell| cell.center().x).collect::<Vec<_>>(),
            "center_y" => valid.iter().map(|cell| cell.center().y).collect::<Vec<_>>(),
            "chem_center_x" => valid.iter().map(|cell| cell.chem_center.x).collect::<Vec<_>>(),
            "chem_center_y" => valid.iter().map(|cell| cell.chem_center.y).collect::<Vec<_>>(),
            "chem_mass" => valid.iter().map(|cell| cell.chem_mass).collect::<Vec<_>>(),
            "is_dividing" => valid.iter().map(|cell| cell.is_dividing()).collect::<Vec<_>>()
        )
    }
}

#[derive(Serialize, Deserialize)]
struct CellIndexNodeLink {
    index: CellIndex,
    node_link: NodeLinkData<GrnGeneType, EdgeWeight, GrnMutParams>
}
