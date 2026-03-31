use crate::cell::Cell;
use crate::evolution::bit_genome::BitGenome;
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
use cellulars_lib::spin::{spin_to_str, Spin};
use image::imageops::{flip_vertical_in_place, FilterType};
use image::{ColorType, GrayImage, ImageReader, RgbaImage};
use polars::prelude::*;
use std::collections::{HashMap, HashSet};
use std::io;
use std::path::{Path, PathBuf};

static IMAGES_PATH: &str = "images";
static CELLS_PATH: &str = "cells";
static CELL_LATTICES_PATH: &str = "lattices";
static CHEM_LATTICES_PATH: &str = "chem_lattices";
static ACT_LATTICES_PATH: &str = "act_lattices";
static ACT_CONTACTS_PATH: &str = "act_contacts";
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

#[derive(Builder)]
pub struct IoManager {
    pub outdir: PathBuf,
    pub image_format: String,
    pub movie_maker: Option<MovieMaker>,
    plots: Box<[Box<dyn Plot>]>,
    image_period: u32,
    cells_period: u32,
    cells_write_period: u32,
    lattices_period: u32,
    act_contacts_period: u32,
    #[builder(default)]
    celldfs: Vec<LazyFrame>,
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
        std::fs::create_dir(self.outdir.join(CELL_LATTICES_PATH))?;
        std::fs::create_dir(self.outdir.join(CHEM_LATTICES_PATH))?;
        std::fs::create_dir(self.outdir.join(ACT_LATTICES_PATH))?;
        std::fs::create_dir(self.outdir.join(ACT_CONTACTS_PATH))
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

    /// Returns the last time step in a simulation directory from which a backup can be restored.
    pub fn find_last_time_step(dir: impl AsRef<Path>) -> anyhow::Result<u32> {
        let dir = dir.as_ref();
        let paths = [CELLS_PATH, ACT_LATTICES_PATH, CELL_LATTICES_PATH, CHEM_LATTICES_PATH];
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
        let resized = layout.resize_exact(
            pond_width as u32,
            pond_height as u32,
            FilterType::Nearest
        );
        let flipped = resized.flipv();
        Ok(flipped.into_luma8())
    }

    fn make_cells_from_data(
        celldf: DataFrame,
        mut_rate: f32,
        genome_length: u8
    ) -> anyhow::Result<CellContainer<Cell>> {
        let mut cells = CellContainer::new();
        // We need this to call replace on cells later
        for _ in 0..=celldf.column("index")?.u32()?.max().ok_or(anyhow::anyhow!("null column"))? {
            cells.push(Cell::new_empty(
                0,
                0,
                BitGenome::new(
                    0,
                    0,
                    0.,
                    1
                ).unwrap(),
            ));
        }

        let cols: HashMap<_, _> = HashMap::from_iter(
            celldf.get_column_names()
                .into_iter()
                .enumerate()
                .map(|(i, name)| (name.as_str(), i))
        );

        for i in 0..celldf.height() {
            let row = celldf.get_row(i)?.0;
            cells.replace(RelCell {
                index: row[cols["index"]].try_extract::<u32>()?,
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
                    neighbors: HashMap::new(),
                    tot_act: 0,
                    tot_kact: 0.,
                    rel_chem: 0.0,
                    genome: BitGenome::new(
                        row[cols["ligands"]].try_extract::<u64>()?,
                        row[cols["receptors"]].try_extract::<u64>()?,
                        mut_rate,
                        genome_length
                    ).ok_or(anyhow!("invalid `genome_length`"))?,
                }
            });
        }
        Ok(cells)
    }

    pub fn read_cells(
        cells_path: impl AsRef<Path>,
        mut_rate: f32,
        genome_length: u8
    ) -> anyhow::Result<CellContainer<Cell>> {
        let cells_path = cells_path.as_ref();
        let file = std::fs::File::open(cells_path).context(format!("while opening {}", cells_path.display()))?;
        let celldf = ParquetReader::new(file).finish()?;
        Self::make_cells_from_data(celldf, mut_rate, genome_length)
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
            .join(format!("{}.parquet", Self::pad_time_step(time_step)))
    }

    pub fn resolve_cell_lattice_path(
        sim_path: impl AsRef<Path>,
        time_step: u32
    ) -> PathBuf {
        sim_path.as_ref()
            .join(CELL_LATTICES_PATH)
            .join(format!("{}.parquet", Self::pad_time_step(time_step)))
    }

    pub fn resolve_chem_lattice_path(
        sim_path: impl AsRef<Path>,
        time_step: u32
    ) -> PathBuf {
        sim_path.as_ref()
            .join(CHEM_LATTICES_PATH)
            .join(format!("{}.parquet", Self::pad_time_step(time_step)))
    }

    pub fn resolve_act_lattice_path(
        sim_path: impl AsRef<Path>,
        time_step: u32
    ) -> PathBuf {
        sim_path.as_ref()
            .join(ACT_LATTICES_PATH)
            .join(format!("{}.parquet", Self::pad_time_step(time_step)))
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
        env: &mut MyEnvironment
    ) -> anyhow::Result<()> {
        self.write_data_if_time(time_step, env)?;
        self.write_image_if_time(time_step, env)
    }

    fn write_data_if_time(
        &mut self,
        time_step: u32,
        env: &mut MyEnvironment
    ) -> anyhow::Result<()> {
        if time_step.is_multiple_of(self.cells_period) {
            env.update_neighbours();
            env.update_act();
            let mut celldf = env
                .to_dataframe()
                .with_context(|| "failed to make data frame from cells")?;
            celldf
                .with_column(Series::new("time".into(), vec![time_step; celldf.height()]))
                .with_context(|| "failed to add time column to cell data frame")?;
            self.celldfs.push(celldf.lazy());
            env.reset_act();
        }

        // We might eventually want to buffer the dataframes into an Option<Vec<DF>>
        // and write it less frequently if the volume of files become a problem
        if time_step.is_multiple_of(self.cells_write_period) {
            let file_path = self.outdir
                .join(CELLS_PATH)
                .join(format!("{}.parquet", Self::pad_time_step(time_step)));
            let file = std::fs::File::create(file_path)?;

            let newdfs = std::mem::take(&mut self.celldfs);
            let mut celldf = concat(newdfs, UnionArgs::default())
                .and_then(|df| df.collect())
                .with_context(|| "failed to concatenate cell data frames")?;

            ParquetWriter::new(file).finish(&mut celldf)?;

            for cell in env.cells.iter_mut() {
                if !cell.is_valid() {
                    continue;
                }
                cell.ancestor = Some(cell.index);
            }
        }

        if time_step.is_multiple_of(self.lattices_period) {
            let file_name = format!("{}.parquet", Self::pad_time_step(time_step));
            let cell_lat_file_path = self.outdir
                .join(CELL_LATTICES_PATH)
                .join(&file_name);
            Self::write_lattice(cell_lat_file_path.as_path(), &env.cell_lattice)?;

            let chem_lat_file_path = self.outdir
                .join(CHEM_LATTICES_PATH)
                .join(&file_name);
            Self::write_lattice_u32(chem_lat_file_path.as_path(), &env.chem_lattice)?;

            let act_lat_file_path = self.outdir
                .join(ACT_LATTICES_PATH)
                .join(&file_name);
            Self::write_lattice_u32(act_lat_file_path.as_path(), &env.act_lattice)?;
        }

        if time_step.is_multiple_of(self.act_contacts_period) {
            let file_name = format!("{}.parquet", Self::pad_time_step(time_step));
            let file_path = self.outdir
                .join(ACT_CONTACTS_PATH)
                .join(&file_name);

            let (spins, neighs, acts, kacts) = Self::act_contacts(env);
            let mut actdf = df![
                "spin" => spins,
                "neigh" => neighs,
                "act" => acts,
                "kact" => kacts,
            ]?;

            let file = std::fs::File::create(file_path)?;
            ParquetWriter::new(file).finish(&mut actdf)?;
        }
        Ok(())
    }

    fn act_contacts(
        env: &MyEnvironment,
    ) -> (Vec<String>, Vec<String>, Vec<u32>, Vec<f64>) {
        let mut spins = Vec::new();
        let mut neighs = Vec::new();
        let mut acts = Vec::new();
        let mut kacts = Vec::new();
        for edge in env.edge_book.iter() {
            let mut spin1 = env.cell_lattice[edge.p1];
            let mut spin2 = env.cell_lattice[edge.p2];
            if !matches!(spin1, Spin::Some(_)) {
                std::mem::swap(&mut spin1, &mut spin2);
            }
            spins.push(spin_to_str(spin1));
            neighs.push(spin_to_str(spin2));
            acts.push(env.act_lattice[edge.p1]);
            kacts.push(env.kact(edge.p1));

            if let Spin::Some(_) = spin2 {
                spins.push(spin_to_str(spin2));
                neighs.push(spin_to_str(spin1));
                acts.push(env.act_lattice[edge.p2]);
                kacts.push(env.kact(edge.p2));
            }
        }
        (spins, neighs, acts, kacts)
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
                    .map(|val| spin_to_str(*val))
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
            time_step.is_multiple_of(mm.frame_period) && mm.window_works()
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

    fn pad_time_step(time_step: u32) -> String {
        format!("{time_step:0>PAD_FILE_LEN$}")
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

impl ToDataFrame for MyEnvironment {
    fn to_dataframe(&self) -> PolarsResult<DataFrame> {
        let valid = self.cells.iter().filter(|cell| cell.is_valid()).collect::<Box<_>>();
        df!(
            "index" => valid.iter().map(|cell| cell.index).collect::<Box<_>>(),
            "ancestor" => valid.iter().map(|cell| cell.ancestor).collect::<Box<_>>(),
            "area" => valid.iter().map(|cell| cell.area()).collect::<Box<_>>(),
            "target_area" => valid.iter().map(|cell| cell.target_area()).collect::<Box<_>>(),
            "perimeter" => valid.iter().map(|cell| cell.perimeter).collect::<Box<_>>(),
            "target_perimeter" => valid.iter().map(|cell| cell.target_perimeter).collect::<Box<_>>(),
            "center_x" => valid.iter().map(|cell| cell.center().x).collect::<Box<_>>(),
            "center_y" => valid.iter().map(|cell| cell.center().y).collect::<Box<_>>(),
            "chem_center_x" => valid.iter().map(|cell| cell.chem_center.x).collect::<Box<_>>(),
            "chem_center_y" => valid.iter().map(|cell| cell.chem_center.y).collect::<Box<_>>(),
            "chem_mass" => valid.iter().map(|cell| cell.chem_mass).collect::<Box<_>>(),
            "ligands" => valid.iter().map(|cell| cell.genome.ligands()).collect::<Box<_>>(),
            "receptors" => valid.iter().map(|cell| cell.genome.receptors()).collect::<Box<_>>(),
            "neighbors" => valid.iter().map(|cell| cell.neighbors.keys().map(|v| spin_to_str(*v)).collect::<Box<[String]>>().join(" ")).collect::<Box<[String]>>(),
            "neighbor_contacts" => valid.iter().map(|cell| cell.neighbors.values().map(|v| v.to_string()).collect::<Box<[String]>>().join(" ")).collect::<Box<[String]>>(),
            "med_neighbor" => valid.iter().map(|cell| cell.neighbors.contains_key(&Spin::Medium)).collect::<Box<_>>(),
            "solid_neighbor" => valid.iter().map(|cell| cell.neighbors.contains_key(&Spin::Solid)).collect::<Box<_>>(),
            "tot_act" => valid.iter().map(|cell| cell.tot_act).collect::<Box<_>>(),
            "tot_kact" => valid.iter().map(|cell| cell.tot_kact).collect::<Box<_>>(),
        )
    }
}
