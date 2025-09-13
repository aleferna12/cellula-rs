use crate::cell::Cell;
use crate::genetics::grn::{EdgeWeight, Grn, GrnGeneType};
use crate::io::movie_maker::MovieMaker;
use crate::io::node_link::{GrnMutParams, NodeLinkData};
use crate::io::parameters::{Parameters, PlotParameters, PlotType};
use crate::io::plot::*;
use crate::pond::Pond;
use anyhow::{anyhow, bail, Context};
use bon::Builder;
use cellulars_lib::basic_cell::{BasicCell, Cellular, RelCell};
use cellulars_lib::cell_container::CellContainer;
use cellulars_lib::constants::Spin;
use cellulars_lib::environment::Habitable;
use cellulars_lib::lattice::Lattice;
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::positional::rect::Rect;
use image::imageops::flip_vertical_in_place;
use image::{GenericImage, RgbaImage};
use polars::prelude::*;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt::Display;
use std::io;
use std::io::{BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

static IMAGES_PATH: &str = "images";
static CELLS_PATH: &str = "cells";
static GENOMES_PATH: &str = "genomes";
static LATTICES_PATH: &str = "lattices";
static CONFIG_COPY_PATH: &str = "config.toml";
static POND_PREFIX: &str = "pond_";

#[derive(Builder)]
pub struct IoManager {
    pub outdir: PathBuf,
    pub image_format: String,
    pub movie_maker: Option<MovieMaker>,
    // TODO: Can we do something smarter about this? Maybe using dynamic dispatch
    pub plots: PlotParameters,
    image_period: u32,
    cell_period: u32,
    genome_period: u32,
    lattice_period: u32,
}

impl IoManager {
    pub fn create_directories(&self, replace_outdir: bool, n_ponds: u32) -> io::Result<()> {
        let outdir_exists = self.outdir.try_exists()?;
        if outdir_exists {
            if replace_outdir {
                log::info!("Cleaning contents of '{}'", self.outdir.display());
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

        for i in 0..n_ponds {
            let pond_path = self.outdir.join(format!("{POND_PREFIX}{i}"));
            std::fs::create_dir(&pond_path)?;
            std::fs::create_dir(pond_path.join(CELLS_PATH))?;
            std::fs::create_dir(pond_path.join(GENOMES_PATH))?;
            std::fs::create_dir(pond_path.join(LATTICES_PATH))?
        }
        Ok(())
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
        genomes: Vec<SpinNodeLink>,
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
            cells.push(Cell::new_empty(0, 0, Grn::empty()), None);
        }

        let spin_map: HashMap<_, _> = HashMap::from_iter(
            celldf.column("spin")?
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
            let spin = genome.spin;
            let grn = Grn::<1, 1>::try_from(genome.node_link)?;
            let row_ix = spin_map
                .get(&spin)
                .ok_or(anyhow!("spin {spin} was found in `genomes` but is missing from `celldf`"))?;
            let row = celldf.get_row(*row_ix)?.0;

            cells.replace(RelCell {
                spin,
                mom: row[cols["mom"]].try_extract::<Spin>()?,
                cell: Cell {
                    basic_cell: BasicCell {
                        target_area: row[cols["target_area"]].try_extract::<u32>()?,
                        newborn_target_area: row[cols["newborn_target_area"]].try_extract::<u32>()?,
                        area: row[cols["area"]].try_extract::<u32>()?,
                        center: Pos::new(
                            row[cols["center_x"]].try_extract::<f32>()?,
                            row[cols["center_y"]].try_extract::<f32>()?,
                        )
                    },
                    divide_area: row[cols["divide_area"]].try_extract::<u32>()?,
                    chem_center: Pos::new(
                        row[cols["chem_center_x"]].try_extract::<f32>()?,
                        row[cols["chem_center_y"]].try_extract::<f32>()?,
                    ),
                    chem_mass: row[cols["chem_mass"]].try_extract::<f32>()?,
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
        time_step: u32,
        pond_i: u32
    ) -> PathBuf {
        sim_path.as_ref()
            .join(format!("{POND_PREFIX}{pond_i}"))
            .join(CELLS_PATH)
            .join(format!("{time_step}.parquet"))
    }

    pub fn resolve_genomes_path(
        sim_path: impl AsRef<Path>,
        time_step: u32,
        pond_i: u32
    ) -> PathBuf {
        sim_path.as_ref()
            .join(format!("{POND_PREFIX}{pond_i}"))
            .join(GENOMES_PATH)
            .join(format!("{time_step}.json"))
    }

    pub fn resolve_lattice_path(
        sim_path: impl AsRef<Path>,
        time_step: u32,
        pond_i: u32
    ) -> PathBuf {
        sim_path.as_ref()
            .join(format!("{POND_PREFIX}{pond_i}"))
            .join(LATTICES_PATH)
            .join(format!("{time_step}.txt"))
    }

    fn read_genomes(file_path: impl AsRef<Path>) -> anyhow::Result<Vec<SpinNodeLink>> {
        let file_path = file_path.as_ref();
        let file = std::fs::File::open(file_path).context(format!("while opening {}", file_path.display()))?;
        let reader = BufReader::new(file);
        Ok(serde_json::from_reader(reader)?)
    }

    pub fn read_lattice(file_path: impl AsRef<Path>, rect: Rect<usize>) -> anyhow::Result<Lattice<Spin>> {
        let file_path = file_path.as_ref();
        let file = std::fs::File::open(file_path).context(format!("while opening {}", file_path.display()))?;
        let buffer = BufReader::new(file);
        let mut numbers = Vec::new();
        let mut current = Vec::new();

        for byte_result in buffer.bytes() {
            let b = byte_result?;
            if b == b' ' {
                if !current.is_empty() {
                    let num = std::str::from_utf8(&current)?
                        .parse::<Spin>()?;
                    numbers.push(num);
                    current.clear();
                }
            } else {
                current.push(b);
            }
        }

        // Handle the last number if no trailing space
        if !current.is_empty() {
            let num = std::str::from_utf8(&current)?
                .parse::<Spin>()?;
            numbers.push(num);
        }

        if rect.area() != numbers.len() {
            bail!(
                "mismatch between `rect` area ({}) and size of the lattice stored in `file_path` ({})",
                rect.area(),
                numbers.len(),
            );
        }
        Ok(Lattice::from_box(
            numbers.into_boxed_slice(),
            rect
        ).unwrap())
    }

    pub fn write_if_time(
        &mut self,
        time_step: u32,
        ponds: &[Pond]
    ) -> anyhow::Result<()> {
        self.write_data_if_time(time_step, ponds)?;
        self.write_image_if_time(time_step, ponds)
    }

    fn write_data_if_time(
        &self,
        time_step: u32,
        ponds: &[Pond]
    ) -> anyhow::Result<()> {
        let time_str = time_step.to_string();
        if time_step % self.cell_period == 0 {
            for (i, pond) in ponds.iter().enumerate() {
                let mut celldf = pond.env.cells.to_dataframe()?;
                let file_path = self.outdir
                    .join(format!("pond_{i}"))
                    .join(CELLS_PATH)
                    .join(format!("{time_str}.parquet"));
                let file = std::fs::File::create(file_path)?;
                ParquetWriter::new(file).finish(&mut celldf)?;
            }
        }

        if time_step % self.genome_period == 0 {
            for (i, pond) in ponds.iter().enumerate() {
                let genomes = Self::make_genome_node_links(pond.env.cells());
                let file_path = self.outdir
                    .join(format!("pond_{i}"))
                    .join(GENOMES_PATH)
                    .join(format!("{time_str}.json"));
                let file = std::fs::File::create(file_path)?;
                serde_json::to_writer(file, &genomes)?;
            }
        }

        if time_step % self.lattice_period == 0 {
            for (i, pond) in ponds.iter().enumerate() {
                let file_path = self.outdir
                    .join(format!("pond_{i}"))
                    .join(LATTICES_PATH)
                    .join(format!("{time_str}.txt"));
                Self::write_lattice(file_path.as_path(), &pond.env.cell_lattice)?;
            }
        }
        Ok(())
    }

    fn make_genome_node_links(
        cells: &CellContainer<Cell>
    ) -> Vec<SpinNodeLink> {
        cells.iter()
            .map(|cell| SpinNodeLink {
                spin: cell.spin,
                node_link: NodeLinkData::from(cell.genome.clone())
            })
            .collect()
    }

    fn write_lattice<T: Display>(file_path: &Path, lattice: &Lattice<T>) -> io::Result<()>{
        let file = std::fs::File::create(file_path)?;
        let mut writer = BufWriter::new(file);
        let mut lat_it = lattice.iter_values();
        if let Some(val) = lat_it.next() {
            write!(writer, "{val}")?;
        }
        for val in lat_it {
            writer.write_all(b" ")?;
            write!(writer, "{val}")?;
        }
        writer.flush()
    }

    fn write_image_if_time(
        &mut self,
        time_step: u32,
        ponds: &[Pond],
    ) -> anyhow::Result<()> {
        // There might be a way to use LazyCell here but i got tired of fighting the borrow checker
        let mut frame = None;
        let movie_update = if let Some(mm) = &self.movie_maker {
            time_step % mm.frame_period == 0 && mm.window_works()
        } else {
            false
        };
        if movie_update {
            frame = Some(self.make_simulation_image(ponds)?);
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
                frame = Some(self.make_simulation_image(ponds)?);
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
        ponds: &[Pond]
    ) -> Result<RgbaImage, HexError> {
        let mut images = vec![];
        for pond in ponds {
            let env = &pond.env;
            images.push(RgbaImage::new(
                env.width() as u32,
                env.height() as u32
            ));
            let image = images.last_mut().unwrap();
            for plot in self.plots.order.clone() {
                match plot {
                    PlotType::Spin => {
                        SpinPlot {
                            env,
                            solid_color: hex_to_srgb(&self.plots.solid_color)?,
                            medium_color: match &self.plots.medium_color {
                                None => None,
                                Some(c) => Some(hex_to_srgb(c)?)
                            }
                        }.plot(image)
                    },
                    PlotType::Center => {
                        CenterPlot {
                            env,
                            color: hex_to_srgb(&self.plots.center_color)?
                        }.plot(image)
                    },
                    PlotType::ChemCenter => {
                        ChemCenterPlot {
                            env,
                            color: hex_to_srgb(&self.plots.chem_center_color)?
                        }.plot(image)
                    },
                    PlotType::Clones => {
                        ClonesPlot {
                            env,
                            clone_pairs: &pond.ca.adhesion.clone_pairs,
                            color: hex_to_srgb(&self.plots.clones_color)?,
                            all_clones: self.plots.all_clones
                        }.plot(image)
                    },
                    PlotType::Area => {
                        AreaPlot{
                            env,
                            min_color: srgb_to_luv(hex_to_srgb(&self.plots.area_min_color)?),
                            max_color: srgb_to_luv(hex_to_srgb(&self.plots.area_max_color)?),
                        }.plot(image)
                    },
                    PlotType::Border => {
                        BorderPlot {
                            env,
                            color: hex_to_srgb(&self.plots.border_color)?
                        }.plot(image)
                    },
                    PlotType::Chem => {
                        ChemPlot {
                            lat: &env.chem_lattice,
                            min_color: srgb_to_luv(hex_to_srgb(&self.plots.chem_min_color)?),
                            max_color: srgb_to_luv(hex_to_srgb(&self.plots.chem_max_color)?)
                        }.plot(image)
                    },
                    PlotType::CellType => {
                        CellTypePlot {
                            env,
                            mig_color: hex_to_srgb(&self.plots.migrating_color)?,
                            div_color: hex_to_srgb(&self.plots.dividing_color)?,
                        }.plot(image)
                    }
                }
            }
        }
        let mut grid = Self::grid_layout(&images).expect("must have at least one pond");
        flip_vertical_in_place(&mut grid);
        Ok(grid)
    }

    fn grid_layout(images: &[RgbaImage]) -> Option<RgbaImage> {
        if images.is_empty() {
            return None;
        }
        let img_width = images[0].width();
        let img_height = images[0].height();
        let n_images = images.len() as u32;
        let cols = (n_images as f32).sqrt().ceil() as u32;
        let rows = n_images.div_ceil(cols);

        let mut result = RgbaImage::new(
            img_width * cols,
            img_height * rows
        );
        for (i, img) in images.iter().enumerate() {
            let i_u32 = i as u32;
            let col = i_u32 % cols;
            let row = i_u32 / cols;
            let x = col * img_width;
            let y = row * img_height;
            result.copy_from(img, x, y).unwrap();
        }

        Some(result)
    }
}

trait ToDataFrame {
    fn to_dataframe(&self) -> PolarsResult<DataFrame>;
}

impl ToDataFrame for CellContainer<Cell> {
    fn to_dataframe(&self) -> PolarsResult<DataFrame> {
        let valid = self.iter().filter(|cell| cell.is_valid()).collect::<Vec<_>>();
        df!(
            "spin" => valid.iter().map(|cell| cell.spin).collect::<Vec<_>>(),
            "mom" => valid.iter().map(|cell| cell.mom).collect::<Vec<_>>(),
            "area" => valid.iter().map(|cell| cell.area()).collect::<Vec<_>>(),
            "target_area" => valid.iter().map(|cell| cell.target_area()).collect::<Vec<_>>(),
            "newborn_target_area" => valid.iter().map(|cell| cell.newborn_target_area).collect::<Vec<_>>(),
            "divide_area" => valid.iter().map(|cell| cell.divide_area).collect::<Vec<_>>(),
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
struct SpinNodeLink {
    spin: Spin,
    node_link: NodeLinkData<GrnGeneType, EdgeWeight, GrnMutParams>
}
