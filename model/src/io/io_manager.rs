use crate::cell::Cell;
use crate::genetics::grn::{EdgeWeight, GrnGeneType};
use crate::io::movie_maker::MovieMaker;
use crate::io::node_link::{NodeLinkData, ToNodeLink};
use crate::io::parameters::{Parameters, PlotParameters, PlotType};
use crate::io::plot::*;
use crate::pond::Pond;
use cellulars_lib::cell_container::CellContainer;
use cellulars_lib::basic_cell::Cellular;
use cellulars_lib::constants::Spin;
use image::imageops::flip_vertical_in_place;
use image::{GenericImage, RgbaImage};
use polars::prelude::*;
use serde::Serialize;
use std::error::Error;
use std::path::{Path, PathBuf};

static IMAGES_PATH: &str = "images";
static CELLS_PATH: &str = "cells";
static GENOMES_PATH: &str = "genomes";
static LATTICES_PATH: &str = "lattices";
static CONFIG_COPY_PATH: &str = "config.toml";

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
    pub fn new(
        outdir: impl AsRef<Path>,
        image_format: String,
        image_period: u32,
        cell_period: u32,
        genome_period: u32,
        lattice_period: u32,
        plots: PlotParameters,
        movie_maker: Option<MovieMaker>,
    ) -> Self {
        Self {
            outdir: outdir.as_ref().to_path_buf(),
            image_format,
            plots,
            movie_maker,
            image_period,
            cell_period,
            genome_period,
            lattice_period,
        }
    }

    pub fn create_directories(&self, replace_outdir: bool) -> std::io::Result<()> {
        let outdir_exists = self.outdir.try_exists()?;
        if outdir_exists {
            if replace_outdir {
                log::info!("Cleaning contents of '{}'", self.outdir.display());
                std::fs::remove_dir_all(&self.outdir)?;
            } else {
                return Err(std::io::Error::new(
                    std::io::ErrorKind::AlreadyExists,
                    "`outdir` already exists and `replace_outdir` is `false`"
                ));
            }
        }
        std::fs::create_dir_all(&self.outdir)?;
        std::fs::create_dir(self.outdir.join(IMAGES_PATH))?;
        std::fs::create_dir(self.outdir.join(CELLS_PATH))?;
        std::fs::create_dir(self.outdir.join(GENOMES_PATH))?;
        std::fs::create_dir(self.outdir.join(LATTICES_PATH))
    }

    pub fn create_parameters_file(&self, parameters: &Parameters) -> Result<(), Box<dyn Error>> {
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

    pub fn try_io(
        &mut self,
        time_step: u32,
        ponds: &[Pond]
    ) -> Result<(), Box<dyn Error>> {
        self.try_data_io(time_step, ponds)?;
        self.try_image_io(time_step, ponds)
    }

    fn try_data_io(
        &self,
        time_step: u32,
        ponds: &[Pond]
    ) -> Result<(), Box<dyn Error>> {
        let time_str = time_step.to_string();
        if time_step % self.cell_period == 0 {
            let mut celldf = self.make_cell_dataframe(ponds)?.collect()?;
            let file_name = format!("{time_str}.parquet");
            let file = std::fs::File::create(self.outdir.join(CELLS_PATH).join(file_name))?;
            ParquetWriter::new(file).finish(&mut celldf)?;
        }

        if time_step % self.genome_period == 0 {
            let genomes = self.make_cell_node_links(ponds);
            let json = serde_json::to_string(&genomes)?;
            let file_name = format!("{time_str}.json");
            std::fs::write(self.outdir.join(GENOMES_PATH).join(file_name), json)?;
        }

        if time_step % self.lattice_period == 0 {
            let lattices = self.make_cell_lattices(ponds);
            let json = serde_json::to_string(&lattices)?;
            let file_name = format!("{time_str}.json");
            std::fs::write(self.outdir.join(LATTICES_PATH).join(file_name), json)?;
        }
        Ok(())
    }

    fn make_cell_dataframe(
        &self,
        ponds: &[Pond],
    ) -> PolarsResult<LazyFrame> {
        let mut dfs = vec![];
        for (i, pond) in ponds.iter().enumerate() {
            let cell_df = pond.env.cells.to_dataframe()?.lazy();
            dfs.push(cell_df.with_column(lit(i as u32).alias("pond")));
        }
        concat(dfs, UnionArgs::default())
    }

    fn make_cell_node_links(
        &self,
        ponds: &[Pond]
    ) -> Vec<NodeLinkData<GrnGeneType, EdgeWeight>> {
        let mut res = vec![];
        for (i, pond) in ponds.iter().enumerate() {
            for cell in pond.env.cells.iter() {
                let mut node_link = cell.genome.to_node_link();
                node_link.graph.insert("pond".to_string(), serde_json::json!(i));
                node_link.graph.insert("spin".to_string(), serde_json::json!(cell.spin));
                res.push(node_link);
            }
        }
        res
    }

    fn make_cell_lattices<'a>(&self, ponds: &'a [Pond]) -> Vec<PondCellLatttice<'a>> {
        let mut res = vec![];
        for (i, pond) in ponds.iter().enumerate() {
            res.push(PondCellLatttice {
                pond: i as u32,
                lattice: pond.env.cell_lattice.as_array()
            })
        }
        res
    }

    fn try_image_io(
        &mut self,
        time_step: u32,
        ponds: &[Pond],
    ) -> Result<(), Box<dyn Error>> {
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
            "divide_area" => valid.iter().map(|cell| cell.divide_area).collect::<Vec<_>>(),
            "center_x" => valid.iter().map(|cell| cell.center().x).collect::<Vec<_>>(),
            "center_y" => valid.iter().map(|cell| cell.center().y).collect::<Vec<_>>(),
            "chem_center_x" => valid.iter().map(|cell| cell.chem_center.x).collect::<Vec<_>>(),
            "chem_center_y" => valid.iter().map(|cell| cell.chem_center.y).collect::<Vec<_>>(),
            "chem_mass" => valid.iter().map(|cell| cell.chem_mass).collect::<Vec<_>>(),
            "cell_type" => valid.iter().map(|cell| cell.is_dividing()).collect::<Vec<_>>()
        )
    }
}

#[derive(Serialize)]
struct PondCellLatttice<'a> {
    pond: u32,
    lattice: &'a[Spin]
}
