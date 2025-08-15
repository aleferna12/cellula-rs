use crate::io::movie_maker::MovieMaker;
use crate::io::parameters::{Parameters, PlotParameters, PlotType};
use crate::io::plot::{hex_to_srgb, srgb_to_luv, AreaPlot, BorderPlot, CellTypePlot, CenterPlot, ChemCenterPlot, ChemPlot, ClonesPlot, HexError, MaternalFactorPlot, Plot, SpinPlot};
use crate::pond::Pond;
use image::imageops::flip_vertical_in_place;
use image::{GenericImage, RgbaImage};
use std::error::Error;
use std::path::{Path, PathBuf};

pub(crate) static IMAGES_PATH: &str = "images";
pub(crate) static CONFIG_COPY_PATH: &str = "config.toml";

pub struct IoManager {
    pub outdir: PathBuf,
    pub image_period: u32,
    pub image_format: String,
    pub movie_maker: Option<MovieMaker>,
    pub plots: PlotParameters
}

impl IoManager {
    pub fn new(
        outdir: impl AsRef<Path>,
        image_period: u32,
        image_format: String,
        plots: PlotParameters,
        movie_maker: Option<MovieMaker>
    ) -> Self {
        Self {
            outdir: outdir.as_ref().to_path_buf(),
            image_period,
            image_format,
            plots,
            movie_maker

        }
    }

    pub fn image_io(
        &mut self,
        time_step: u32,
        ponds: &Vec<Pond>,
    ) -> Result<(), Box<dyn Error>> {
        let mut frame = None;
        let movie_update = if let Some(mm) = &self.movie_maker {
            time_step % mm.frame_period == 0 && mm.window_works()
        } else {
            false
        };
        if movie_update {
            frame = Some(self.simulation_image(ponds)?);
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
                frame = Some(self.simulation_image(ponds)?);
            }
            frame.unwrap().save(
                &self.outdir
                    .join(IMAGES_PATH)
                    .join(format!("{time_step}.{}", self.image_format.to_lowercase())
            ))?;
        }
        Ok(())
    }

    pub fn simulation_image(
        &self, 
        ponds: &Vec<Pond>
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
                            space: &env.space,
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
                            lat: &env.space.chem_lattice,
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
                    PlotType::MaternalFactor => {
                        MaternalFactorPlot {
                            env,
                        }.plot(image)
                    }
                }
            }
        }
        let mut grid = Self::grid_layout(&images).expect("must have at least one pond");
        flip_vertical_in_place(&mut grid);
        Ok(grid)
    }

    pub fn grid_layout(images: &[RgbaImage]) -> Option<RgbaImage> {
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
        std::fs::create_dir(self.outdir.join(IMAGES_PATH))
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
}