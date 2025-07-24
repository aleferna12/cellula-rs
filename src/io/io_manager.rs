use crate::environment::Environment;
use crate::io::movie_maker::MovieMaker;
use crate::io::parameters::{IoParameters, MovieParameters, Parameters, PlotParameters, PlotType};
use crate::io::plot::{hex_to_srgb, srgb_to_luv, AreaPlot, BorderPlot, CellTypePlot, CenterPlot, ClonesPlot, LightCenterPlot, LightPlot, Plot, SpinPlot};
use crate::spin_table::SpinTable;
use image::imageops::flip_vertical_in_place;
use image::RgbaImage;
use std::error::Error;
use std::path::PathBuf;

pub(crate) static IMAGES_PATH: &str = "images";
pub(crate) static CONFIG_COPY_PATH: &str = "config.toml";

pub struct IoManager {
    pub outdir: PathBuf,
    pub replace_outdir: bool,
    pub image_period: u32,
    pub image_format: String,
    pub movie_maker: Option<MovieMaker>,
    pub plots: PlotParameters
}

impl IoManager {
    pub fn image_io(
        &mut self,
        time_step: u32,
        env: &Environment,
        clone_pairs: &SpinTable<bool>
    ) -> Result<(), Box<dyn Error>> {
        let mut frame = None;
        let movie_update = if let Some(mm) = &self.movie_maker {
            time_step % mm.frame_period == 0 && mm.window_works()
        } else {
            false
        };
        if movie_update {
            frame = Some(self.simulation_image(env, clone_pairs)?);
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
                frame = Some(self.simulation_image(env, clone_pairs)?);
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
        env: &Environment, 
        clone_pairs: &SpinTable<bool>
    ) -> Result<RgbaImage, Box<dyn Error>> {
        let mut image = RgbaImage::new(
            env.width() as u32,
            env.height() as u32
        );
        for plot in self.plots.order.clone() {
            match plot {
                PlotType::Spin => {
                    SpinPlot {
                        env,
                        solid_color: hex_to_srgb(
                            &self.plots.solid_color
                        ).expect("`solid-color` is not a valid rgb"),
                        medium_color: match &self.plots.medium_color { 
                            None => None,
                            Some(c) => Some(hex_to_srgb(c)?)
                        }
                    }.plot(&mut image)
                },
                PlotType::Center => {
                    CenterPlot {
                        env,
                        color: hex_to_srgb(&self.plots.center_color)?
                    }.plot(&mut image)
                },
                PlotType::LightCenter => {
                    LightCenterPlot {
                        env,
                        color: hex_to_srgb(&self.plots.light_center_color)?
                    }.plot(&mut image)
                },
                PlotType::Clones => {
                    ClonesPlot {
                        env,
                        clone_pairs,
                        color: hex_to_srgb(&self.plots.clones_color)?,
                        all_clones: self.plots.all_clones
                    }.plot(&mut image)
                },
                PlotType::Area => {
                    AreaPlot{ 
                        env,
                        min_color: srgb_to_luv(hex_to_srgb(&self.plots.area_min_color)?),
                        max_color: srgb_to_luv(hex_to_srgb(&self.plots.area_max_color)?),
                    }.plot(&mut image)
                },
                PlotType::Border => {
                    BorderPlot {
                        env,
                        color: hex_to_srgb(&self.plots.border_color)?
                    }.plot(&mut image)
                },
                PlotType::Light => {
                    LightPlot {
                        env,
                        min_color: srgb_to_luv(hex_to_srgb(&self.plots.light_min_color)?),
                        max_color: srgb_to_luv(hex_to_srgb(&self.plots.light_max_color)?)
                    }.plot(&mut image)
                },
                PlotType::CellType => {
                    CellTypePlot {
                        env,
                        mig_color: hex_to_srgb(&self.plots.migrating_color)?,
                        div_color: hex_to_srgb(&self.plots.dividing_color)?,
                    }.plot(&mut image)
                }
            }
        }
        flip_vertical_in_place(&mut image);
        Ok(image)
    }

    pub fn create_directories(&self) -> std::io::Result<()> {
        let outdir_exists = self.outdir.try_exists()?;
        if outdir_exists {
            if self.replace_outdir {
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

impl TryFrom<IoParameters> for IoManager {
    type Error = <MovieMaker as TryFrom<MovieParameters>>::Error;
    fn try_from(params: IoParameters) -> Result<Self, Self::Error> {
        Ok(Self {
            outdir: PathBuf::from(params.outdir),
            replace_outdir: params.replace_outdir,
            image_period: params.image_period,
            image_format: params.image_format,
            movie_maker: if params.movie.show {
                match MovieMaker::try_from(params.movie.clone()) {
                    Ok(mm) => {
                        log::info!("Creating window for real-time movie display");
                        Some(mm)
                    },
                    Err(e) => {
                        log::warn!("Failed to initialise movie maker with error `{e}`");
                        None
                    }
                }
            } else {
                None
            },
            plots: params.plots
        })
    }
}