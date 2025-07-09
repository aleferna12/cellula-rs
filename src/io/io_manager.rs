use std::collections::HashSet;
use std::path::PathBuf;
use image::RgbImage;
use crate::constants::Spin;
use crate::environment::Environment;
use crate::io::movie_maker::MovieMaker;
use crate::io::plotter::{CellCenterPlotter, ClonesPlotter, Plotter, SpinPlotter};
use crate::parameters::{IoParameters, MovieParameters, Parameters};

pub(crate) static IMAGES_PATH: &str = "images";
pub(crate) static CONFIG_COPY_PATH: &str = "config.toml";

pub struct IoManager {
    pub outdir: PathBuf,
    pub replace_outdir: bool,
    pub show_attached_cells: bool,
    pub show_cell_centers: bool,
    pub image_period: u32,
    pub image_format: String,
    pub movie_maker: Option<MovieMaker>
}

impl IoManager {
    pub fn image_io(
        &mut self,
        time_step: u32,
        env: &Environment,
        clone_pairs: &HashSet<(Spin, Spin)>
    ) {
        let mut frame = None;
        let movie_update = if let Some(mm) = &self.movie_maker {
            time_step % mm.frame_period == 0 && mm.window_works()
        } else {
            false
        };
        if movie_update {
            frame = Some(self.simulation_image(env, clone_pairs));
            let mm = self.movie_maker.as_mut().unwrap();
            let resized = image::imageops::resize(
                frame.as_ref().unwrap(),
                mm.width,
                mm.height,
                image::imageops::Nearest,
            );
            if let Err(e) = mm.update(&resized) {
                log::warn!(
                "Failed to display simulation frame at time step {} with error `{}`",
                time_step,
                e
            );
            }
        }

        if time_step % self.image_period == 0 {
            if frame.is_none() {
                frame = Some(self.simulation_image(env, clone_pairs));
            }
            let saved = frame.unwrap().save(
                &self.outdir
                    .join(IMAGES_PATH)
                    .join(format!("{time_step}.{}", self.image_format.to_lowercase()))
            );
            if let Err(e) = saved {
                log::warn!("Failed to save simulation frame at time step {} with error `{}`", time_step, e);
            }
        }
    }

    pub fn simulation_image(&self, env: &Environment, clone_pairs: &HashSet<(Spin, Spin)>) -> RgbImage {
        let mut image = RgbImage::from_pixel(
            env.width() as u32,
            env.height() as u32,
            [255, 255, 255].into()
        );
        SpinPlotter{}.plot(&mut image, env);
        if self.show_cell_centers {
            CellCenterPlotter{}.plot(&mut image, env);
        }
        if self.show_attached_cells {
            ClonesPlotter { clone_pairs }.plot(&mut image, env);
        }
        image
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
    
    pub fn create_parameters_file(&self, parameters: &Parameters) -> Result<(), Box<dyn std::error::Error>> {
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
            show_attached_cells: params.show_attached_cells,
            show_cell_centers: params.show_cell_centers,
            image_period: params.image_period,
            image_format: params.image_format,
            movie_maker: if params.movie.show {
                match MovieMaker::try_from(params.movie.clone()) {
                    Ok(mm) => {
                        log::info!("Creating window for real-time movie display");
                        Some(mm)
                    },
                    Err(e) => {
                        log::warn!("Failed to initialise movie maker with error `{}`", e);
                        None
                    }
                }
            } else {
                None
            },
        })
    }
}