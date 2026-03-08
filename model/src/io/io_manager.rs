//! Contains logic related to [`IoManager`].

use crate::io::parameters::Parameters;
use crate::my_cell::MyCell;
use crate::my_environment::MyEnvironment;
use anyhow::Context;
use bon::Builder;
use cellulars::io::read::parquet_reader::ParquetReader;
use cellulars::io::read::read_trait::Read;
#[cfg(feature = "movie-io")]
use cellulars::io::write::image::movie_window::MovieWindow;
use cellulars::io::write::image::plot::Plot;
use cellulars::io::write::parquet_writer::ParquetWriter;
use cellulars::io::write::write_trait::Write;
use cellulars::lattice::Lattice;
use cellulars::prelude::CellContainer;
use cellulars::spin::Spin;
use image::imageops::{flip_vertical_in_place, FilterType};
use image::{ColorType, GrayImage, ImageReader, RgbaImage};
use std::collections::HashSet;
use std::fs::File;
use std::io;
use std::path::{Path, PathBuf};

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
    /// Period with which to save an image of the simulation.
    pub image_period: u32,
    /// Period with which to save cell data.
    pub cells_period: u32,
    /// Period with which to save the cell lattice.
    pub lattice_period: u32,
    /// Used to update the simulation video when it's time.
    #[cfg(feature = "movie-io")]
    pub movie_module: Option<MovieModule>,
    plots: Box<[Box<dyn Plot<MyEnvironment>>]>,
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

    fn pad_time_step(time_step: u32) -> String {
        format!("{time_step:0>PAD_FILE_LEN$}")
    }

    /// Given a path to the main folder of a simulation, resolve the path to the file
    /// containing the simulation parameters.
    fn parameters_path(sim_path: impl AsRef<Path>) -> PathBuf {
        sim_path.as_ref().join(CONFIG_COPY_PATH)
    }

    /// Given a path to the main folder of a simulation, resolve the path to the cell data file
    /// that was saved at `time_step`.
    fn cells_path(
        sim_path: impl AsRef<Path>,
        time_step: u32
    ) -> PathBuf {
        sim_path.as_ref()
            .join(CELLS_PATH)
            .join(format!("{}.parquet", Self::pad_time_step(time_step)))
    }

    /// Reads the parameter file of a simulation at `sim_path`.
    pub fn read_parameters(sim_path: impl AsRef<Path>) -> anyhow::Result<Parameters> {
        Parameters::parse(IoManager::parameters_path(sim_path))
    }

    /// Reads a cell [`CellContainer`] from the simulation at `sim_path` at `time_step`.
    pub fn read_cells(sim_path: impl AsRef<Path>, time_step: u32) -> anyhow::Result<CellContainer<MyCell>> {
        ParquetReader { reader: File::open(Self::cells_path(sim_path, time_step))? }
            .read()
            .map_err(anyhow::Error::new)
    }

    /// Reads a cell [`Lattice`] from the simulation at `sim_path` at `time_step`.
    pub fn read_cell_lattice(sim_path: impl AsRef<Path>, time_step: u32) -> anyhow::Result<Lattice<Spin>> {
        ParquetReader{ reader: File::open(Self::lattice_path(sim_path, time_step))? }
            .read()
            .map_err(anyhow::Error::new)
    }

    /// Reads a layout PNG file at `layout_path` for a pond with dimensions
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
        Ok(layout.resize_exact(pond_width as u32, pond_height as u32, FilterType::Nearest).flipv().into_luma8())
    }

    /// Given a path to the main folder of a simulation, resolve the path to the lattice file
    /// that was saved at `time_step`.
    fn lattice_path(
        sim_path: impl AsRef<Path>,
        time_step: u32
    ) -> PathBuf {
        sim_path.as_ref()
            .join(LATTICES_PATH)
            .join(format!("{}.parquet", Self::pad_time_step(time_step)))
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
            let file_path = self.outdir
                .join(CELLS_PATH)
                .join(format!("{time_str}.parquet"));
            ParquetWriter { writer: File::create(file_path)?, overwrites: vec![] }.write(&env.env.cells)?;
        }

        if time_step.is_multiple_of(self.lattice_period) {
            let file_path = self.outdir
                .join(LATTICES_PATH)
                .join(format!("{time_str}.parquet"));
            ParquetWriter { writer: File::create(file_path)?, overwrites: vec![] }.write(&env.env.cell_lattice)?;
        }
        Ok(())
    }

    fn write_image_if_time(
        &mut self,
        time_step: u32, 
        env: &MyEnvironment
    ) -> anyhow::Result<()> {
        // This looks like it should be a LazyCell but that doesnt work (i tried)
        let mut frame = None;

        #[cfg(feature = "movie-io")]
        let movie_update = if let Some(mm) = &self.movie_module {
            time_step.is_multiple_of(mm.frame_period) && mm.movie_window.is_open()
        } else {
            false
        };
        #[cfg(feature = "movie-io")]
        if movie_update {
            frame = Some(self.make_simulation_image(env));
            let mm = self.movie_module.as_mut().unwrap();
            let resized = image::imageops::resize(
                frame.as_ref().unwrap(),
                mm.movie_window.width.try_into()?,
                mm.movie_window.height.try_into()?,
                image::imageops::Nearest,
            );
            mm.movie_window.update(&resized)?
        }

        if time_step.is_multiple_of(self.image_period) {
            if frame.is_none() {
                frame = Some(self.make_simulation_image(env));
            }
            frame.unwrap().save(
                &self.outdir
                    .join(IMAGES_PATH)
                    .join(format!(
                        "{}.webp",
                        Self::pad_time_step(time_step),
                    ))
            )?;
        }
        Ok(())
    }

    /// Makes a new frame of the simulation by drawing a succession of plots.
    pub fn make_simulation_image(
        &self, 
        env: &MyEnvironment
    ) -> RgbaImage {
        let mut image = RgbaImage::new(
            env.env.width() as u32,
            env.env.height() as u32
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

#[cfg(feature = "movie-io")]
/// Groups together a movie window and its associated frame rate.
pub struct MovieModule {
    /// Movie window used to display the simulation.
    pub movie_window: MovieWindow,
    /// Period with which to update the movie.
    pub frame_period: u32
}
