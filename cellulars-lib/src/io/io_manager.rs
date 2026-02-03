use crate::io::pad_file::{pad_file_name, U32_STR_LEN};
use std::ffi::OsString;
use std::path::PathBuf;

pub struct IoManager {
    root_dir: PathBuf,
    pad_files: bool
}

impl IoManager {
    fn file_path(
        &self,
        subfolder: &str,
        ext: &str,
        time_step: u32
    ) -> Option<PathBuf> {
        let stem = format!("{time_step}.{ext}");
        let file_name = if self.pad_files {
            pad_file_name(
                &stem,
                U32_STR_LEN
            )?
        } else {
            OsString::from(stem)
        };
        Some(self.root_dir.join(subfolder).join(file_name))
    }
}

#[cfg(feature = "image-io")]
mod images {
    use crate::io::io_manager::IoManager;
    use crate::io::writer::{Write, Writer};
    use image::{ImageError, RgbaImage};
    use std::path::PathBuf;

    impl IoManager {
        pub fn write_images(&self, image: &RgbaImage, time_step: u32) -> Result<PathBuf, ImageError> {
            let file_path = self.file_path(
                "images",
                "webp",
                time_step
            ).expect("failed to pad time step when saving image");  // This should never fail
            Writer {}.write(image, &file_path).map(|_| file_path)
        }
    }
}

#[cfg(feature = "data-io")]
mod data {
    use crate::io::io_manager::IoManager;
    use crate::io::writer::data::CellsWriteError;
    use crate::io::writer::{Write, Writer};
    use crate::prelude::{CellContainer, Cellular, Lattice, Spin};
    use parquet::errors::ParquetError;
    use serde::{Deserialize, Serialize};
    use std::path::PathBuf;

    impl IoManager {
        pub fn write_lattice(
            &self,
            lattice: &Lattice<Spin>,
            time_step: u32
        ) -> Result<PathBuf, ParquetError> {
            let file_path = self.file_path(
                "lattices",
                "parquet",
                time_step
            ).expect("failed to pad time step when saving cell lattice");  // This should never fail
            Writer {}.write(lattice, &file_path).map(|_| file_path)
        }

        pub fn write_cells<'de, C: Cellular + Serialize + Deserialize<'de>>(
            &self,
            cells: &CellContainer<C>,
            time_step: u32
        ) -> Result<PathBuf, CellsWriteError> {
            let file_path = self.file_path(
                "cells",
                "parquet",
                time_step
            ).expect("failed to pad time step when saving cells");  // This should never fail
            Writer {}.write(cells, &file_path).map(|_| file_path)
        }
    }
}