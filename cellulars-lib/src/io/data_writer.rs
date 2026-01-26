use std::path::PathBuf;
use image::{ImageError, RgbaImage};
use crate::io::file::{pad_file_path, U32_STR_LEN};

pub trait WriteData<D, E> {
    fn write(&mut self, data: &mut D, time_step: u32) -> Result<PathBuf, E>;
}

pub struct DataWriter {
    outdir: PathBuf
}

impl WriteData<RgbaImage, ImageError> for DataWriter {
    fn write(&mut self, data: &mut RgbaImage, time_step: u32) -> Result<PathBuf, ImageError> {
        let file_path = pad_file_path(
            &self.outdir.join("images").join(format!("{time_step}.webp")),
            U32_STR_LEN
        ).expect("failed to pad time step when saving image");  // This should never fail
        data.save(&file_path).map(|_| file_path)
    }
}