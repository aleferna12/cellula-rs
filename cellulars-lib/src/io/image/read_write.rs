use crate::io::file::file_path;
use crate::io::writer::{Write, Writer};
use image::{ImageError, RgbaImage};
use std::path::PathBuf;

impl Write<RgbaImage, ImageError> for Writer {
    fn write(&mut self, data: &RgbaImage, time_step: u32) -> Result<PathBuf, ImageError> {
        let file_path = file_path(
            self.outdir.as_path(),
            "images",
            "webp",
            time_step
        ).expect("failed to pad time step when saving image");  // This should never fail
        data.save(&file_path).map(|_| file_path)
    }
}