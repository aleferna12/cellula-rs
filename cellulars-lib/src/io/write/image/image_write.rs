
use image::{ImageError, RgbaImage};
use std::path::Path;
use crate::io::write::writer::{Write, Writer};

impl Write<RgbaImage, ImageError> for Writer {
    fn write(&mut self, data: &RgbaImage, path: impl AsRef<Path>) -> Result<(), ImageError> {
        data.save(path).map(|_| ())
    }
}