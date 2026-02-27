use std::io;
use crate::io::write::write::Write;
use image::{ImageError, ImageFormat, RgbaImage};

/// WEBP writer.
#[derive(Clone, Debug)]
pub struct WebpWriter<W> {
    file: W
}

impl<W> WebpWriter<W> {
    pub fn new(file: W) -> Self {
        Self { file }
    }
}

impl<W: io::Write + io::Seek> Write<RgbaImage, ImageError> for WebpWriter<W> {
    fn write(mut self, data: &RgbaImage) -> Result<(), ImageError> {
        data.write_to(&mut self.file, ImageFormat::WebP)
    }
}