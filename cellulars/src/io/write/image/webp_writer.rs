//! Contains logic associated with [`WebpWriter`].

use std::io;
use crate::io::write::r#trait::Write;
use image::{ImageError, ImageFormat, RgbaImage};

/// WEBP writer.
#[derive(Clone, Debug)]
pub struct WebpWriter<W> {
    /// Object responsible for writing the image.
    pub writer: W
}

impl<W: io::Write + io::Seek> Write<RgbaImage, ImageError> for WebpWriter<W> {
    fn write(mut self, data: &RgbaImage) -> Result<(), ImageError> {
        data.write_to(&mut self.writer, ImageFormat::WebP)
    }
}