//! Contains logic associated with [`MovieWindow`].

use image::{EncodableLayout, RgbaImage};
use minifb::{Window, WindowOptions};

/// Movie window used to display the simulation in real time.
pub struct MovieWindow {
    /// Width of the window.
    pub width: usize,
    /// Height of the window.
    pub height: usize,
    /// Backend [`Window`] where the movie is displayed.
    window: Window,
}

impl MovieWindow {
    /// Makes a new movie window.
    pub fn new(
        width: usize,
        height: usize,
    ) -> Result<Self, minifb::Error> {
        let window = Window::new(
            "Cellulars Model",
            width,
            height,
            WindowOptions::default()
        )?;
        Ok(Self {
            width,
            height,
            window
        })
    }

    /// Tries to update the window with a new frame `image`.
    pub fn update(&mut self, image: &RgbaImage) -> minifb::Result<()> {
        let buffer: Box<_> = image
            .as_bytes()
            .chunks_exact(4)
            .map(|rgba| {
                u32::from_le_bytes([rgba[2], rgba[1], rgba[0], rgba[3]])
            })
            .collect();
        self.window.update_with_buffer(&buffer, self.width, self.height)
    }
    
    /// Returns whether the window is open.
    pub fn is_open(&self) -> bool {
        self.window.is_open() && !self.window.is_key_down(minifb::Key::Escape)
    }
}