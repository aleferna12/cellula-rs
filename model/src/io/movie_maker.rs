//! Contains logic to display the movie of the simulation.

use image::{EncodableLayout, RgbaImage};

/// This struct manages a [minifb::Window] to display simulation frames.
pub struct MovieMaker {
    /// Width of the screen in pixels.
    pub width: u32,
    /// Height of the screen in pixels.
    pub height: u32,
    /// Minimum period with which frames are updated.
    pub frame_period: u32,
    /// Window where the movie is displayed.
    pub window: minifb::Window,
}

impl MovieMaker {
    /// Make a new movie manager with specified parameters.
    pub fn new(
        width: u32,
        height: u32,
        frame_period: u32
    ) -> Result<Self, minifb::Error> {
        let window = minifb::Window::new(
            "Cellulars Model",
            width as usize,
            height as usize,
            minifb::WindowOptions::default()
        )?;
        Ok(Self {
            width,
            height,
            frame_period,
            window
        })
    }

    /// Whether the window is open and not being actively closed.
    pub fn window_works(&self) -> bool {
        self.window.is_open() && !self.window.is_key_down(minifb::Key::Escape)
    }

    /// Tries to update the window with a new frame `image`.
    pub fn update(&mut self, image: &RgbaImage) -> minifb::Result<()> {
        let buffer: Vec<_> = image
            .as_bytes()
            .chunks_exact(4)
            .map(|rgba| {
                u32::from_le_bytes([rgba[2], rgba[1], rgba[0], rgba[3]])
            })
            .collect();
        self.window.update_with_buffer(&buffer, self.width as usize, self.height as usize)
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn test_window() {
        let window = minifb::Window::new(
            "test",
            100,
            100,
            minifb::WindowOptions::default()
        );
        assert!(window.is_ok());
    }
}