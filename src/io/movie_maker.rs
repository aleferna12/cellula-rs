use image::{EncodableLayout, RgbImage};
use crate::parameters::MovieParameters;

pub struct MovieMaker {
    pub width: u32,
    pub height: u32,
    pub frame_period: u32,
    pub window: minifb::Window,
}

impl MovieMaker {
    pub fn window_works(&self) -> bool {
        self.window.is_open() && !self.window.is_key_down(minifb::Key::Escape)
    }

    pub fn update(&mut self, image: &RgbImage) -> minifb::Result<()> {
        let buffer: Vec<_> = image
            .as_bytes()
            .chunks_exact(3)
            .map(|rgb| {
                u32::from_le_bytes([rgb[2], rgb[1], rgb[0], 255])
            })
            .collect();
        self.window.update_with_buffer(&buffer, self.width as usize, self.height as usize)
    }
}

impl TryFrom<MovieParameters> for MovieMaker {
    type Error = minifb::Error;

    fn try_from(params: MovieParameters) ->  Result<Self, Self::Error> {
        let window = minifb::Window::new(
            "Evo-CPM",
            params.width as usize,
            params.height as usize,
            minifb::WindowOptions::default()
        )?;
        Ok(Self {
            width: params.width,
            height: params.height,
            frame_period: params.frame_period,
            window
        })
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
        window.unwrap().update();
    }
}