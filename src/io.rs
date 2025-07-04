use std::hash::{DefaultHasher, Hash, Hasher};
use std::path::Path;
use image::{EncodableLayout, RgbImage};
use std::io;
use minifb::{Window, WindowOptions};
use crate::boundary::Boundary;
use crate::cell::Cell;
use crate::constants::Spin;
use crate::environment::{Environment, LatticeEntity};
use crate::pos::Pos2D;

pub(crate) static IMAGES_PATH: &str = "images";
pub(crate) static CONFIG_COPY_PATH: &str = "config.toml";

pub fn create_directories(outpath: impl AsRef<Path>, replace_outdir: bool) -> io::Result<()> {
    let outpath = outpath.as_ref();
    let outdir_exists = outpath.try_exists()?;
    if outdir_exists {
        if replace_outdir {
            log::info!("Cleaning contents of '{}'", outpath.display());
            std::fs::remove_dir_all(outpath)?;
        } else {
            return Err(io::Error::new(
                io::ErrorKind::AlreadyExists,
                "`outdir` already exists and `replace_outdir` is `false`"
            ));
        }
    }
    std::fs::create_dir_all(outpath)?;
    std::fs::create_dir(outpath.join(IMAGES_PATH))
}

pub fn simulation_image(env: &Environment) -> RgbImage {
    let spins: Vec<_> = env
        .cell_lattice
        .iter_values()
        .flat_map(spin_to_rgb)
        .collect();
    
    let mut image = RgbImage::from_vec(
        env.width() as u32,
        env.height() as u32,
        spins
    ).unwrap();
    
    for cell in &env.cell_vec {
        let center = env.cell_lattice.bound.valid_pos(Pos2D::new(
            cell.center.pos().x as isize,
            cell.center.pos().y as isize,
        ));
        if let Some(pos) = center {
            image.put_pixel(pos.x as u32, pos.y as u32, [0, 255, 0].into());
        }
    }
    image
}

/// Converts a spin into a unique color.
///
/// This method guarantees 5232 unique colors, starting from this spin the colors will repeat.
fn spin_to_rgb(spin: Spin) -> [u8; 3] {
    if spin == LatticeEntity::Medium::<&Cell>.spin() {
        return [255, 255, 255];
    } else if spin == LatticeEntity::Solid::<&Cell>.spin() {
        return [0, 0, 0]
    }

    let mut hasher = DefaultHasher::new();
    spin.hash(&mut hasher);
    let hashed = hasher.finish();
    [
        (hashed & 0xFF).try_into().unwrap(),
        (hashed >> 8 & 0xFF) as u8,
        (hashed >> 16 & 0xFF) as u8,
    ]
}

pub struct MovieMaker {
    pub width: u32,
    pub height: u32,
    pub window: Window,
}

impl MovieMaker {
    pub fn new(width: u32, height: u32) -> minifb::Result<Self> {
        let window = Window::new(
            "Evo-CPM",
            width as usize,
            height as usize,
            WindowOptions::default()
        )?;
        Ok(Self {
            width,
            height,
            window
        })
    }
    
    pub fn window_works(&self) -> bool {
        self.window.is_open() && !self.window.is_key_down(minifb::Key::Escape)
    }

    pub fn update(&mut self, image: &RgbImage) -> minifb::Result<()> {
        let buffer: Vec<_> = image
            .as_bytes()
            .chunks_exact(3)
            .map(|rgb| {
                u32::from_le_bytes([rgb[0], rgb[1], rgb[2], 255])
            })
            .collect();
        self.window.update_with_buffer(&buffer, self.width as usize, self.height as usize)
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;
    use super::*;

    #[test]
    fn test_spin_to_rgb() {
        let mut tested = HashSet::<[u8; 3]>::default();
        // We can guarantee 5232 unique colors with this method, after that colors repeat
        for i in 0..5232 as Spin {
            let rgb = spin_to_rgb(i);
            assert!(!tested.contains(&rgb));
            tested.insert(rgb);
        }
    }
}
