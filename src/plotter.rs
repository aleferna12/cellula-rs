use std::hash::{DefaultHasher, Hash, Hasher};
use image::RgbImage;
use crate::cell::RelCell;
use crate::constants::Spin;
use crate::environment::{Environment, LatticeEntity};

pub trait Plotter {
    fn plot(&self, image: &mut RgbImage, env: &Environment);
}

struct SpinPlotter {

}

impl SpinPlotter {
    /// Converts a spin into a unique color.
    ///
    /// This method guarantees 5232 unique colors, starting from this spin the colors will repeat.
    fn spin_to_rgb(spin: Spin) -> [u8; 3] {
        if spin == LatticeEntity::Medium::<&RelCell>.spin() {
            return [255, 255, 255];
        } else if spin == LatticeEntity::Solid::<&RelCell>.spin() {
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
}

impl Plotter for SpinPlotter {
    fn plot(&self, image: &mut RgbImage, env: &Environment) {
        for pos in env.cell_lattice.iter_positions() {
            let spin = env.cell_lattice[pos];
            if spin < LatticeEntity::first_cell_spin() {
                continue;
            }
            image.put_pixel(pos.x as u32, pos.y as u32, Self::spin_to_rgb(spin).into());
        }
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;
    use crate::constants::Spin;
    use super::*;

    #[test]
    fn test_spin_to_rgb() {
        let mut tested = HashSet::<[u8; 3]>::default();
        // We can guarantee 5232 unique colors with this method, after that colors repeat
        for i in 0..5232 as Spin {
            let rgb = SpinPlotter::spin_to_rgb(i);
            assert!(!tested.contains(&rgb));
            tested.insert(rgb);
        }
    }
}