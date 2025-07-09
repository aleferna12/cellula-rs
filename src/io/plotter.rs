use std::collections::HashSet;
use std::hash::{DefaultHasher, Hash, Hasher};
use image::RgbImage;
use imageproc::drawing::{draw_cross_mut, draw_line_segment_mut};
use crate::cell::RelCell;
use crate::constants::Spin;
use crate::environment::{Environment, LatticeEntity};
use crate::positional::boundary::Boundary;
use crate::positional::pos::Pos2D;

pub trait Plotter {
    fn plot(&self, image: &mut RgbImage, env: &Environment);
}

pub struct SpinPlotter {}

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

pub struct CellCenterPlotter {}

impl Plotter for CellCenterPlotter {
    fn plot(&self, image: &mut RgbImage, env: &Environment) {
        for cell in &env.cells {
            let center = env.cell_lattice.bound.valid_pos(Pos2D::new(
                cell.center.pos().x as isize,
                cell.center.pos().y as isize,
            ));
            if let Some(pos) = center {
                draw_cross_mut(image, [0, 255, 0].into(), pos.x as i32, pos.y as i32);
            }
        }
    }
}

pub struct ClonesPlotter<'a> {
    pub(crate) clone_pairs: &'a HashSet<(Spin, Spin)>
}

impl Plotter for ClonesPlotter<'_> {
    fn plot(&self, image: &mut RgbImage, env: &Environment) {
        for (spin1, spin2) in self.clone_pairs.iter().copied() {
            let message = "non-cell stored as clone";
            let center1 = env.cells.get_entity(spin1).expect_cell(message).center.pos;
            let center2 = env.cells.get_entity(spin2).expect_cell(message).center.pos;
            draw_line_segment_mut(
                image,
                (center1.x, center1.y),
                (center2.x, center2.y),
                [255, 0, 0].into()
            )
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