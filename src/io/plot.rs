use image::{Pixel, Rgb, RgbaImage};
use imageproc::drawing::{draw_cross_mut, draw_line_segment_mut};
use std::collections::HashSet;
use std::error::Error;
use std::hash::{DefaultHasher, Hash, Hasher};
use crate::constants::Spin;
use crate::environment::{Environment, LatticeEntity};
use crate::positional::boundary::Boundary;
use crate::positional::pos::Pos2D;

pub trait Plot {
    fn plot(&self, image: &mut RgbaImage);
}

pub struct SpinPlot<'e> {
    env: &'e Environment,
    solid_color: Rgb<u8>,
    medium_color: Option<Rgb<u8>>
}

impl<'e> SpinPlot<'e> {
    pub fn new(env: &'e Environment, solid_color: Rgb<u8>, medium_color: Option<Rgb<u8>>) -> Self {
        Self { env, solid_color, medium_color}
    }

    fn spin_to_rgb(spin: Spin) -> Rgb<u8> {
        let mut hasher = DefaultHasher::new();
        spin.hash(&mut hasher);
        let hashed = hasher.finish();
        [
            (hashed & 0xFF).try_into().unwrap(),
            (hashed >> 8 & 0xFF) as u8,
            (hashed >> 16 & 0xFF) as u8
        ].into()
    }
}

impl Plot for SpinPlot<'_> {
    fn plot(&self, image: &mut RgbaImage) {
        for pos in self.env.cell_lattice.iter_positions() {
            let spin = self.env.cell_lattice[pos];
            let rgb = if spin > LatticeEntity::first_cell_spin() {
                Some(Self::spin_to_rgb(spin))
            } else if spin == LatticeEntity::Solid.spin() {
                Some(self.solid_color)
            } else {
                self.medium_color
            };
            if let Some(color) = rgb {
                image.put_pixel(pos.x as u32, pos.y as u32, color.to_rgba());
            }
        }
    }
}

pub struct CenterPlot<'e> {
    env: &'e Environment,
    color: Rgb<u8>
}

impl<'e> CenterPlot<'e> {
    pub fn new(env: &'e Environment, color: Rgb<u8>) -> Self {
        Self { env, color }
    }
}

impl Plot for CenterPlot<'_> {
    fn plot(&self, image: &mut RgbaImage) {
        for cell in &self.env.cells {
            let center = self.env.cell_lattice.bound.valid_pos(Pos2D::new(
                cell.center.pos().x as isize,
                cell.center.pos().y as isize,
            ));
            if let Some(pos) = center {
                draw_cross_mut(image, self.color.to_rgba(), pos.x as i32, pos.y as i32);
            }
        }
    }
}

pub struct ClonesPlot<'a> {
    env: &'a Environment,
    clone_pairs: &'a HashSet<(Spin, Spin)>,
    color: Rgb<u8>
}

impl<'a> ClonesPlot<'a> {
    pub fn new(env: &'a Environment, clone_pairs: &'a HashSet<(Spin, Spin)>, color: Rgb<u8>) -> Self {
        Self { env, clone_pairs, color }
    }
}

impl Plot for ClonesPlot<'_> {
    fn plot(&self, image: &mut RgbaImage) {
        for (spin1, spin2) in self.clone_pairs.iter().copied() {
            let message = "non-cell stored as clone";
            let center1 = self.env.cells.get_entity(spin1).expect_cell(message).center.pos;
            let center2 = self.env.cells.get_entity(spin2).expect_cell(message).center.pos;
            draw_line_segment_mut(
                image,
                (center1.x, center1.y),
                (center2.x, center2.y),
                self.color.to_rgba()
            )
        }
    }
}

pub struct AreaPlot<'e> {
    env: &'e Environment
}

impl<'e> AreaPlot<'e> {
    pub fn new(env: &'e Environment) -> Self {
        Self { env }
    }
}

impl Plot for AreaPlot<'_> {
    fn plot(&self, image: &mut RgbaImage) {
        let mut min = u32::MAX;
        let mut max = 0;
        for cell in &self.env.cells {
            if cell.area < min {
                min = cell.area
            }
            if cell.area > max {
                max = cell.area
            }
        }
        for pos in self.env.cell_lattice.iter_positions() {
            let entity = self.env.cells.get_entity(self.env.cell_lattice[pos]);
            if let LatticeEntity::SomeCell(cell) = entity {
                let frac = lerp(cell.area as f32, min as f32, max as f32);
                let gray = (255. * frac) as u8;
                image.put_pixel(pos.x as u32, pos.y as u32, [gray, gray, gray, 255].into());
            }
        }
    }
}

pub fn lerp(value: f32, min: f32, max: f32) -> f32 {
    let num = value - min;
    if num == 0.0 { num } else { num / (max - min) }
}

pub fn hex_to_rgb(hex: &str) -> Result<Rgb<u8>, Box<dyn Error>> {
    if !hex.starts_with("#") {
        return Err("`hex` must start with `#`".into());
    }
    if hex.len() != 7 {
        return Err("`hex` must be six characters long, excluding `#`".into());
    }
    let hexu32 = hex.replace("#", "00");
    let bytes = u32::from_str_radix(&hexu32, 16)?.to_be_bytes();
    Ok([bytes[1], bytes[2], bytes[3]].into())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::Spin;
    use std::collections::HashSet;

    #[test]
    fn test_spin_to_rgb() {
        let mut tested = HashSet::<Rgb<u8>>::default();
        for i in 0..5232 as Spin {
            let rgb = SpinPlot::spin_to_rgb(i);
            assert!(!tested.contains(&rgb));
            tested.insert(rgb);
        }
    }

    #[test]
    fn test_hex_to_rgb() {
        assert_eq!(hex_to_rgb("#ff00ff").unwrap(), [255, 0, 255].into());
    }
}
