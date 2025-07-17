use crate::constants::Spin;
use crate::environment::{Environment, LatticeEntity};
use crate::positional::boundary::Boundary;
use crate::positional::neighbourhood::Neighbourhood;
use crate::positional::pos::Pos;
use crate::spin_table::SpinTable;
use image::{Rgba, RgbaImage};
use imageproc::drawing::{draw_cross_mut, draw_line_segment_mut};
use palette::{FromColor, IntoColor, Luv, Mix, Srgb, WithAlpha};
use std::error::Error;
use std::fmt::Debug;
use std::hash::{DefaultHasher, Hash, Hasher};
use thiserror::Error;

pub trait Plot {
    fn plot(&self, image: &mut RgbaImage);
}

pub trait ContinuousPlot: Plot {
    fn min_color(&self) -> &Srgb<u8>;
    fn max_color(&self) -> &Srgb<u8>;
    fn lerp(&self, value: f32, min: f32, max: f32) -> Result<Srgb<u8>, LerpError> {
        if max < min {
            return Err(LerpError::NegativeRange);
        }
        if value < min {
            return Err(LerpError::ValueTooSmall);
        }
        if value > max {
            return Err(LerpError::ValueTooLarge);
        }

        let min_luv = Luv::from_color(self.min_color().into_linear::<f32>());
        let max_luv = Luv::from_color(self.max_color().into_linear::<f32>());
        let p = if min == max { 0.5 } else { (value - min) / (max - min) };
        let blended = min_luv.mix(max_luv, p);
        Ok(Srgb::from_linear(blended.into_color()))
    }
}

#[derive(Debug)]
pub enum LerpError {
    ValueTooSmall,
    ValueTooLarge,
    NegativeRange
}

pub struct SpinPlot<'e> {
    pub env: &'e Environment,
    pub solid_color: Srgb<u8>,
    pub medium_color: Option<Srgb<u8>>
}

impl<'e> SpinPlot<'e> {
    fn spin_to_rgb(spin: Spin) -> Srgb<u8> {
        let mut hasher = DefaultHasher::new();
        spin.hash(&mut hasher);
        let hashed = hasher.finish();
        Srgb::new(
            (hashed & 0xFF) as u8,
            (hashed >> 8 & 0xFF) as u8,
            (hashed >> 16 & 0xFF) as u8
        )
    }
}

impl Plot for SpinPlot<'_> {
    fn plot(&self, image: &mut RgbaImage) {
        for pos in self.env.space.cell_lattice.iter_positions() {
            let spin = self.env.space.cell_lattice[pos];
            let rgb = if spin > LatticeEntity::first_cell_spin() {
                Some(Self::spin_to_rgb(spin))
            } else if spin == LatticeEntity::Solid.spin() {
                Some(self.solid_color)
            } else {
                self.medium_color
            };
            if let Some(color) = rgb {
                image.put_pixel(pos.x as u32, pos.y as u32, srgb_to_rgba(color));
            }
        }
    }
}

pub struct CenterPlot<'e> {
    pub env: &'e Environment,
    pub color: Srgb<u8>
}

impl Plot for CenterPlot<'_> {
    fn plot(&self, image: &mut RgbaImage) {
        for cell in &self.env.cells {
            let center = self.env.space.lat_bound.valid_pos(Pos::new(
                cell.center.pos().x as isize,
                cell.center.pos().y as isize,
            ));
            if let Some(pos) = center {
                draw_cross_mut(image, srgb_to_rgba(self.color), pos.x as i32, pos.y as i32);
            }
        }
    }
}

pub struct ClonesPlot<'a> {
    pub env: &'a Environment,
    pub clone_pairs: &'a SpinTable<bool>,
    pub color: Srgb<u8>,
    pub all_clones: bool
}

impl Plot for ClonesPlot<'_> {
    fn plot(&self, image: &mut RgbaImage) {
        for (spin1, spin2) in self.clone_pairs.iter_pairs() {
            let message = "non-cell stored as clone";
            let cell1 = self.env.cells.get_entity(spin1).expect_cell(message);
            let cell2 = self.env.cells.get_entity(spin2).expect_cell(message);
            let center1 = cell1.center.pos;
            let center2 = cell2.center.pos;
            let mut draw = true;
            if !self.all_clones {
                let dist2 = (center1.x - center2.x).powf(2.) + (center1.y - center2.y).powf(2.);
                let mean_area2 = (cell1.area + cell2.area).pow(2);
                if dist2 > mean_area2 as f32 {
                    draw = false
                }
            }
            if draw {
                draw_line_segment_mut(
                    image,
                    (center1.x, center1.y),
                    (center2.x, center2.y),
                    srgb_to_rgba(self.color)
                )
            }
        }
    }
}

pub struct AreaPlot<'e> {
    pub env: &'e Environment,
    pub min_color: Srgb<u8>,
    pub max_color: Srgb<u8>
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
        // TODO: might be faster to iterate cells instead
        for pos in self.env.space.cell_lattice.iter_positions() {
            let entity = self.env.cells.get_entity(self.env.space.cell_lattice[pos]);
            if let LatticeEntity::SomeCell(cell) = entity {
                let color = self.lerp(
                    cell.area as f32,
                    min as f32,
                    max as f32
                );
                match color {
                    Ok(c) => image.put_pixel(
                        pos.x as u32,
                        pos.y as u32,
                        srgb_to_rgba(c)
                    ),
                    Err(e) => log::warn!("Failed to plot area for pos `{pos:?}` with error `{e:?}`")
                };
            }
        }
    }
}

impl ContinuousPlot for AreaPlot<'_> {
    fn min_color(&self) -> &Srgb<u8> {
        &self.min_color
    }

    fn max_color(&self) -> &Srgb<u8> {
        &self.max_color
    }
}

pub struct BorderPlot<'e> {
    pub env: &'e Environment,
    pub color: Srgb<u8>
}

impl Plot for BorderPlot<'_> {
    fn plot(&self, image: &mut RgbaImage) {
        for pos in self.env.space.cell_lattice.iter_positions() {
            let spin = self.env.space.cell_lattice[pos];
            if spin < LatticeEntity::first_cell_spin() {
                continue
            }
            let is_border = self.env
                .space
                .lat_bound
                .valid_positions(self.env.neighbourhood.neighbours(Pos::from(pos)))
                .any(|neigh| {
                    let neigh_spin = self.env.space.cell_lattice[Pos::from(neigh)];
                    neigh_spin != spin
                });
            if is_border {
                image.put_pixel(pos.x as u32, pos.y as u32, srgb_to_rgba(self.color));
            }
        }
    }
}

pub struct LightPlot<'e> {
    pub(crate) env: &'e Environment,
    pub(crate) min_color: Srgb<u8>,
    pub(crate) max_color: Srgb<u8>
}

impl Plot for LightPlot<'_> {
    fn plot(&self, image: &mut RgbaImage) {
        for pos in self.env.space.light_lattice.iter_positions() {
            let light = self.env.space.light_lattice[pos];
            let color = self.lerp(
                light as f32,
                0.,
                self.env.height() as f32
            );
            match color { 
                Ok(c) => image.put_pixel(pos.x as u32, pos.y as u32, srgb_to_rgba(c)),
                Err(e) => log::warn!("Failed to plot light for pos `{pos:?}` with error `{e:?}`")
            }
        }
    }
}

impl ContinuousPlot for LightPlot<'_> {
    fn min_color(&self) -> &Srgb<u8> {
        &self.min_color
    }

    fn max_color(&self) -> &Srgb<u8> {
        &self.max_color
    }
}

pub fn srgb_to_rgba(color: Srgb<u8>) -> Rgba<u8> {
    let arr: [u8; 4] = color.with_alpha(255).into();
    Rgba::from(arr)
}

pub fn hex_to_srgb(hex: &str) -> Result<Srgb<u8>, Box<dyn Error>> {
    if !hex.starts_with("#") {
        return Err(HexError::MissingHashtag.into());
    }
    if hex.len() != 7 {
        return Err(HexError::WrongLength.into());
    }
    let hexu32 = hex.replace("#", "00");
    let bytes = u32::from_str_radix(&hexu32, 16)?.to_be_bytes();
    Ok([bytes[1], bytes[2], bytes[3]].into())
}

#[derive(Error, Debug)]
pub enum HexError {
    #[error("Missing `#` in the color name")]
    MissingHashtag,
    #[error("`hex` must be six characters long, excluding `#`")]
    WrongLength
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::Spin;
    use std::collections::HashSet;

    #[test]
    fn test_spin_to_rgb() {
        let mut tested = HashSet::<[u8; 3]>::default();
        for i in 0..5232 as Spin {
            let rgb: [u8; 3] = SpinPlot::spin_to_rgb(i).into();
            assert!(!tested.contains(&rgb));
            tested.insert(rgb);
        }
    }

    #[test]
    fn test_hex_to_rgb() {
        assert_eq!(hex_to_srgb("#ff00ff").unwrap(), Srgb::new(255, 0, 255));
    }
}
