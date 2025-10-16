use crate::chem_environment::ChemEnvironment;
use crate::io::parameters::{PlotParameters, PlotType};
use crate::io::plot::HexError::ParseU8Error;
use cellulars_lib::basic_cell::Cellular;
use cellulars_lib::constants::CellIndex;
use cellulars_lib::positional::boundaries::Boundary;
use cellulars_lib::positional::neighbourhood::Neighbourhood;
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::spin::Spin;
use image::{Rgba, RgbaImage};
use imageproc::drawing::{draw_cross_mut, draw_line_segment_mut};
use palette::{FromColor, IntoColor, Luv, Mix, Srgb, WithAlpha};
use std::fmt::Debug;
use std::hash::{DefaultHasher, Hash, Hasher};
use thiserror::Error;

pub trait Plot {
    fn plot(&self, env: &ChemEnvironment, image: &mut RgbaImage);
}

pub trait ContinuousPlot: Plot {
    fn min_color(&self) -> Luv;
    fn max_color(&self) -> Luv;
    fn lerp(&self, value: f32, min: f32, max: f32) -> Result<Luv, LerpError> {
        if max < min {
            return Err(LerpError::NegativeRange);
        }
        if value < min {
            return Err(LerpError::ValueTooSmall);
        }
        if value > max {
            return Err(LerpError::ValueTooLarge);
        }

        let p = if min == max { 0.5 } else { (value - min) / (max - min) };
        let blended = self.min_color().mix(self.max_color(), p);
        Ok(blended)
    }
}

#[derive(Debug)]
pub enum LerpError {
    ValueTooSmall,
    ValueTooLarge,
    NegativeRange
}

pub struct SpinPlot {
    pub solid_color: Srgb<u8>,
    pub medium_color: Option<Srgb<u8>>
}

impl SpinPlot {
    fn cell_index_to_rgb(index: CellIndex) -> Srgb<u8> {
        let mut hasher = DefaultHasher::new();
        index.hash(&mut hasher);
        let hashed = hasher.finish();
        Srgb::new(
            (hashed & 0xFF) as u8,
            (hashed >> 8 & 0xFF) as u8,
            (hashed >> 16 & 0xFF) as u8
        )
    }
}

impl Plot for SpinPlot {
    fn plot(&self, env: &ChemEnvironment, image: &mut RgbaImage) {
        for pos in env.cell_lattice.iter_positions() {
            let spin = env.cell_lattice[pos];
            let rgb = match spin {
                Spin::Some(cell_index) => Some(Self::cell_index_to_rgb(cell_index)),
                Spin::Solid => Some(self.solid_color),
                Spin::Medium => self.medium_color
            };
            if let Some(color) = rgb {
                image.put_pixel(pos.x as u32, pos.y as u32, srgb_to_rgba(color));
            }
        }
    }
}

pub struct CenterPlot {
    pub color: Srgb<u8>
}

impl Plot for CenterPlot {
    fn plot(&self, env: &ChemEnvironment, image: &mut RgbaImage) {
        for cell in env.cells.iter() {
            if !cell.is_valid() {
                continue;
            }
            let center = env.bounds.lattice_boundary.valid_pos(Pos::new(
                cell.center().x as isize,
                cell.center().y as isize,
            ));
            if let Some(pos) = center {
                draw_cross_mut(image, srgb_to_rgba(self.color), pos.x as i32, pos.y as i32);
            }
        }
    }
}

pub struct ChemCenterPlot {
    pub color: Srgb<u8>
}

impl Plot for ChemCenterPlot {
    fn plot(&self, env: &ChemEnvironment, image: &mut RgbaImage) {
        for cell in env.cells.iter() {
            if !cell.is_valid() {
                continue;
            }
            let center = env.bounds.lattice_boundary.valid_pos(Pos::new(
                cell.chem_center().x as isize,
                cell.chem_center().y as isize,
            ));
            if let Some(pos) = center {
                draw_cross_mut(image, srgb_to_rgba(self.color), pos.x as i32, pos.y as i32);
            }
        }
    }
}

pub struct ClonesPlot {
    pub color: Srgb<u8>,
    pub all_clones: bool
}

impl Plot for ClonesPlot {
    fn plot(&self, env: &ChemEnvironment, image: &mut RgbaImage) {
        let clones = &env.clones_table;
        for index_pair in clones.iter_index_pairs(None, None) {
            if !clones[index_pair] {
                continue;
            }
            let cell1 = env.cells.get_cell(index_pair.0 as CellIndex);
            let cell2 = env.cells.get_cell(index_pair.1 as CellIndex);
            if !cell1.is_valid() || !cell2.is_valid() {
                continue;
            }
            let center1 = cell1.center();
            let center2 = cell2.center();
            if !self.all_clones {
                let dist2 = (center1.x - center2.x).powf(2.) + (center1.y - center2.y).powf(2.);
                let diag = env.width() * env.width() + env.height() * env.height();
                if dist2 > diag as f32 / 4. {
                    continue;
                }
            }
            draw_line_segment_mut(
                image,
                (center1.x, center1.y),
                (center2.x, center2.y),
                srgb_to_rgba(self.color)
            )
        }
    }
}

pub struct BorderPlot {
    pub color: Srgb<u8>
}

impl Plot for BorderPlot {
    fn plot(&self, env: &ChemEnvironment, image: &mut RgbaImage) {
        for pos in env.cell_lattice.iter_positions() {
            let spin = env.cell_lattice[pos];
            if !matches!(spin, Spin::Some(_)) {
                continue
            }
            let is_border = env
                .bounds
                .lattice_boundary
                .valid_positions(env.neighbourhood.neighbours(pos.to_isize()))
                .any(|neigh| {
                    let neigh_spin = env.cell_lattice[neigh.to_usize()];
                    neigh_spin != spin
                });
            if is_border {
                image.put_pixel(pos.x as u32, pos.y as u32, srgb_to_rgba(self.color));
            }
        }
    }
}

pub struct CellTypePlot {
    pub mig_color: Srgb<u8>,
    pub div_color: Srgb<u8>
}

impl Plot for CellTypePlot {
    fn plot(&self, env: &ChemEnvironment, image: &mut RgbaImage) {
        for pos in env.cell_lattice.iter_positions() {
            let spin = env.cell_lattice[pos];
            if let Spin::Some(cell_index) = spin {
                let cell = env.cells.get_cell(cell_index);
                let color = if cell.is_dividing() { self.div_color } else { self.mig_color };
                image.put_pixel(
                    pos.x as u32,
                    pos.y as u32,
                    srgb_to_rgba(color)
                )
            }
        }
    }
}

pub struct AreaPlot {
    pub min_color: Luv,
    pub max_color: Luv
}

impl Plot for AreaPlot {
    fn plot(&self, env: &ChemEnvironment, image: &mut RgbaImage) {
        let mut min = u32::MAX;
        let mut max = 0;
        for cell in env.cells.iter() {
            if !cell.is_valid() {
                continue;
            }
            if cell.area() < min {
                min = cell.area()
            }
            if cell.area() > max {
                max = cell.area()
            }
        }

        for pos in env.cell_lattice.iter_positions() {
            if let Spin::Some(cell_index) = env.cell_lattice[pos] {
                let cell = env.cells.get_cell(cell_index);
                let color = self.lerp(
                    cell.area() as f32,
                    min as f32,
                    max as f32
                );
                match color {
                    Ok(c) => image.put_pixel(
                        pos.x as u32,
                        pos.y as u32,
                        srgb_to_rgba(Srgb::from_linear(c.into_color()))
                    ),
                    Err(e) => log::warn!("Failed to plot area for pos `{pos:?}` with error `{e:?}`")
                };
            }
        }
    }
}

impl ContinuousPlot for AreaPlot {
    fn min_color(&self) -> Luv {
        self.min_color
    }
    fn max_color(&self) -> Luv {
        self.max_color
    }
}

pub struct FoodPlot {
    pub min_color: Luv,
    pub max_color: Luv
}

impl Plot for FoodPlot {
    fn plot(&self, env: &ChemEnvironment, image: &mut RgbaImage) {
        let mut min = u32::MAX;
        let mut max = 0;
        for cell in env.cells.iter() {
            if !cell.is_valid() {
                continue;
            }
            if cell.food < min {
                min = cell.food
            }
            if cell.food > max {
                max = cell.food
            }
        }

        for pos in env.cell_lattice.iter_positions() {
            if let Spin::Some(cell_index) = env.cell_lattice[pos] {
                let cell = env.cells.get_cell(cell_index);
                let color = self.lerp(
                    cell.food as f32,
                    min as f32,
                    max as f32
                );
                match color {
                    Ok(c) => image.put_pixel(
                        pos.x as u32,
                        pos.y as u32,
                        srgb_to_rgba(Srgb::from_linear(c.into_color()))
                    ),
                    Err(e) => log::warn!("Failed to plot food for pos `{pos:?}` with error `{e:?}`")
                };
            }
        }
    }
}

impl ContinuousPlot for FoodPlot {
    fn min_color(&self) -> Luv {
        self.min_color
    }
    fn max_color(&self) -> Luv {
        self.max_color
    }
}

pub struct ChemPlot {
    pub min_color: Luv,
    pub max_color: Luv
}

impl Plot for ChemPlot {
    fn plot(&self, env: &ChemEnvironment, image: &mut RgbaImage) {
        let lat = &env.chem_lattice;
        for pos in lat.iter_positions() {
            let chem = lat[pos];
            let color = self.lerp(
                chem as f32,
                0.,
                lat.height() as f32
            );
            match color { 
                Ok(c) => image.put_pixel(
                    pos.x as u32, 
                    pos.y as u32,
                    srgb_to_rgba(Srgb::from_linear(c.into_color()))
                ),
                Err(e) => log::warn!("Failed to plot chem for pos `{pos:?}` with error `{e:?}`")
            }
        }
    }
}

impl ContinuousPlot for ChemPlot {
    fn min_color(&self) -> Luv {
        self.min_color
    }

    fn max_color(&self) -> Luv {
        self.max_color
    }
}

pub fn srgb_to_rgba(color: Srgb<u8>) -> Rgba<u8> {
    let arr: [u8; 4] = color.with_alpha(255).into();
    Rgba::from(arr)
}

impl TryFrom<PlotParameters> for Vec<Box<dyn Plot>> {
    type Error = HexError;

    fn try_from(params: PlotParameters) -> Result<Self, HexError> {
        let mut plots = Vec::with_capacity(params.order.len());
        for plot_type in params.order {
            let plot: Box<dyn Plot> = match plot_type {
                PlotType::Spin => Box::new(SpinPlot {
                    solid_color: hex_to_srgb(&params.solid_color)?,
                    medium_color: match &params.medium_color {
                        None => None,
                        Some(c) => Some(hex_to_srgb(c)?)
                    }
                }),
                PlotType::Center => Box::new(CenterPlot {
                    color: hex_to_srgb(&params.center_color)?
                }),
                PlotType::ChemCenter => Box::new(ChemCenterPlot {
                    color: hex_to_srgb(&params.chem_center_color)?
                }),
                PlotType::Clones => Box::new(ClonesPlot {
                    color: hex_to_srgb(&params.clones_color)?,
                    all_clones: params.all_clones
                }),
                PlotType::Area => Box::new(AreaPlot{
                    min_color: srgb_to_luv(hex_to_srgb(&params.area_min_color)?),
                    max_color: srgb_to_luv(hex_to_srgb(&params.area_max_color)?),
                }),
                PlotType::Food => Box::new(FoodPlot{
                    min_color: srgb_to_luv(hex_to_srgb(&params.food_min_color)?),
                    max_color: srgb_to_luv(hex_to_srgb(&params.food_max_color)?),
                }),
                PlotType::Border => Box::new(BorderPlot {
                    color: hex_to_srgb(&params.border_color)?
                }),
                PlotType::Chem => Box::new(ChemPlot {
                    min_color: srgb_to_luv(hex_to_srgb(&params.chem_min_color)?),
                    max_color: srgb_to_luv(hex_to_srgb(&params.chem_max_color)?)
                }),
                PlotType::CellType => Box::new(CellTypePlot {
                    mig_color: hex_to_srgb(&params.migrating_color)?,
                    div_color: hex_to_srgb(&params.dividing_color)?,
                })
            };
            plots.push(plot);
        }
        Ok(plots)
    }
}

pub fn srgb_to_luv(srgb: Srgb<u8>) -> Luv {
    Luv::from_color(srgb.into_linear::<f32>())
}

pub fn hex_to_srgb(hex: &str) -> Result<Srgb<u8>, HexError> {
    if !hex.starts_with("#") {
        return Err(HexError::MissingHashtag);
    }
    if hex.len() != 7 {
        return Err(HexError::WrongLength);
    }
    let hexu32 = hex.replace("#", "00");
    let bytes = u32::from_str_radix(&hexu32, 16).map_err(ParseU8Error)?.to_be_bytes();
    Ok([bytes[1], bytes[2], bytes[3]].into())
}

#[derive(Error, Debug)]
pub enum HexError {
    #[error("missing `#` in the color name")]
    MissingHashtag,
    #[error("`hex` must be six characters long, excluding `#`")]
    WrongLength,
    #[error("failed to parse the string as a hex `u8`: {0}")]
    ParseU8Error(#[from] std::num::ParseIntError),
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    #[test]
    fn test_index_to_rgb() {
        let mut tested = HashSet::<[u8; 3]>::default();
        for i in 0..5232 as CellIndex {
            let rgb: [u8; 3] = SpinPlot::cell_index_to_rgb(i).into();
            assert!(!tested.contains(&rgb));
            tested.insert(rgb);
        }
    }

    #[test]
    fn test_hex_to_rgb() {
        assert_eq!(hex_to_srgb("#ff00ff").unwrap(), Srgb::new(255, 0, 255));
    }
}
