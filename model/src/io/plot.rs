//! Contains logic for plotting data about the simulation.

use crate::cell::CellType;
use crate::environment::Environment;
use crate::io::parameters::{PlotParameters, PlotType};
use crate::io::plot::HexError::ParseU8Error;
use cellulars_lib::constants::{CellIndex, FloatType};
use cellulars_lib::spin::Spin;
use cellulars_lib::traits::cellular::Cellular;
use image::{Rgba, RgbaImage};
use imageproc::drawing::draw_cross_mut;
use palette::{FromColor, IntoColor, Luv, Mix, Srgb, WithAlpha};
use std::fmt::Debug;
use std::hash::{DefaultHasher, Hash, Hasher};
use thiserror::Error;

/// A trait to plot information about the environment.
pub trait Plot {
    /// Plots the information in `env` by drawing on `image`.
    fn plot(&self, env: &Environment, image: &mut RgbaImage);
}

/// [`Plot`]s that can display continuous variables.
pub trait ContinuousPlot: Plot {
    /// Color for when `value == min`.
    fn min_color(&self) -> Luv;
    /// Color for when `value == max`.
    fn max_color(&self) -> Luv;
    /// Linearly interpolates `value` between `min` and `max`.
    fn lerp(&self, value: FloatType, min: FloatType, max: FloatType) -> Result<Luv, LerpError> {
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
        let blended = self.min_color().mix(self.max_color(), p as f32);
        Ok(blended)
    }
}

/// Error thrown when linear interpolation fails.
#[derive(Debug)]
pub enum LerpError {
    /// Value falls outside the range because it's too small.
    ValueTooSmall,
    /// Value falls outside the range because it's too large.
    ValueTooLarge,
    /// Minimum value passed is larger than maximum.
    NegativeRange
}

/// Plots the spin of cells in random colors (except for the solid and medium spin colors, which can be chosen).
pub struct SpinPlot {
    /// Color used for [`Spin::Solid`].
    pub solid_color: Srgb<u8>,
    /// Color used for [`Spin::Medium`].
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
    fn plot(&self, env: &Environment, image: &mut RgbaImage) {
        for pos in env.base_env.cell_lattice.iter_positions() {
            let spin = env.base_env.cell_lattice[pos];
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

/// Plots the center of cells.
pub struct CenterPlot {
    /// Color of the cell center.
    pub color: Srgb<u8>
}

impl Plot for CenterPlot {
    fn plot(&self, env: &Environment, image: &mut RgbaImage) {
        for rel_cell in env.base_env.cells.iter() {
            if !rel_cell.cell.is_empty() {
                continue;
            }
            let center = rel_cell
                .cell
                .center()
                .round();
            draw_cross_mut(image, srgb_to_rgba(self.color), center.x as i32, center.y as i32);
        }
    }
}

/// Plots the perceived chemical center of cells.
pub struct ChemCenterPlot {
    /// Color of the cell chemical center.
    pub color: Srgb<u8>
}

impl Plot for ChemCenterPlot {
    fn plot(&self, env: &Environment, image: &mut RgbaImage) {
        for rel_cell in env.base_env.cells.iter() {
            if !rel_cell.cell.is_empty() {
                continue;
            }
            let center = rel_cell
                .cell
                .chem_center()
                .round();
            draw_cross_mut(image, srgb_to_rgba(self.color), center.x as i32, center.y as i32);
        }
    }
}

/// Plots the border of cells.
pub struct BorderPlot {
    /// Color of the border.
    pub color: Srgb<u8>
}

impl Plot for BorderPlot {
    fn plot(&self, env: &Environment, image: &mut RgbaImage) {
        for pos in env.base_env.cell_lattice.iter_positions() {
            let spin = env.base_env.cell_lattice[pos];
            let Spin::Some(cell_index) = spin else {
                continue;
            };

            let is_border = env
                .base_env
                .valid_neighbours(pos)
                .any(|neigh| {
                    let neigh_spin = env.base_env.cell_lattice[neigh];
                    match neigh_spin {
                        Spin::Some(neigh_index) => cell_index < neigh_index,
                        _ => true
                    }
                });
            if is_border {
                image.put_pixel(pos.x as u32, pos.y as u32, srgb_to_rgba(self.color));
            }
        }
    }
}

/// Plots cells according to their cell type.
pub struct CellTypePlot {
    /// Color for the migrating cells.
    pub mig_color: Srgb<u8>,
    /// Color for the dividing cells.
    pub div_color: Srgb<u8>
}

impl Plot for CellTypePlot {
    fn plot(&self, env: &Environment, image: &mut RgbaImage) {
        for pos in env.base_env.cell_lattice.iter_positions() {
            let spin = env.base_env.cell_lattice[pos];
            if let Spin::Some(cell_index) = spin {
                let rel_cell = &env.base_env.cells[cell_index];
                let color = match rel_cell.cell.cell_type {
                    CellType::Migrating => self.mig_color,
                    CellType::Dividing => self.div_color
                };
                image.put_pixel(
                    pos.x as u32,
                    pos.y as u32,
                    srgb_to_rgba(color)
                )
            }
        }
    }
}

/// Plots cell area.
pub struct AreaPlot {
    /// Color used to display the smallest value of the plot.
    pub min_color: Luv,
    /// Color used to display the largest value of the plot.
    pub max_color: Luv
}

impl Plot for AreaPlot {
    fn plot(&self, env: &Environment, image: &mut RgbaImage) {
        let mut min = u32::MAX;
        let mut max = 0;
        for rel_cell in env.base_env.cells.iter() {
            if !rel_cell.cell.is_empty() {
                continue;
            }
            if rel_cell.cell.area() < min {
                min = rel_cell.cell.area()
            }
            if rel_cell.cell.area() > max {
                max = rel_cell.cell.area()
            }
        }

        for pos in env.base_env.cell_lattice.iter_positions() {
            if let Spin::Some(cell_index) = env.base_env.cell_lattice[pos] {
                let rel_cell = &env.base_env.cells[cell_index];
                let color = self.lerp(
                    rel_cell.cell.area() as FloatType,
                    min as FloatType,
                    max as FloatType
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

/// Plots the chemical lattice.
pub struct ChemPlot {
    /// Color used to display the smallest value of the plot.
    pub min_color: Luv,
    /// Color used to display the largest value of the plot.
    pub max_color: Luv
}

impl Plot for ChemPlot {
    fn plot(&self, env: &Environment, image: &mut RgbaImage) {
        let lat = &env.chem_lattice;
        for pos in lat.iter_positions() {
            let chem = lat[pos];
            let color = self.lerp(
                chem as FloatType,
                0.,
                lat.height() as FloatType
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

/// Adds an alpha = 255 component to the `color`.
pub fn srgb_to_rgba(color: Srgb<u8>) -> Rgba<u8> {
    let arr: [u8; 4] = color.with_alpha(255).into();
    Rgba::from(arr)
}

impl TryFrom<PlotParameters> for Box<[Box<dyn Plot>]> {
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
                PlotType::Area => Box::new(AreaPlot{
                    min_color: srgb_to_luv(hex_to_srgb(&params.area_min_color)?),
                    max_color: srgb_to_luv(hex_to_srgb(&params.area_max_color)?),
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
        Ok(plots.into_boxed_slice())
    }
}

/// Converts [`Srgb<u8>`] to [`Luv`].
pub fn srgb_to_luv(srgb: Srgb<u8>) -> Luv {
    Luv::from_color(srgb.into_linear::<f32>())
}

/// Parses a hex string as an [`Srgb<u8>`].
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

/// Error thrown when a string could not be parsed into a [`Srgb<u8>`]
#[derive(Error, Debug)]
pub enum HexError {
    /// Hex string is missing "#" in the beginning.
    #[error("missing \"#\" in the color name")]
    MissingHashtag,
    /// Hex string has the wrong length.
    #[error("`hex` must be six characters long, excluding `#`")]
    WrongLength,
    /// Failed to parse the hex string as a number.
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
