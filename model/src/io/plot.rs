//! Contains logic for plotting data about the simulation.

use crate::io::parameters::{PlotParameters, PlotType};
use crate::io::plot::HexError::ParseU8Error;
use crate::my_cell::CellType;
use crate::my_environment::MyEnvironment;
use cellulars_lib::constants::FloatType;
use cellulars_lib::io::lerper::Lerper;
use cellulars_lib::io::plot::{srgba_to_rgba, AreaPlot, BorderPlot, CenterPlot, Plot, SpinPlot};
use cellulars_lib::spin::Spin;
use cellulars_lib::traits::cellular::Cellular;
use image::RgbaImage;
use imageproc::drawing::draw_cross_mut;
use palette::{Clamp, FromColor, Oklab, Srgba};
use std::fmt::Debug;
use thiserror::Error;

/// Plots the perceived chemical center of cells.
pub struct ChemCenterPlot {
    /// Color of the cell chemical center.
    pub color: Srgba<FloatType>,
}

impl Plot<MyEnvironment> for ChemCenterPlot {
    fn plot(&self, env: &MyEnvironment, image: &mut RgbaImage) {
        for rel_cell in env.env.cells.iter() {
            if !rel_cell.cell.is_empty() {
                continue;
            }
            let center = rel_cell
                .cell
                .chem_center()
                .round();
            draw_cross_mut(image, srgba_to_rgba(self.color), center.x as i32, center.y as i32);
        }
    }
}

/// Plots cells according to their cell type.
pub struct CellTypePlot {
    /// Color for the migrating cells.
    pub mig_color: Srgba<FloatType>,
    /// Color for the dividing cells.
    pub div_color: Srgba<FloatType>
}

impl Plot<MyEnvironment> for CellTypePlot {
    fn plot(&self, env: &MyEnvironment, image: &mut RgbaImage) {
        for pos in env.env.cell_lattice.iter_positions() {
            let spin = env.env.cell_lattice[pos];
            if let Spin::Some(cell_index) = spin {
                let rel_cell = &env.env.cells[cell_index];
                let color = match rel_cell.cell.cell_type {
                    CellType::Migrating => self.mig_color,
                    CellType::Dividing => self.div_color
                };
                image.put_pixel(
                    pos.x as u32,
                    pos.y as u32,
                    srgba_to_rgba(color)
                )
            }
        }
    }
}

/// Plots the chemical lattice.
pub struct ChemPlot {
    pub lerper: Lerper<Oklab<FloatType>>,
}

impl Plot<MyEnvironment> for ChemPlot {
    fn plot(&self, env: &MyEnvironment, image: &mut RgbaImage) {
        let lat = &env.chem_lattice;
        for pos in lat.iter_positions() {
            let chem = lat[pos];
            let color = self.lerper.lerp(
                chem as FloatType,
                0.,
                lat.height() as FloatType
            );
            match color {
                Ok(c) => {
                    image.put_pixel(
                        pos.x as u32,
                        pos.y as u32,
                        srgba_to_rgba(oklab_to_srgba(c)),
                    )
                },
                Err(e) => log::warn!("Failed to plot chem for pos `{pos:?}` with error `{e:?}`")
            }
        }
    }
}

impl TryFrom<PlotParameters> for Box<[Box<dyn Plot<MyEnvironment>>]> {
    type Error = HexError;

    fn try_from(params: PlotParameters) -> Result<Self, HexError> {
        let mut plots = Vec::with_capacity(params.order.len());
        for plot_type in params.order {
            let plot: Box<dyn Plot<MyEnvironment>> = match plot_type {
                PlotType::Spin => Box::new(SpinPlot {
                    solid_color: hex_to_srgba(&params.solid_color)?,
                    medium_color: match &params.medium_color {
                        None => None,
                        Some(c) => Some(hex_to_srgba(c)?)
                    }
                }),
                PlotType::Center => Box::new(CenterPlot {
                    color: hex_to_srgba(&params.center_color)?
                }),
                PlotType::ChemCenter => Box::new(ChemCenterPlot {
                    color: hex_to_srgba(&params.chem_center_color)?
                }),
                PlotType::Area => Box::new(AreaPlot::<Oklab<FloatType>> {
                    lerper: Lerper {
                        min_color: srgba_to_oklab(hex_to_srgba(&params.area_min_color)?),
                        max_color: srgba_to_oklab(hex_to_srgba(&params.area_max_color)?),
                    }
                }),
                PlotType::Border => Box::new(BorderPlot {
                    color: hex_to_srgba(&params.border_color)?
                }),
                PlotType::Chem => Box::new(ChemPlot {
                    lerper: Lerper {
                        min_color: srgba_to_oklab(hex_to_srgba(&params.chem_min_color)?),
                        max_color: srgba_to_oklab(hex_to_srgba(&params.chem_max_color)?)
                    }
                }),
                PlotType::CellType => Box::new(CellTypePlot {
                    mig_color: hex_to_srgba(&params.migrating_color)?,
                    div_color: hex_to_srgba(&params.dividing_color)?,
                })
            };
            plots.push(plot);
        }
        Ok(plots.into_boxed_slice())
    }
}

fn oklab_to_srgba(color: Oklab<FloatType>) -> Srgba<FloatType> {
    Srgba::from_color(color).clamp()
}

fn srgba_to_oklab(color: Srgba<FloatType>) -> Oklab<FloatType> {
    Oklab::from_color(color).clamp()
}

/// Parses a hex string as an [`Srgba<FloatType>`].
fn hex_to_srgba(hex: &str) -> Result<Srgba<FloatType>, HexError> {
    if !hex.starts_with("#") {
        return Err(HexError::MissingHashtag);
    }
    if hex.len() != 7 {
        return Err(HexError::WrongLength);
    }
    let hexu32 = hex.replace("#", "00");
    let bytes = u32::from_str_radix(&hexu32, 16).map_err(ParseU8Error)?.to_be_bytes();
    Ok(Srgba::new(bytes[1], bytes[2], bytes[3], 255).into_format())
}

/// Error thrown when a string could not be parsed into a [`Srgba<FloatType>`]
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

    #[test]
    fn test_hex_to_rgb() {
        assert_eq!(
            hex_to_srgba("#ff00ff").unwrap().into_format(),
            Srgba::new(255u8, 0, 255, 255)
        );
    }
}
