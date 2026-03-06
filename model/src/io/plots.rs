//! Contains logic for plotting data about the simulation.

use crate::io::parameters::{PlotParameters, PlotType};
use crate::my_cell::CellType;
use crate::my_environment::MyEnvironment;
use anyhow::{anyhow, bail};
use cellulars::constants::FloatType;
use cellulars::empty_cell::Empty;
use cellulars::io::write::image::lerper::Lerper;
use cellulars::io::write::image::plot::{srgba_to_rgba, AreaPlot, BorderPlot, CenterPlot, Plot, SpinPlot};
use cellulars::prelude::{Boundary, FixedBoundary, HasCenter, Pos};
use cellulars::spin::Spin;
use image::RgbaImage;
use imageproc::drawing::draw_cross_mut;
use palette::Clamp;
use palette::{FromColor, Oklab, Srgba};

/// Plots the perceived chemical center of cells.
pub struct ChemCenterPlot {
    /// Color of the cell chemical center.
    pub color: Srgba<FloatType>,
}

impl Plot<MyEnvironment> for ChemCenterPlot {
    fn plot(&self, env: &MyEnvironment, image: &mut RgbaImage) {
        let color = srgba_to_rgba(self.color);
        for rel_cell in env.env.cells.iter() {
            if rel_cell.cell.is_empty() {
                continue;
            }
            let center = rel_cell
                .cell
                .chem_center()
                .round();
            draw_cross_mut(image, color, center.x as i32, center.y as i32);
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
    /// Interpolator struct.
    pub lerper: Lerper<Oklab<FloatType>>,
}

impl Plot<MyEnvironment> for ChemPlot {
    fn plot(&self, env: &MyEnvironment, image: &mut RgbaImage) {
        let lat = &env.chem_lattice;
        for pos in lat.iter_positions() {
            let chem = lat[pos];
            let color = self.lerper.lerp(
                chem as FloatType / lat.height() as FloatType,
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

pub struct DivisionAxisPlot {
    color: Srgba<FloatType>,
    length: i32
}

impl Plot<MyEnvironment> for DivisionAxisPlot {
    fn plot(&self, env: &MyEnvironment, image: &mut RgbaImage) {
        let color = srgba_to_rgba(self.color);
        for rel_cell in env.env.cells.iter() {
            let div_axis = env.find_division_axis(rel_cell);
            let center = rel_cell.cell.center();
            let fixed_bound = FixedBoundary::new(env.env.bounds.boundary.rect().clone());

            for t in -self.length / 2..self.length / 2 {
                let x = center.x + div_axis.0 * t as FloatType;
                let y = center.y + div_axis.1 * t as FloatType;

                let Some(valid_pos) = fixed_bound.valid_pos(Pos::new(x, y)) else {
                    continue;
                };
                let lat_pos = valid_pos.cast_as();
                image.put_pixel(lat_pos.x, lat_pos.y, color);
            }
        }
    }
}

impl TryFrom<PlotParameters> for Box<[Box<dyn Plot<MyEnvironment>>]> {
    type Error = anyhow::Error;

    fn try_from(params: PlotParameters) -> anyhow::Result<Self> {
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
                }),
                PlotType::DivisionAxis => Box::new(DivisionAxisPlot {
                    color: hex_to_srgba(&params.division_axis_color)?,
                    length: params.division_axis_length
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
fn hex_to_srgba(hex: &str) -> anyhow::Result<Srgba<FloatType>> {
    if hex.len() != 7 {
        bail!("Hex string must have length 7");
    }
    let hexu32 = hex
        .strip_prefix("#")
        .ok_or(anyhow!("Missing starting `#` character in hex string"))?;
    let bytes = u32::from_str_radix(hexu32, 16)
        .map_err(anyhow::Error::new)?
        .to_be_bytes();
    Ok(Srgba::new(bytes[1], bytes[2], bytes[3], 255).into_format())
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
        assert_eq!(
            hex_to_srgba("#79933b").unwrap().into_format(),
            Srgba::new(121u8, 147, 59, 255)
        );
    }
}