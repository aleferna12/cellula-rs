//! Contains logic used to plot information to an [`RgbaImage`].

// TODO! parallelize this and all other lattice-wide operations

use crate::constants::{CellIndex, FloatType};
use crate::empty_cell::Empty;
use crate::io::write::image::lerper::Lerper;
use crate::prelude::{Cellular, Habitable, HasCenter, Spin};
use image::{Rgba, RgbaImage};
use imageproc::drawing::draw_cross_mut;
use palette::{IntoColor, Mix, Srgba};
use std::hash::{DefaultHasher, Hash, Hasher};

/// A trait to plot information about the simulation.
pub trait Plot<P> {
    /// Plots the information in `plottable` by drawing on `image`.
    fn plot(&self, plottable: &P, image: &mut RgbaImage);
}

/// Plots the spin of cells in random colors (except for the solid and medium spin colors, which can be chosen).
#[derive(Clone, PartialEq, Debug)]
pub struct SpinPlot {
    /// Color used for [`Spin::Solid`].
    pub solid_color: Srgba<FloatType>,
    /// Color used for [`Spin::Medium`].
    pub medium_color: Option<Srgba<FloatType>>
}

impl SpinPlot {
    fn cell_index_to_rgba(index: CellIndex) -> Srgba<FloatType> {
        let mut hasher = DefaultHasher::new();
        index.hash(&mut hasher);
        let hashed = hasher.finish();
        Srgba::new(
            (hashed & 0xFF) as u8,
            (hashed >> 8 & 0xFF) as u8,
            (hashed >> 16 & 0xFF) as u8,
            255
        ).into_format()
    }
}

impl<E: Habitable> Plot<E> for SpinPlot {
    fn plot(&self, env: &E, image: &mut RgbaImage) {
        for pos in env.env().cell_lattice.iter_positions() {
            let spin = env.env().cell_lattice[pos];
            let rgba = match spin {
                Spin::Some(cell_index) => Some(Self::cell_index_to_rgba(cell_index)),
                Spin::Solid => Some(self.solid_color),
                Spin::Medium => self.medium_color,
            };
            if let Some(color) = rgba {
                image.put_pixel(pos.x as u32, pos.y as u32, srgba_to_rgba(color));
            }
        }
    }
}

/// Plots the center of cells.
#[derive(Clone, PartialEq, Debug)]
pub struct CenterPlot {
    /// Color of the cell center.
    pub color: Srgba<FloatType>
}

impl<E: Habitable> Plot<E> for CenterPlot
where 
    E::Cell: HasCenter + Empty {
    fn plot(&self, env: &E, image: &mut RgbaImage) {
        for rel_cell in env.env().cells.iter() {
            if rel_cell.cell.is_empty() {
                continue;
            }
            let center = rel_cell
                .cell
                .center()
                .round();
            draw_cross_mut(image, srgba_to_rgba(self.color), center.x as i32, center.y as i32);
        }
    }
}

/// Plots the border of cells.
#[derive(Clone, PartialEq, Debug)]
pub struct BorderPlot {
    /// Color of the border.
    pub color: Srgba<FloatType>
}

impl<E: Habitable> Plot<E> for BorderPlot {
    fn plot(&self, env: &E, image: &mut RgbaImage) {
        for pos in env.env().cell_lattice.iter_positions() {
            let spin = env.env().cell_lattice[pos];
            let Spin::Some(cell_index) = spin else {
                continue;
            };

            let is_border = env
                .env()
                .valid_neighbors(pos)
                .any(|neigh| {
                    let neigh_spin = env.env().cell_lattice[neigh];
                    match neigh_spin {
                        Spin::Some(neigh_index) => cell_index < neigh_index,
                        _ => true
                    }
                });
            if is_border {
                image.put_pixel(pos.x as u32, pos.y as u32, srgba_to_rgba(self.color));
            }
        }
    }
}

/// Plots cell area.
#[derive(Clone, PartialEq, Debug, Eq, Hash)]
pub struct AreaPlot<C> {
    /// Lerper used for interpolation.
    pub lerper: Lerper<C>
}

impl<E, C> Plot<E> for AreaPlot<C>
where
    E: Habitable,
    E::Cell: Empty,
    C: Mix<Scalar = FloatType> + Clone + IntoColor<Srgba<FloatType>> {
    fn plot(&self, env: &E, image: &mut RgbaImage) {
        let mut min = u32::MAX;
        let mut max = 0;
        for rel_cell in env.env().cells.iter() {
            if rel_cell.cell.is_empty() {
                continue;
            }
            if rel_cell.cell.area() < min {
                min = rel_cell.cell.area()
            }
            if rel_cell.cell.area() > max {
                max = rel_cell.cell.area()
            }
        }

        for pos in env.env().cell_lattice.iter_positions() {
            if let Spin::Some(cell_index) = env.env().cell_lattice[pos] {
                let rel_cell = &env.env().cells[cell_index];
                let color = self.lerper.lerp(
                    (rel_cell.cell.area() as FloatType - min as FloatType) / (max as FloatType - min as FloatType),
                );
                match color {
                    Ok(c) => image.put_pixel(
                        pos.x as u32,
                        pos.y as u32,
                        srgba_to_rgba(c.into_color())
                    ),
                    Err(_e) => {
                        #[cfg(feature = "log")]
                        log::warn!("Failed to plot area for pos `{pos:?}` with error `{_e:?}`")
                    }
                };
            }
        }
    }
}

/// Converts a [`palette`]s [`Srgba`] into an [`image`]s [`Rgba`].
pub fn srgba_to_rgba(srgba: Srgba<FloatType>) -> Rgba<u8> {
    let srgba_u8 = srgba.into_format();
    Rgba([srgba_u8.red, srgba_u8.green, srgba_u8.blue, srgba_u8.alpha])
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::CellIndex;
    use std::collections::HashSet;

    #[test]
    fn test_index_to_rgb() {
        let mut tested = HashSet::<[u8; 4]>::default();
        for i in 0..5232 as CellIndex {
            let rgb: [u8; 4] = SpinPlot::cell_index_to_rgba(i).into_format().into();
            assert!(!tested.contains(&rgb));
            tested.insert(rgb);
        }
    }
}