//! Contains logic associated with [`PerimeterConstraint`].

use crate::constants::FloatType;
use crate::positional::pos::Pos;
use crate::spin::Spin;
use crate::traits::habitable::Habitable;
use crate::traits::track_perimeter::{perimeter_deltas, TrackPerimeter};

/// Penalises deviations of a cell's perimeter from its target perimeter, which keeps cells compact and
/// discourages the long protrusions that migration biases tend to create.
///
/// Cells of the environment must [`TrackPerimeter`], since the perimeters used here are the ones the cells
/// report (they are not measured from the lattice).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PerimeterConstraint {
    /// Scaling constant associated with the penalty given to perimeter deviations.
    pub lambda: FloatType
}

impl PerimeterConstraint {
    /// Returns the energy differential resulting from shifting the perimeter of a cell by `delta_perimeter`.
    pub fn perimeter_energy_diff(
        &self,
        delta_perimeter: i32,
        perimeter: u32,
        target_perimeter: u32
    ) -> FloatType {
        let delta = delta_perimeter as FloatType;
        // Starting from diff = lambda * (p + delta - t) ** 2 - lambda * (p - t) ** 2
        self.lambda * delta * (2. * (perimeter as FloatType - target_perimeter as FloatType) + delta)
    }

    /// Returns the total energy differential of the perimeter constraint if the spin at `pos_source` were
    /// to be copied into `pos_target`.
    ///
    /// <div class="warning">
    /// This must be called before the cell lattice is updated.
    /// </div>
    pub fn energy_diff<H>(
        &self,
        pos_source: Pos<usize>,
        pos_target: Pos<usize>,
        env: &H
    ) -> FloatType
    where
        H: Habitable,
        H::Cell: TrackPerimeter {
        let env = env.env();
        let spin_source = env.cell_lattice[pos_source];
        // Recomputing these when the copy is granted is cheaper than keeping them around,
        // since most copy attempts are rejected
        let deltas = perimeter_deltas(env, pos_target, spin_source);
        let mut energy = 0.;
        if let Spin::Some(cell_index) = spin_source {
            let cell = &env.cells[cell_index].cell;
            energy += self.perimeter_energy_diff(deltas.source, cell.perimeter(), cell.target_perimeter());
        }
        if let Spin::Some(cell_index) = env.cell_lattice[pos_target] {
            let cell = &env.cells[cell_index].cell;
            energy += self.perimeter_energy_diff(deltas.target, cell.perimeter(), cell.target_perimeter());
        }
        energy
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn growing_towards_the_target_perimeter_is_favourable() {
        let constraint = PerimeterConstraint { lambda: 1. };
        // A cell below its target perimeter gains energy by shrinking and loses energy by growing
        assert!(constraint.perimeter_energy_diff(1, 10, 20) < 0.);
        assert!(constraint.perimeter_energy_diff(-1, 10, 20) > 0.);
        // And the other way around once it is above its target
        assert!(constraint.perimeter_energy_diff(1, 30, 20) > 0.);
        assert!(constraint.perimeter_energy_diff(-1, 30, 20) < 0.);
    }

    #[test]
    fn energy_is_quadratic_on_the_deviation() {
        let constraint = PerimeterConstraint { lambda: 0.5 };
        // (12 - 10) ** 2 - (10 - 10) ** 2, scaled by lambda
        assert_eq!(constraint.perimeter_energy_diff(2, 10, 10), 0.5 * 4.);
    }

    #[test]
    fn no_shift_costs_nothing() {
        let constraint = PerimeterConstraint { lambda: 3. };
        assert_eq!(constraint.perimeter_energy_diff(0, 17, 20), 0.);
    }
}
