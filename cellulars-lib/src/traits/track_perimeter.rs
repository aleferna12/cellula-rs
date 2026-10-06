//! Contains logic associated with [`TrackPerimeter`].

use crate::base::environment::Environment;
use crate::positional::boundaries::ToLatticeBoundary;
use crate::positional::neighbourhood::Neighbourhood;
use crate::positional::pos::Pos;
use crate::spin::Spin;

/// Cells that keep track of their perimeter.
///
/// The perimeter of a cell is the number of edges between the positions it owns and positions owned by a
/// different [`Spin`], so it is measured in units of the environment's
/// [`Neighbourhood`] (a lone cell in a [`MooreNeighbourhood`](crate::positional::neighbourhood::MooreNeighbourhood)
/// of radius 1 has a perimeter of 8).
///
/// Cells cannot maintain this quantity on their own (they know nothing about the lattice around them), so
/// implementors are expected to have their perimeter shifted by whoever owns the lattice, with the
/// differentials computed by [`perimeter_deltas()`].
pub trait TrackPerimeter {
    /// Returns the current perimeter of the cell.
    fn perimeter(&self) -> u32;

    /// Returns the target perimeter of the cell.
    fn target_perimeter(&self) -> u32;

    /// Shifts the perimeter of the cell by `delta`.
    fn shift_perimeter(&mut self, delta: i32);
}

/// Perimeter differentials of the two spins involved in a
/// [`grant_position()`](crate::traits::habitable::Habitable::grant_position()).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Hash)]
pub struct PerimeterDeltas {
    /// Perimeter differential of the spin that gains the position.
    pub source: i32,
    /// Perimeter differential of the spin that loses the position.
    pub target: i32
}

/// Returns by how much the perimeters of the spins involved change if `to` were to be copied into `pos`.
///
/// Both differentials are 0 if `to` already owns `pos`.
///
/// Only the spins gaining and losing `pos` can have their perimeter changed by the copy: an edge between
/// `pos` and a third spin is either counted before and after the copy or neither.
///
/// <div class="warning">
/// This must be called before the cell lattice is updated.
/// </div>
pub fn perimeter_deltas<C, N: Neighbourhood, B: ToLatticeBoundary>(
    env: &Environment<C, N, B>,
    pos: Pos<usize>,
    to: Spin
) -> PerimeterDeltas {
    let from = env.cell_lattice[pos];
    if from == to {
        return PerimeterDeltas::default();
    }
    let mut deltas = PerimeterDeltas::default();
    for neigh in env.valid_neighbours(pos) {
        let neigh_spin = env.cell_lattice[neigh];
        // Edges with the new owner become internal, all others become part of its perimeter
        deltas.source += if neigh_spin == to { -1 } else { 1 };
        // The opposite happens to the previous owner
        deltas.target += if neigh_spin == from { 1 } else { -1 };
    }
    deltas
}

/// Measures the perimeter of every cell of `env` straight from its cell lattice, indexed by
/// [`CellIndex`](crate::constants::CellIndex).
///
/// Perimeters are normally maintained incrementally with [`perimeter_deltas()`], so this is meant for when
/// only the lattice is known (when restoring a simulation from a back-up, for instance) or to validate
/// that the tracked perimeters have not drifted from the lattice.
pub fn measure_perimeters<C, N: Neighbourhood, B: ToLatticeBoundary>(
    env: &Environment<C, N, B>
) -> Vec<u32> {
    let mut perimeters = vec![0; env.cells.n_cells() as usize];
    for pos in env.cell_lattice.iter_positions() {
        let Spin::Some(cell_index) = env.cell_lattice[pos] else {
            continue;
        };
        for neigh in env.valid_neighbours(pos) {
            if env.cell_lattice[neigh] != Spin::Some(cell_index) {
                perimeters[cell_index as usize] += 1;
            }
        }
    }
    perimeters
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::cell::Cell;
    use crate::constants::FloatType;
    use crate::positional::boundaries::{Boundaries, UnsafePeriodicBoundary};
    use crate::positional::neighbourhood::MooreNeighbourhood;
    use crate::positional::rect::Rect;
    use crate::traits::habitable::Habitable;

    fn make_test_env() -> Environment<Cell, MooreNeighbourhood, UnsafePeriodicBoundary<FloatType>> {
        let mut env = Environment::new_empty(
            MooreNeighbourhood::new(1),
            Boundaries::new(UnsafePeriodicBoundary::new(
                Rect::new(Pos::new(0., 0.), Pos::new(10., 10.))
            ))
        );
        env.cells.push(Cell::new_empty(4));
        env.cells.push(Cell::new_empty(4));
        env
    }

    #[test]
    fn lone_position_has_full_perimeter() {
        let mut env = make_test_env();
        let pos = Pos::new(5, 5);
        assert_eq!(perimeter_deltas(&env, pos, Spin::Some(0)).source, 8);

        env.grant_position(pos, Spin::Some(0));
        // Removing the only position of the cell brings its perimeter back to 0
        assert_eq!(perimeter_deltas(&env, pos, Spin::Medium).target, -8);
    }

    #[test]
    fn growing_along_a_side_shares_edges() {
        let mut env = make_test_env();
        env.grant_position(Pos::new(5, 5), Spin::Some(0));
        // Only 1 of the 8 neighbours of the new position belongs to the cell, so it shares one edge with it
        let deltas = perimeter_deltas(&env, Pos::new(6, 5), Spin::Some(0));
        assert_eq!(deltas.source, 7 - 1);
    }

    #[test]
    fn taking_a_position_from_another_cell() {
        let mut env = make_test_env();
        env.grant_position(Pos::new(5, 5), Spin::Some(0));
        env.grant_position(Pos::new(6, 5), Spin::Some(1));

        let deltas = perimeter_deltas(&env, Pos::new(6, 5), Spin::Some(0));
        // Cell 0 gains the position and shares one edge with it
        assert_eq!(deltas.source, 7 - 1);
        // Cell 1 is left with nothing
        assert_eq!(deltas.target, -8);
    }

    #[test]
    fn no_op_copy_changes_nothing() {
        let mut env = make_test_env();
        env.grant_position(Pos::new(5, 5), Spin::Some(0));
        assert_eq!(perimeter_deltas(&env, Pos::new(5, 5), Spin::Some(0)), PerimeterDeltas::default());
    }

    #[test]
    fn measured_perimeters_match_the_deltas() {
        let mut env = make_test_env();
        // A 2x2 square of cell 0 next to a lone position of cell 1
        for pos in [Pos::new(5, 5), Pos::new(6, 5), Pos::new(5, 6), Pos::new(6, 6)] {
            env.grant_position(pos, Spin::Some(0));
        }
        env.grant_position(Pos::new(8, 8), Spin::Some(1));

        let perimeters = measure_perimeters(&env);
        // Each of the 4 positions of the square has 8 neighbours, 3 of which are owned by the cell
        assert_eq!(perimeters[0], 4 * (8 - 3));
        assert_eq!(perimeters[1], 8);
    }
}
