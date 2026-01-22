//! Contains logic associated with [`Habitable`].

use crate::base::environment::{Environment, EdgesUpdate};
use crate::cell_container::RelCell;
use crate::positional::boundaries::ToLatticeBoundary;
use crate::positional::neighbourhood::Neighbourhood;
use crate::positional::pos::Pos;
use crate::spin::Spin;
use crate::traits::cellular::{Cellular, EmptyCell};

/// This trait asserts that a type is habitable,
/// which is to say that it can contain active cells.
///
/// Overriding methods of this trait (especially [Habitable::grant_position()])
/// allows for custom logic of how to update the simulation.
pub trait Habitable {
    /// Cell type of the environment associated with this trait.
    type Cell: Cellular;

    /// Returns a reference to the environment where cells live.
    fn env(&self) -> &Environment<Self::Cell, impl Neighbourhood, impl ToLatticeBoundary>;

    /// Returns a mutable reference to the environment  where cells live.
    fn env_mut(&mut self) -> &mut Environment<Self::Cell, impl Neighbourhood, impl ToLatticeBoundary>;

    /// Grants position `pos` to the entity represented by spin `to`.
    ///
    /// This method should be used whenever a position on the environment changes ownership.
    ///
    /// Assumes that `pos` is a valid position in the environment's lattice.
    fn grant_position(&mut self, pos: Pos<usize>, to: Spin) -> EdgesUpdate;

    /// Spawns a cell by progressively granting `empty_cell` a series of `positions` with [Habitable::grant_position()].
    ///
    /// Assumes that all `positions` are valid.
    fn spawn_cell(
        &mut self,
        empty_cell: EmptyCell<Self::Cell>,
        positions: impl IntoIterator<Item = Pos<usize>>
    ) -> &RelCell<Self::Cell> {
        let cell_index = self.env_mut().cells.add(empty_cell).index;
        let new_spin = Spin::Some(cell_index);
        for pos in positions {
            self.grant_position(pos, new_spin);
        }
        &self.env().cells[cell_index]
    }

    /// Spawns a [Spin::Solid] at each position in `positions`.
    ///
    /// Assumes that all `positions` are valid.
    fn spawn_solid(&mut self, positions: impl Iterator<Item = Pos<usize>>) {
        for pos in positions {
            self.grant_position(pos, Spin::Solid);
        }
    }
}