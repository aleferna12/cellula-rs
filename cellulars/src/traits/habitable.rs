//! Contains logic associated with [`Habitable`].

use crate::environment::{EdgesUpdate, Environment};
use crate::cell_container::RelCell;
use crate::empty_cell::{Empty, EmptyCell};
use crate::positional::boundaries::ToLatticeBoundary;
use crate::positional::pos::Pos;
use crate::prelude::{Cellular, Neighborhood};
use crate::spin::Spin;

/// Types that can be inhabited by cells, since they can be downcast to
/// [`Environment`], where cells live, and know how to [`TransferPosition`] between cells.
pub trait Habitable: TransferPosition + AsEnv {}

impl<H: TransferPosition + AsEnv<Cell = C>, C> Habitable for H {}

/// This trait asserts that a type can spawn cells.
pub trait Spawn: Habitable
where
    Self::Cell: Cellular + Empty {
    /// Spawns a cell by progressively granting `empty_cell` a series of `positions` with [`TransferPosition::transfer_position()`].
    ///
    /// # Panics
    ///
    /// If any position in `positions` is not valid.
    fn spawn_cell(
        &mut self,
        empty_cell: EmptyCell<Self::Cell>,
        positions: impl IntoIterator<Item=Pos<usize>>
    ) -> &RelCell<Self::Cell> {
        let cell_index = self.env_mut().cells.add(empty_cell).index;
        let new_spin = Spin::Some(cell_index);
        for pos in positions {
            self.transfer_position(pos, new_spin);
        }
        &self.env().cells[cell_index]
    }

    /// Spawns a [`Spin::Solid`] at each position in `positions`.
    ///
    /// # Panics
    ///
    /// If any position in `positions` is not valid.
    fn spawn_solid(&mut self, positions: impl Iterator<Item = Pos<usize>>) {
        for pos in positions {
            self.transfer_position(pos, Spin::Solid);
        }
    }
}

/// Types that can transfer ownership of their positions between [`Spin`]s.
pub trait TransferPosition {
    /// Transfers ownership of position `pos` to the entity represented by spin `to`.
    ///
    /// # Panics
    ///
    /// If `pos` is not a valid position in the environment's lattice.
    fn transfer_position(&mut self, pos: Pos<usize>, to: Spin) -> EdgesUpdate;
}

/// Types that can cheaply downcast to a reference to an [`Environment`].
pub trait AsEnv {
    /// Cell type of the environment.
    type Cell;

    /// Returns a reference to the environment where cells live.
    fn env(&self) -> &Environment<Self::Cell, impl Neighborhood, impl ToLatticeBoundary>;

    /// Returns a mutable reference to the environment where cells live.
    fn env_mut(&mut self) -> &mut Environment<Self::Cell, impl Neighborhood, impl ToLatticeBoundary>;
}