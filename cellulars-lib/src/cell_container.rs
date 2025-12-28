//! Contains logic associated with [CellContainer].

use crate::constants::CellIndex;
use crate::traits::cellular::{Alive, Cellular};

/// This is a vector type containing cell instances that can be accessed with their respective unique [CellIndex]es.
#[derive(Clone, Debug, PartialEq)]
pub struct CellContainer<C> {
    vec: Vec<RelCell<C>>
}

impl<C> CellContainer<C> {
    /// Creates an empty cell container.
    pub fn new() -> Self {
        Self {
            vec: vec![],
        }
    }
    
    /// Returns the total number of cells in the cell container, including invalid cells (see [Cellular::is_valid()]).
    pub fn n_cells(&self) -> CellIndex {
        self.vec.len().try_into().expect("there are more cells than supported by the type `CellIndex`")
    }
    
    /// Replaces the cell at `rel_cell.index` with `rel_cell`.
    pub fn replace(&mut self, rel_cell: RelCell<C>) -> RelCell<C> {
        std::mem::replace(&mut self.vec[rel_cell.index as usize], rel_cell)
    }

    /// Returns a reference to a cell using its unique cell index.
    /// 
    /// The cell might be invalid (see [Cellular::is_valid()]).
    pub fn get_cell(&self, index: CellIndex) -> &RelCell<C> {
        &self.vec[index as usize]
    }

    /// Returns a mutable reference to a cell using its unique cell index.
    ///
    /// The cell might be invalid (see [Cellular::is_valid()]).
    pub fn get_cell_mut(&mut self, index: CellIndex) -> &mut RelCell<C> {
        &mut self.vec[index as usize]
    }

    /// Removes all cells from the ccell container, returning it to a clean-slate state.
    pub fn wipe_out(&mut self) {
        self.vec.clear()
    }

    /// Returns an iterator to references of all cells (including invalid).
    pub fn iter(&self) -> impl Iterator<Item=&RelCell<C>> {
        self.vec.iter()
    }

    /// Returns an iterator to mutable references of all cells (including invalid).
    pub fn iter_mut(&mut self) -> impl Iterator<Item=&mut RelCell<C>> {
        self.vec.iter_mut()
    }
}

impl<C: Cellular> CellContainer<C> {
    /// Returns the number of cells that are in a valid state (see [Cellular::is_valid()]).
    pub fn n_valid(&self) -> CellIndex {
        self.vec
            .iter()
            .filter(|rel_cell| rel_cell.cell.is_valid())
            .count() as CellIndex
    }
    
    /// Returns the number of cells that are alive (see [Alive::is_alive()]).
    pub fn n_alive(&self) -> CellIndex
    where C: Alive {
        self.vec
            .iter()
            .filter(|rel_cell| rel_cell.cell.is_alive())
            .count() as CellIndex
    }

    fn next_index(&self) -> CellIndex {
        self.vec
            .iter()
            .find(|rel_cell| !rel_cell.cell.is_valid())
            .map(|rel_cell| rel_cell.index)
            .unwrap_or(self.n_cells())
    }

    /// Add a cell by replacing the first invalid cell (see [Cellular::is_valid()]).
    pub fn add(&mut self, cell: C) -> &mut RelCell<C> {
        let new_index = self.next_index();
        let rel_cell = RelCell {
            index: new_index,
            cell
        };

        if new_index == self.n_cells() {
            self.vec.push(rel_cell);
        } else {
            self.replace(rel_cell);
        }
        &mut self.vec[new_index as usize]
    }

    /// Add a cell to the end of the cell container without replacing any invalid cells. 
    pub fn push(&mut self, cell: C) -> &mut RelCell<C> {
        let new_index = self.n_cells();
        let rel_cell = RelCell {
            index: new_index,
            cell
        };
        self.vec.push(rel_cell);
        &mut self.vec[new_index as usize]
    }
}

impl<C> Default for CellContainer<C> {
    fn default() -> Self {
        Self::new()
    }
}

/// Represents a cell that is bound to an [Environment](crate::base::base_environment::BaseEnvironment).
///
/// Functions that do not need information about the cell's `index` relational operators should take
/// the inner cell type `C` directly.
///
/// Implements [Deref<Target = C>].
#[derive(Clone, Debug, PartialEq)]
pub struct RelCell<C> {
    /// Relational cell index that is unique to this cell in its
    /// [Environment](crate::base::base_environment::BaseEnvironment).
    pub index: CellIndex,
    /// Inner cell instance.
    pub cell: C
}

impl<C> RelCell<C> {
    /// Creates a mock cell with index and mom = 0 for testing.
    pub fn mock(cell: C) -> Self {
        RelCell {
            index: 0,
            cell
        }
    }
}
