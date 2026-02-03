//! Contains logic associated with [CellContainer].

use crate::constants::CellIndex;
use crate::traits::cellular::{Alive, Cellular, EmptyCell};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::ops::{Index, IndexMut};

/// This is a vector type containing cell instances that can be accessed with their respective unique [CellIndex]es.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
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

    /// Constructs a new, empty [`CellContainer`] with at least the specified `capacity`.
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            vec: Vec::with_capacity(capacity)
        }
    }
    
    /// Returns the total number of cells in the cell container, including empty cells (see [Cellular::is_empty()]).
    pub fn n_cells(&self) -> CellIndex {
        self.vec.len().try_into().expect("there are more cells than supported by the type `CellIndex`")
    }
    
    /// Replaces the cell at `rel_cell.index` with `rel_cell`.
    pub fn replace(&mut self, rel_cell: RelCell<C>) -> RelCell<C> {
        std::mem::replace(&mut self.vec[rel_cell.index as usize], rel_cell)
    }

    /// Removes all cells from the cell container, returning it to a clean-slate state.
    pub fn wipe_out(&mut self) {
        self.vec.clear()
    }

    /// Returns an iterator of all cells (including empty cells).
    pub fn iter(&self) -> impl Iterator<Item = &RelCell<C>> {
        self.vec.iter()
    }

    /// Returns an iterator of mutable references to all cells (including empty cells).
    pub fn iter_mut(&mut self) -> impl Iterator<Item=&mut RelCell<C>> {
        self.vec.iter_mut()
    }
}

impl<C: Cellular> CellContainer<C> {
    /// Returns the number of cells that are not empty (see [Cellular::is_empty()]).
    pub fn n_non_empty(&self) -> CellIndex {
        self.vec
            .iter()
            .filter(|rel_cell| !rel_cell.cell.is_empty())
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
            .find(|rel_cell| rel_cell.cell.is_empty())
            .map(|rel_cell| rel_cell.index)
            .unwrap_or(self.n_cells())
    }

    /// Add a cell by replacing the first empty cell (see [Cellular::is_empty()]).
    pub fn add(&mut self, cell: EmptyCell<C>) -> &mut RelCell<C> {
        let new_index = self.next_index();
        let rel_cell = RelCell {
            index: new_index,
            cell: cell.into_cell()
        };

        if new_index == self.n_cells() {
            self.vec.push(rel_cell);
        } else {
            self.replace(rel_cell);
        }
        &mut self.vec[new_index as usize]
    }

    /// Add a cell to the end of the cell container without replacing any empty cells.
    pub fn push(&mut self, cell: EmptyCell<C>) -> &mut RelCell<C> {
        let new_index = self.n_cells();
        let rel_cell = RelCell {
            index: new_index,
            cell: cell.into_cell()
        };
        self.vec.push(rel_cell);
        &mut self.vec[new_index as usize]
    }

    /// Returns an iterator of non-empty cells.
    pub fn iter_non_empty(&self) -> impl Iterator<Item = &RelCell<C>> {
        self.iter().filter(|rel_cell| !rel_cell.cell.is_empty())
    }

    /// Returns an iterator of mutable references to non-empty cells.
    pub fn iter_non_empty_mut(&mut self) -> impl Iterator<Item = &mut RelCell<C>> {
        self.iter_mut().filter(|rel_cell| !rel_cell.cell.is_empty())
    }

    /// Gets a reference to a cell using its unique cell index.
    ///
    /// Returns [None] if the index points to an empty cell.
    pub fn get(&self, index: CellIndex) -> Option<&RelCell<C>> {
        self.vec.get(index as usize).filter(|rel_cell| !rel_cell.cell.is_empty())
    }

    /// Gets a mutable reference to a cell using its unique cell index.
    ///
    /// Returns [None] if the index points to an empty cell.
    pub fn get_mut(&mut self, index: CellIndex) -> Option<&mut RelCell<C>> {
        self.vec.get_mut(index as usize).filter(|rel_cell| !rel_cell.cell.is_empty())
    }
}

impl<C> Default for CellContainer<C> {
    fn default() -> Self {
        Self::new()
    }
}

impl<C> IntoIterator for CellContainer<C> {
    type Item = RelCell<C>;
    type IntoIter = std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.vec.into_iter()
    }
}

impl<'c, C> IntoIterator for &'c CellContainer<C> {
    type Item = &'c RelCell<C>;
    type IntoIter = std::slice::Iter<'c, RelCell<C>>;

    fn into_iter(self) -> Self::IntoIter {
        self.vec.iter()
    }
}

impl<'c, C> IntoIterator for &'c mut CellContainer<C> {
    type Item = &'c mut RelCell<C>;
    type IntoIter = std::slice::IterMut<'c, RelCell<C>>;

    fn into_iter(self) -> Self::IntoIter {
        self.vec.iter_mut()
    }
}

impl<C> Index<CellIndex> for CellContainer<C> {
    type Output = RelCell<C>;

    fn index(&self, index: CellIndex) -> &Self::Output {
        &self.vec[index as usize]
    }
}

impl<C> IndexMut<CellIndex> for CellContainer<C> {
    fn index_mut(&mut self, index: CellIndex) -> &mut Self::Output {
        &mut self.vec[index as usize]
    }
}

/// Represents a cell that is bound to an [Environment](crate::base::environment::Environment).
///
/// Functions that do not need information about the cell's `index` relational operators should take
/// the inner cell type `C` directly.
///
/// Implements [Deref<Target = C>].
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct RelCell<C> {
    /// Relational cell index that is unique to this cell in its
    /// [Environment](crate::base::environment::Environment).
    pub index: CellIndex,
    /// Inner cell instance.
    pub cell: C
}

/// Initializes a [`CellContainer`] with [`EmptyCell`]s.
///
/// Syntax is the same as `vec![]`.
///
/// # Examples
///
/// ```
/// use cellulars_lib::cell_container;
/// use cellulars_lib::prelude::{Cell, CellContainer};
///
/// let cc1 = cell_container![Cell::new_empty(10), Cell::new_empty(10)];
/// let cc2 = cell_container![Cell::new_empty(10); 2];
/// assert_eq!(cc1, cc2);
/// ```
#[macro_export]
macro_rules! cell_container {
    () => {{
        CellContainer::new()
    }};
    ( $($cell:expr),+ $(,)? ) => {{
        let mut cc = CellContainer::new();
        $(
            cc.push($cell);
        )*
        cc
    }};
    ( $cell:expr; $number:expr ) => {{
        let mut cc = CellContainer::with_capacity($number);
        for _ in 0..$number {
            cc.push($cell.clone());
        }
        cc
    }};
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::prelude::Cell;

    fn cell_container() -> CellContainer<Cell> {
        let mut cc = cell_container![Cell::new_empty(10); 2];
        cc.replace(RelCell{ index: 0, cell: Cell::new_ready(10, 10, (0., 0.).into()) });
        cc
    }

    #[test]
    fn test_get() {
        let cc = cell_container();
        assert!(cc.get(0).is_some());
        assert!(cc.get(1).is_none());
    }

    #[test]
    fn test_index() {
        let cc = cell_container();
        assert!(!cc[0].cell.is_empty());
        assert!(cc[1].cell.is_empty());
    }

    #[test]
    fn test_add() {
        let mut cc = cell_container();
        let cell = cc.add(Cell::new_empty(10));
        assert_eq!(1, cell.index);
        assert_eq!(cc.n_cells(), 2);
    }

    #[test]
    fn test_push() {
        let mut cc = cell_container();
        cc.push(Cell::new_empty(10));
        assert_eq!(cc.n_cells(), 3);
        assert!(cc.get(2).is_none());
    }

    #[test]
    fn test_replace() {
        let mut cc = cell_container();
        let cell = cc.replace(RelCell{ index: 0, cell: Cell::new_empty(10).into_cell() });
        assert_eq!(0, cell.index);
        assert_eq!(cc.n_cells(), 2);
        assert!(cc.get(0).is_none());
    }
}
