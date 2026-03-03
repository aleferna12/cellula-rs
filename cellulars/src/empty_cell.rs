//! Contains logic associated with [`EmptyCell`].

use crate::prelude::Cellular;

/// A cell who is guaranteed to be empty.
///
/// Empty cells own no sites in a lattice and can be added to the
/// [`CellContainer`](crate::prelude::CellContainer) of an [`Environment`](crate::prelude::Environment).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct EmptyCell<C>(C);

impl<C> EmptyCell<C> {
    /// Returns the inner cell, which is guaranteed to be [Empty::`is_empty()`].
    pub fn into_cell(self) -> C {
        self.0
    }

    /// Returns a reference to the inner cell, which is guaranteed to be [`Empty::is_empty()`].
    pub fn as_cell(&self) -> &C {
        &self.0
    }

    /// Wraps the `cell` in an [`EmptyCell`] without checking if the cell is actually empty.
    pub fn new_unchecked(cell: C) -> Self {
        EmptyCell(cell)
    }
}

impl<C: Empty> EmptyCell<C>
where
    C: Cellular {
    /// Returns `Some(cell)` if `cell` is [`Empty::is_empty()`] and [`None`] otherwise.
    pub fn new(cell: C) -> Option<Self> {
        if cell.is_empty() {
            return Some(EmptyCell(cell))
        }
        None
    }
}

/// Cell types that can be empty, which makes them compatible with [`CellContainer`](crate::prelude::CellContainer).
pub trait Empty where Self: Sized {
    /// Returns a default empty cell which is used by [`CellContainer`](crate::prelude::CellContainer)s constructor.
    fn empty_default() -> EmptyCell<Self>;
    /// Returns whether this cell is empty or not.
    fn is_empty(&self) -> bool;
}