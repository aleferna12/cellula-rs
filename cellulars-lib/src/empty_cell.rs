use crate::prelude::Cellular;

/// A cell who is guaranteed to be empty (see [`Cellular::is_empty()`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct EmptyCell<C>(C);

impl<C> EmptyCell<C> {
    /// Returns the inner cell, which is guaranteed to be [Cellular::`is_empty()`].
    pub fn into_cell(self) -> C {
        self.0
    }

    /// Returns a reference to the inner cell, which is guaranteed to be [`Cellular::is_empty()`].
    pub fn as_cell(&self) -> &C {
        &self.0
    }
    
    pub fn new_unchecked(cell: C) -> Self {
        EmptyCell(cell)
    }
}

impl<C: Empty> EmptyCell<C>
where
    C: Cellular {
    /// Returns `Some(cell)` if `cell` is [`Cellular::is_empty()`] and [`None`] otherwise.
    pub fn new(cell: C) -> Option<Self> {
        if cell.is_empty() {
            return Some(EmptyCell(cell))
        }
        None
    }
}

pub trait Empty where Self: Sized {
    fn empty_default() -> EmptyCell<Self>;
    fn is_empty(&self) -> bool;
}