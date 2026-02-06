//! Contains cell traits.

use crate::constants::FloatType;
use crate::positional::boundaries::Boundary;
use crate::positional::com::ShiftError;
use crate::positional::pos::Pos;

/// Types that can be simulated like a cell.
pub trait Cellular {
    /// Returns the target area of the cell.
    fn target_area(&self) -> u32;
    /// Returns the area of the cell.
    fn area(&self) -> u32;
    /// Returns whether the cell is empty or not.
    ///
    /// Empty cells cannot recover from this state, and can effectively be ignored by the simulation algorithm.
    ///
    /// A cell that has been validated to be empty is an [`EmptyCell`].
    fn is_empty(&self) -> bool;
    /// Shifts the area of the cell by adding (`add == true`)
    /// or removing (`add == false`) a position from it.
    fn shift_position(
        &mut self,
        pos: Pos<usize>,
        add: bool,
        bound: &impl Boundary<Coord = FloatType>
    ) -> Result<(), ShiftError>;
}

pub trait HasCenter {
    fn center(&self) -> Pos<FloatType>;
}

/// This trait indicates that a [`Cellular`] can be killed.
pub trait Alive: Cellular + Sized {
    /// Returns whether the cell is alive or not.
    fn is_alive(&self) -> bool;
    /// Kills the cell.
    fn apoptosis(&mut self);
    /// Returns a new cell that inherits properties from `self` but is empty and can be filled with
    /// [`Habitable::grant_position()`](crate::traits::habitable::Habitable::grant_position).
    fn birth(&self) -> EmptyCell<Self>;
}

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
}

impl<C> EmptyCell<C>
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