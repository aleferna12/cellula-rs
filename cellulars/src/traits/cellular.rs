//! Contains cell traits.

use crate::constants::FloatType;
use crate::empty_cell::EmptyCell;
use crate::positional::boundaries::Boundary;
use crate::positional::com::ShiftError;
use crate::positional::pos::Pos;

/// Types that can be simulated like a cell.
pub trait Cellular {
    /// Returns the target area of the cell.
    fn target_area(&self) -> u32;

    /// Returns the area of the cell.
    fn area(&self) -> u32;

    /// Shifts the area of the cell by adding (`add == true`)
    /// or removing (`add == false`) a position from it.
    fn shift_position(
        &mut self,
        pos: Pos<usize>,
        add: bool,
        bound: &impl Boundary<Coord = FloatType>
    ) -> Result<(), ShiftError>;
}

/// Indicates that this cell keeps track of its center.
pub trait HasCenter {
    /// Returns the center of the cell.
    fn center(&self) -> Pos<FloatType>;
}

/// This trait indicates that a [`Cellular`] can be killed.
pub trait Alive: Cellular + Sized {
    /// Returns whether the cell is alive or not.
    fn is_alive(&self) -> bool;

    /// Kills the cell.
    fn apoptosis(&mut self);

    /// Returns a new cell that inherits properties from `self` but is empty and can be filled with
    /// [`TransferPosition::transfer_position()`](crate::prelude::TransferPosition::transfer_position()).
    fn birth(&self) -> EmptyCell<Self>;
}