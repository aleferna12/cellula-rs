//! Contains cell traits.

use crate::positional::boundaries::Boundary;
use crate::positional::pos::Pos;

/// Types that can be used a cell in a simulation.
pub trait Cellular {
    /// Returns the target area of the cell.
    fn target_area(&self) -> u32;
    /// Returns the area of the cell.
    fn area(&self) -> u32;
    /// Returns the center of mass of the cell.
    fn center(&self) -> Pos<f32>;
    /// Returns whether the cell is still valid or not.
    ///
    /// Invalid cells cannot recover from this state, and can effectively be ignored by the simulation algorithm. 
    fn is_valid(&self) -> bool;
    /// Shifts the center and area of the cell by granting (`add == true`) 
    /// or stealing (`add == false`) a position from it.
    fn shift_position(
        &mut self,
        pos: Pos<usize>,
        add: bool,
        bound: &impl Boundary<Coord = f32>
    );
}

/// This trait indicates that a [Cellular] can be killed.
pub trait Alive: Cellular {
    /// Returns whether the cell is alive or not.
    fn is_alive(&self) -> bool;
    /// Kills the cell.
    fn apoptosis(&mut self);
    /// Returns a new cell that inherits properties from `self` but is empty and can be filled with 
    /// [Habitable::grant_position()](crate::habitable::Habitable::grant_position).
    fn birth(&self) -> Self;
}