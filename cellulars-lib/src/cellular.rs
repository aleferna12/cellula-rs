//! Contains cell traits.

use crate::constants::CellIndex;
use crate::positional::boundaries::Boundary;
use crate::positional::pos::Pos;
use std::ops::{Deref, DerefMut};

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

/// Represents a cell that is bound to an [Environment](crate::environment::Environment).
///
/// Functions that do not need information about the cell's `index` relational operators should take 
/// the inner cell type `C` directly.
///
/// Implements [Deref<Target = C>].
#[derive(Debug, Clone)]
pub struct RelCell<C> {
    /// Relational cell index that is unique to this cell in its 
    /// [Environment](crate::environment::Environment).
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

impl<C> Deref for RelCell<C> {
    type Target = C;

    fn deref(&self) -> &Self::Target {
        &self.cell
    }
}

impl<C> DerefMut for RelCell<C> {
    fn deref_mut(&mut self) -> &mut <Self as Deref>::Target {
        &mut self.cell
    }
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