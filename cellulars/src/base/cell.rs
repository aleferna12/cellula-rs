//! Contains logic associated with [`Cell`].

use crate::positional::boundaries::Boundary;
use crate::positional::com::{Com, ShiftError};
use crate::prelude::{Alive, Cellular, FloatType, HasCenter, Pos};
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use crate::empty_cell::{Empty, EmptyCell};

/// Minimum components required to simulate a cell.
///
/// Comparisons with [`PartialEq`] always return true if both cells [`Empty::is_empty()`].
#[derive(Clone, Debug)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Cell {
    /// Cell's current target area.
    pub target_area: u32,
    /// Center of mass of the cell.
    com: Com
}

impl Cell {
    /// Makes a new, ready-to-go cell from a pre-existing state.
    ///
    /// Useful to initialize a cell from a backup.
    /// For most use cases, use [`Cell::new_empty()`] instead.
    pub fn new_ready(
        area: u32,
        target_area: u32,
        center: Pos<FloatType>
    ) -> Self {
        Self {
            com: Com { pos: center, mass: area },
            target_area
        }
    }

    /// Creates a new [`EmptyCell`].
    pub fn new_empty(target_area: u32) -> EmptyCell<Self> {
        EmptyCell::new_unchecked(Self::new_ready(0, target_area, Pos::new(0., 0.)))
    }
}

impl Cellular for Cell {
    fn target_area(&self) -> u32 {
        self.target_area
    }

    fn area(&self) -> u32 {
        self.com.mass
    }

    fn shift_position(
        &mut self,
        pos: Pos<usize>,
        adding: bool,
        boundary: &impl Boundary<Coord = FloatType>
    ) -> Result<(), ShiftError> {
        let shifted = self.com.shift(
            Com { pos: pos.cast_as(), mass: 1 },
            adding,
            boundary
        );
        shifted.map(|new_com| {
            self.com = new_com;
        })
    }
}

impl HasCenter for Cell {
    fn center(&self) -> Pos<FloatType> {
        self.com.pos
    }
}

impl Alive for Cell {
    fn is_alive(&self) -> bool {
        !self.is_empty() && self.target_area() > 0
    }

    fn apoptosis(&mut self) {
        self.target_area = 0
    }

    fn birth(&self) -> EmptyCell<Cell> {
        let mut newborn = self.clone();
        newborn.com.mass = 0;
        EmptyCell::new(newborn).expect("cell is not empty")
    }
}

impl PartialEq for Cell {
    fn eq(&self, other: &Self) -> bool {
        if self.is_empty() && other.is_empty() {
            return true;
        }
        self.target_area == other.target_area
            && self.com == other.com
    }
}

impl Empty for Cell {
    fn empty_default() -> EmptyCell<Self> {
        EmptyCell::new_unchecked(Self {
            target_area: 0,
            com: Com { pos: Pos::new(0., 0.), mass: 0 }
        })
    }

    fn is_empty(&self) -> bool {
        self.area() == 0
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_cell() {
        let empty_cell = Cell::empty_default();
        let cell = empty_cell.as_cell();
        assert!(cell.is_empty());
        assert_eq!(cell.area(), 0);
        assert!(!cell.is_alive())
    }
}