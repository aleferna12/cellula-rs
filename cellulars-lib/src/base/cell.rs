//! Contains logic associated with [`Cell`].

use crate::positional::boundaries::Boundary;
use crate::positional::com::{Com, ShiftError};
use crate::prelude::{Alive, Cellular, FloatType, HasCenter, Pos};
use crate::traits::cellular::EmptyCell;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Minimum components required to simulate a cell.
#[derive(Clone, Debug, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Cell {
    /// Cell's current target area.
    pub target_area: u32,
    /// Center of mass of the cell.
    com: Com
}

impl Cell {
    /// Returns an empty cell to be filled by methods like 
    /// [`Habitable::spawn_cell()`](crate::traits::habitable::Habitable::spawn_cell()).
    pub fn new_empty(target_area: u32) -> EmptyCell<Self> {
        EmptyCell::new(Self {
            target_area,
            com: Com { pos: Pos::new(0., 0.), mass: 0 }
        }).unwrap()
    }

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
}

impl Cellular for Cell {
    fn target_area(&self) -> u32 {
        self.target_area
    }

    fn area(&self) -> u32 {
        self.com.mass
    }

    // Experimented with encoding this using typestate pattern but it was not helpful nor ergonomic
    fn is_empty(&self) -> bool {
        self.com.mass == 0
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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_empty_cell() {
        let empty_cell = Cell::new_empty(0);
        let cell = empty_cell.as_cell();
        assert!(cell.is_empty());
        assert_eq!(cell.area(), 0);
        assert!(!cell.is_alive())
    }
}