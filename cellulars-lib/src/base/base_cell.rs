//! Contains logic associated with [BaseCell].

use crate::positional::boundaries::Boundary;
use crate::positional::pos::Pos;
use crate::traits::cellular::{Alive, Cellular};
use thiserror::Error;

/// Minimum components required to simulate a cell.
#[derive(Clone, Debug)]
pub struct BaseCell {
    /// Cell's current target area.
    pub target_area: u32,
    /// Cell's area.
    area: u32,
    /// Center of mass of the cell.
    center: Pos<f32>,
}

impl BaseCell {
    /// Returns an empty cell to be filled by methods like 
    /// [Habitable::spawn_cell()](crate::traits::habitable::Habitable::spawn_cell())
    pub fn new_empty(target_area: u32) -> Self {
        Self {
            target_area,
            area: 0,
            center: Pos::new(0., 0.,)
        }
    }

    /// Makes a new, ready-to-go cell from a pre-existing state.
    ///
    /// Useful to initialize a cell from a backup.
    /// For most use cases, use [BaseCell::new_empty()] instead.
    pub fn new_ready(
        area: u32,
        center: Pos<f32>,
        target_area: u32
    ) -> Self {
        Self {
            area,
            center,
            target_area
        }
    }

    /// Returns the cell's area.
    pub fn area(&self) -> u32 {
        self.area
    }

    /// Returns the center of mass of the cell.
    pub fn center(&self) -> Pos<f32> {
        self.center
    }
}

impl Cellular for BaseCell {
    fn target_area(&self) -> u32 {
        self.target_area
    }

    fn area(&self) -> u32 {
        self.area
    }

    fn center(&self) -> Pos<f32> {
        self.center
    }

    // Experimented with encoding this using typestate pattern but it was not helpful nor ergonomic
    fn is_valid(&self) -> bool {
        self.area > 0
    }

    fn shift_position(
        &mut self,
        pos: Pos<usize>,
        add: bool,
        bound: &impl Boundary<Coord = f32>
    ) {
        let shift = if add { 1 } else { -1 };
        // The order here matters (area is last), be careful
        let shifted = shifted_com(
            self.center,
            pos,
            self.area as f32,
            1.,
            shift,
            bound
        );
        match shifted {
            Ok(new_center) => self.center = new_center,
            Err(e) => log::warn!("Failed to shift center of mass: {e}")
        }
        self.area = self.area.checked_add_signed(shift).expect("overflow in `shift_position`");
    }
}

impl Alive for BaseCell {
    fn is_alive(&self) -> bool {
        self.is_valid() && self.target_area() > 0
    }

    fn apoptosis(&mut self) {
        self.target_area = 0
    }

    fn birth(&self) -> Self {
        let mut newborn = self.clone();
        newborn.area = 0;
        newborn
    }
}

/// Shifts a center of mass (`com`) with associated `mass` by `pos`.
pub fn shifted_com<B: Boundary<Coord = f32>>(
    com: Pos<f32>,
    pos: Pos<usize>,
    com_mass: f32,
    pos_mass: f32,
    shift: i32,
    bound: &B
) -> Result<Pos<f32>, ShiftError> {
    let shift = shift as f32;
    let added_mass = shift * pos_mass;
    let new_mass = com_mass + added_mass;
    if new_mass == 0. {
        return Ok(com)
    } else if new_mass < 0. {
        return Err(ShiftError::NegativeMass(new_mass));
    }
    let (dx, dy) = bound.displacement(com, Pos::new(pos.x as f32, pos.y as f32));
    let new_com = Pos::new(
        com.x + dx * added_mass / new_mass,
        com.y + dy * added_mass / new_mass,
    );
    // We call this to rewrap the position if necessary
    bound.valid_pos(new_com).ok_or(ShiftError::OutOfBounds(new_com))
}

#[derive(Error, Debug)]
/// Error thrown when a [shifted_com()] operation fails.
pub enum ShiftError {
    /// Shifting resulted in a negative mass.
    #[error("shifted COM has negative mass {0}")]
    NegativeMass(f32),
    /// Shifting resulted in position out of bounds.
    #[error("shifted COM `{0:?}` is out of bounds")]
    OutOfBounds(Pos<f32>),
}