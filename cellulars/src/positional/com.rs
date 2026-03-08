//! Contains logic associated with [`Com`].

use crate::constants::FloatType;
use crate::positional::boundaries::Boundary;
use crate::prelude::Pos;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// A center of mass of a cell that is shifted throughout the simulation.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct Com {
    /// Position of the [`Com`].
    pub pos: Pos<FloatType>,
    /// Mass of the [`Com`].
    pub mass: u32,
}

impl Com {
    /// Shifts this [`Com`] by `other`, weighting their relative masses.
    ///
    /// Set `adding` to `false` to instead do the inverse operation,
    /// removing the relative influence of `other` from `self`.
    ///
    /// This operation is prone to accumulating floating-point errors, so be careful with those.
    pub fn shift(
        &self,
        other: Com,
        adding: bool,
        bound: &impl Boundary<Coord = FloatType>
    ) -> Result<Com, ShiftError> {
        let shift = if adding { 1 } else { -1 };
        let added_mass = shift * other.mass as i32;
        let Some(new_mass) = self.mass.checked_add_signed(added_mass) else {
            return Err(ShiftError::NegativeMass(self.mass as i32 + added_mass));
        };
        if new_mass == 0 {
            return Ok(Com { pos: self.pos, mass: new_mass });
        }

        let (dx, dy) = bound.displacement(self.pos, other.pos);
        let new_pos = Pos::new(
            self.pos.x + dx * added_mass as FloatType / new_mass as FloatType,
            self.pos.y + dy * added_mass as FloatType / new_mass as FloatType,
        );
        let valid_pos = bound.valid_pos(new_pos);
        match valid_pos {
            Some(pos) => Ok(Com { pos, mass: new_mass }),
            None => Err(ShiftError::OutOfBounds(new_pos))
        }
    }
}

#[derive(thiserror::Error, Debug)]
/// Error thrown when shifting a position fails.
pub enum ShiftError {
    /// Shifting resulted in a negative mass.
    #[error("shifted COM has negative mass {0}")]
    NegativeMass(i32),
    /// Shifting resulted in position out of bounds.
    #[error("shifted COM `{0:?}` is out of bounds")]
    OutOfBounds(Pos<FloatType>),
}