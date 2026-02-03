//! Contains logic associated with [`Spin`].

use crate::constants::CellIndex;
#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};
use std::fmt::Debug;

/// This enum represents anything that can own a position in the cell lattice.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash, Default)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub enum Spin {
    /// Represents a position owned by the medium.
    #[default]
    Medium,
    /// Represents a position owned by a solid, immutable object.
    Solid,
    /// Represents a position owned by a cell identified by a [`CellIndex`].
    Some(CellIndex),
}