use crate::constants::CellIndex;
use std::fmt::Debug;

/// This enum represents anything that can be on the cell lattice.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash, Default)]
pub enum Spin {
    #[default]
    Medium,
    Solid,
    Some(CellIndex),
}