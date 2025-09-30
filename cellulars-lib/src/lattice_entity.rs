use std::fmt::Debug;
use crate::constants::CellIndex;

pub type Spin = LatticeEntity<CellIndex>;

/// This enum represents anything that can be on the cell lattice.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash)]
pub enum LatticeEntity<C> {
    Solid,
    Medium,
    SomeCell(C),
}

impl<C> Default for LatticeEntity<C> {
    fn default() -> Self {
        LatticeEntity::Medium
    }
}