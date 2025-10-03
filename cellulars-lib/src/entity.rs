use crate::constants::CellIndex;
use std::fmt::Debug;

pub type Spin = Entity<CellIndex>;

/// This enum represents anything that can be on the cell lattice.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash, Default)]
pub enum Entity<C> {
    #[default]
    Medium,
    Solid,
    Some(C),
}

impl<C> Entity<C> {
    pub fn map<D, F: FnOnce(C) -> D>(self, f: F) -> Entity<D> {
        match self {
            Entity::Some(cell) => Entity::Some(f(cell)),
            Entity::Medium => Entity::Medium,
            Entity::Solid => Entity::Solid,
        }
    }
}