use crate::basic_cell::RelCell;
use crate::constants::Spin;
use crate::lattice_entity::LatticeEntity::*;
use std::fmt::Debug;

/// This enum represents anything that can be on the cell lattice.
#[derive(Debug, Copy, Clone)]
pub enum LatticeEntity<C> {
    Solid,
    Medium,
    SomeCell(C),
}

impl<C> LatticeEntity<C> {
    pub fn map<D, F: FnOnce(C) -> D>(self, f: F) -> LatticeEntity<D> {
        match self {
            SomeCell(c) => SomeCell(f(c)),
            Medium => Medium,
            Solid => Solid,
        }
    }
}

impl<C> LatticeEntity<&RelCell<C>> {
    /// Maps the `LatticeEntity` to a unique spin value.
    ///
    /// If you need the spin of a non-cell variant use `discriminant()` instead.
    pub fn spin(&self) -> Spin {
        match self {
            SomeCell(cell) => cell.spin,
            Medium => Medium.discriminant(),
            Solid => Solid.discriminant()
        }
    }
}

impl<C> LatticeEntity<C> {
    pub fn unwrap_cell(self) -> C
    where C: Debug {
        match self {
            SomeCell(cell) => cell,
            _ => panic!("called `LatticeEntity::unwrap_cell()` on a `{self:?}` value")
        }
    }

    pub fn expect_cell(self, message: &str) -> C {
        match self {
            SomeCell(cell) => cell,
            _ => panic!("{}", message)
        }
    }
}

impl LatticeEntity<()> {
    // TODO!: What if this is just Enum::length and then we subtract one for each discriminant?
    /// Returns the first spin that corresponds to a cell.
    ///
    /// This is required to be larger than `Medium.discriminant()` and `Solid.discriminant()`.
    pub fn first_cell_spin() -> Spin {
        2
    }

    pub fn discriminant(&self) -> Spin {
        match self {
            SomeCell(_) => Self::first_cell_spin(),
            Medium => 0,
            Solid => 1
        }
    }
}