use crate::constants::Spin;
use crate::lattice_entity::LatticeEntity;
use crate::positional::pos::Pos;
use std::ops::{Deref, DerefMut};

pub trait Cellular {
    fn target_area(&self) -> u32;
    fn area(&self) -> u32;
    fn center(&self) -> Pos<f32>;
    fn is_alive(&self) -> bool;
}

/// Represents a cell that is bound to an `Environment`.
///
/// Functions that do not need information about a cell's relational operators 
/// (`spin` and `mom`) should take `&Cell` as an argument instead.
///
/// Implements `Deref<Cell>`.
#[derive(Debug, Clone)]
pub struct RelCell<C> {
    pub spin: Spin,
    pub mom: Spin,
    pub cell: C
}

impl<C> RelCell<C> {
    /// Creates a mock cell with spin and mom = `LatticeEntity<()>::first_cell_spin()` for testing.
    pub fn mock(cell: C) -> Self {
        RelCell {
            spin: LatticeEntity::first_cell_spin(),
            mom: LatticeEntity::first_cell_spin(),
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

pub struct _TestCell {}

impl Cellular for _TestCell {
    fn target_area(&self) -> u32 {
        todo!()
    }

    fn area(&self) -> u32 {
        todo!()
    }

    fn center(&self) -> Pos<f32> {
        todo!()
    }

    fn is_alive(&self) -> bool {
        todo!()
    }
}