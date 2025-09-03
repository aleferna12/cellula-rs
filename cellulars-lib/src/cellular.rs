use std::ops::{Deref, DerefMut};
use crate::constants::Spin;
use crate::lattice_entity::LatticeEntity;
use crate::positional::boundary::Boundary;
use crate::positional::pos::Pos;

pub trait Cellular {
    fn target_area(&self) -> u32;
    fn set_target_area(&mut self, value: u32);
    fn area(&self) -> u32;
    fn center(&self) -> Pos<f32>;
    fn shift_position<B: Boundary<Coord = f32>>(&mut self, pos: Pos<usize>, add: bool, bound: &B);
    fn update(&mut self);
    fn birth(&self) -> Self;
    fn die(&mut self);
    fn is_alive(&self) -> bool;
    fn is_valid(&self) -> bool;
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