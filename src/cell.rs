use crate::boundary::LatticeBoundary;
use crate::constants::Spin;
use crate::environment::LatticeEntity;
use crate::pos::{Pos2D, WrappedPos};
use std::ops::{Deref, DerefMut};

/// Represents a cell that is bound to an `Environment`.
///
/// Functions that do not need information about a cell's relational operators 
/// (`spin` and `mom`) should take `&Cell` as an argument instead.
///
/// Implements `Deref<Cell>`.
#[derive(Debug, Clone)]
pub struct RelCell {
    pub spin: Spin,
    pub mom: Spin,
    pub(crate) cell: Cell
}

impl RelCell {
    /// Creates a mock cell with spin and mom = `LatticeEntity<()>::first_cell_spin()` for testing.
    pub fn mock(cell: Cell) -> Self {
        RelCell {
            spin: LatticeEntity::first_cell_spin(),
            mom: LatticeEntity::first_cell_spin(),
            cell
        }
    }
}

impl Deref for RelCell {
    type Target = Cell;

    fn deref(&self) -> &Self::Target {
        &self.cell
    }
}

impl DerefMut for RelCell {
    fn deref_mut(&mut self) -> &mut <Self as Deref>::Target {
        &mut self.cell
    }
}

#[derive(Clone, Debug)]
pub struct Cell {
    pub area: u32,
    pub target_area: u32,
    pub center: WrappedPos
}

impl Cell {
    pub fn new(area: u32, target_area: u32, mut center: WrappedPos) -> Self {
        // Weights projection with current cell area
        center.projection.x_sin *= area as f32;
        center.projection.x_cos *= area as f32;
        center.projection.y_sin *= area as f32;
        center.projection.y_cos *= area as f32;
        Self {
            area,
            target_area,
            center
        }
    }

    pub fn shift_position<B: LatticeBoundary>(
        &mut self,
        pos: Pos2D<usize>,
        width: usize,
        height: usize,
        add: bool
    ) {
        // The order here matters, be careful
        self.shift_area(add);
        B::shift_cell_center(self, pos, width, height, add);
    }

    pub fn shift_area(&mut self, add: bool) {
        if add {
            self.area += 1;
        } else {
            self.area = self.area.saturating_sub(1);
        }
    }
    
}