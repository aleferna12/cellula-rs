use std::ops::{Deref, DerefMut};
use crate::boundary::LatticeBoundary;
use crate::constants::Spin;
use crate::environment::LatticeEntity;
use crate::pos::{AngularProjection, Pos2D};

/// Represents a cell that is bound to an `Environment`.
///
/// Functions that do not need information about a cell's relational operators 
/// (`spin` and `mom`) should take `Cell` as an argument instead.
///
/// Derefs into `Cell`.
#[derive(Debug, Clone)]
pub struct RelCell {
    pub spin: Spin,
    pub mom: Spin,
    pub(crate) free_cell: Cell
}

impl RelCell {
    /// Creates a mock cell with spin and mom = `LatticeEntity<()>::first_cell_spin()` for testing.
    pub fn mock(free_cell: Cell) -> Self {
        RelCell {
            spin: LatticeEntity::first_cell_spin(),
            mom: LatticeEntity::first_cell_spin(),
            free_cell
        }
    }
}

#[derive(Clone, Debug)]
pub struct Cell {
    pub area: u32,
    pub target_area: u32,
    pub center: CellCenter
}

impl Cell {
    pub fn new(area: u32, target_area: u32, mut center: CellCenter) -> Self {
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

impl Deref for RelCell {
    type Target = Cell;

    fn deref(&self) -> &Self::Target {
        &self.free_cell
    }
}

impl DerefMut for RelCell {
    fn deref_mut(&mut self) -> &mut <Self as Deref>::Target { 
        &mut self.free_cell
    }
}

#[derive(Debug, Clone)]
pub struct CellCenter {
    pub(crate) pos: Pos2D<f32>,
    pub(crate) projection: AngularProjection
}

impl CellCenter {
    pub fn new(pos: Pos2D<f32>, width: usize, height: usize) -> Self {
        Self {
            pos,
            projection: AngularProjection::from_pos(pos, width, height)
        }
    }

    /// Represents the origin of the lattice, at 0, 0.
    pub fn origin() -> Self {
        Self {
            pos: (0., 0.).into(),
            projection: AngularProjection {
                x_sin: 0.,
                x_cos: 1.,
                y_sin: 0.,
                y_cos: 1.,
            }
        }
    }
    
    pub fn pos(&self) -> Pos2D<f32> {
        self.pos
    }
}