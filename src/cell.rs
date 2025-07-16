use crate::constants::Spin;
use crate::environment::LatticeEntity;
use crate::positional::boundary::Boundary;
use crate::positional::pos::{AngularProjection, Pos, WrappedPos};
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
    pub fn new(area: u32, target_area: u32, center: WrappedPos) -> Self {
        Self {
            area,
            target_area,
            center
        }
    }

    pub fn shift_position<B: Boundary<Coord = f32>>(
        &mut self,
        pos: Pos<usize>,
        add: bool,
        bound: &B
    ) {
        // The order here matters, be careful
        self.shift_center(pos, add, bound);
        self.shift_area(add);
    }

    pub fn shift_area(&mut self, add: bool) {
        if add {
            self.area += 1;
        } else {
            self.area = self.area.saturating_sub(1);
        }
    }

    pub fn shift_center<B: Boundary<Coord = f32>>(&mut self, pos: Pos<usize>, add: bool, bound: &B) -> bool {
        let shift = if add { 1. } else { -1. };
        let new_mass = self.area as f32 + shift;
        if new_mass <= 0.0 {
            return false;
        }

        let center_pos = self.center.pos;
        let (dx, dy) = bound.displacement(center_pos, Pos::new(pos.x as f32, pos.y as f32));
        // TODO! work this into Boundary
        let x = (center_pos.x + dx * shift / new_mass).rem_euclid(bound.rect().width());
        let y = (center_pos.y + dy * shift / new_mass).rem_euclid(bound.rect().height());
        let new_center = Pos::new(x, y);
        self.center = WrappedPos::new(new_center, AngularProjection::from_pos(
            new_center,
            bound.rect().width() as usize,
            bound.rect().height() as usize
        ));
        true
    }
}