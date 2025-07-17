use crate::constants::Spin;
use crate::environment::LatticeEntity;
use crate::positional::boundary::Boundary;
use crate::positional::pos::Pos;
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
    pub center: Pos<f32>,
    pub light_center: Pos<f32>,
    pub light_mass: u32
}

impl Cell {
    /// Initialises an empty `Cell` to be filled progressively with `shift_position()`.
    pub fn new(target_area: u32) -> Self {
        Self {
            area: 0,
            target_area,
            center: Pos::new(0., 0.),
            light_center: Pos::new(0., 0.),
            light_mass: 0
        }
    }

    pub fn shift_position<B: Boundary<Coord = f32>>(
        &mut self,
        pos: Pos<usize>,
        light_at_pos: u32,
        add: bool,
        bound: &B
    ) {
        let shift = if add { 1 } else { -1 };
        // The order here matters (area is last), be careful
        if let Some(new_center) = shifted_com(
            self.center,
            pos, 
            self.area as f32,
            1.,
            shift, 
            bound
        ) {
            self.center = new_center;
            self.area = self.area.saturating_add_signed(shift);
        }
        if let Some(new_chem) = shifted_com(
            self.light_center,
            pos,
            self.light_mass as f32,
            light_at_pos as f32,
            shift,
            bound
        ) {
            self.light_center = new_chem;
            self.light_mass = self.light_mass.saturating_add_signed(shift * light_at_pos as i32);
        }
    }
}

/// Shifts a center of mass (`com`) with associated `mass` by `pos`.
fn shifted_com<B: Boundary<Coord = f32>>(
    com: Pos<f32>,
    pos: Pos<usize>,
    com_mass: f32,
    pos_mass: f32,
    shift: i32,
    bound: &B
) -> Option<Pos<f32>> {
    let shift = shift as f32;
    let added_mass = shift * pos_mass;
    let new_mass = com_mass + added_mass;
    if new_mass <= 0. { 
        return None;
    }

    let (dx, dy) = bound.displacement(com, Pos::new(pos.x as f32, pos.y as f32));
    // We call this to rewrap the position if necessary
    let shifted_pos = Pos::new(com.x + dx, com.y + dy);
    let new_com = Pos::new(
        (com.x * com_mass + shifted_pos.x * added_mass) / new_mass,
        (com.y * com_mass + shifted_pos.y * added_mass) / new_mass,
    );
    bound.valid_pos(new_com).expect("Shifted the center to outside the available space").into()
}