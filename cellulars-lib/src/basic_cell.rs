use crate::constants::CellIndex;
use crate::positional::boundaries::Boundary;
use crate::positional::pos::Pos;
use std::ops::{Deref, DerefMut};

#[derive(Clone, Debug)]
pub struct BasicCell {
    pub target_area: u32,
    pub area: u32,
    pub center: Pos<f32>,
}

impl BasicCell {
    pub fn new_empty(target_area: u32) -> Self {
        Self {
            target_area,
            area: 0,
            center: Pos::new(0., 0.,)
        }
    }

    pub fn set_target_area(&mut self, value: u32) {
        self.target_area = value;
    }
}

impl Cellular for BasicCell {
    fn target_area(&self) -> u32 {
        self.target_area
    }

    fn area(&self) -> u32 {
        self.area
    }

    fn center(&self) -> Pos<f32> {
        self.center
    }

    // TODO!: This should be type-encoded (CellContainer should be Vec<[State<CellT>]>) where State = {Invalid, Valid}
    fn is_valid(&self) -> bool {
        self.area > 0
    }

    fn shift_position(
        &mut self,
        pos: Pos<usize>,
        add: bool,
        bound: &impl Boundary<Coord = f32>
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
        }
        self.area = self.area.checked_add_signed(shift).expect("overflow in `shift_position`");
    }
}

impl Alive for BasicCell {
    fn is_alive(&self) -> bool {
        self.target_area() > 0
    }

    fn apoptosis(&mut self) {
        self.set_target_area(0)
    }

    fn birth(&self) -> Self {
        let mut newborn = self.clone();
        newborn.area = 0;
        newborn
    }
}

pub trait Cellular {
    fn target_area(&self) -> u32;
    fn area(&self) -> u32;
    fn center(&self) -> Pos<f32>;
    fn is_valid(&self) -> bool;
    fn shift_position(
        &mut self,
        pos: Pos<usize>,
        add: bool,
        bound: &impl Boundary<Coord = f32>
    );
}

// TODO!: mom should not be here and should instead be a symmetric table on ChemEnvironment or Environment
/// Represents a cell that is bound to an `Environment`.
///
/// Functions that do not need information about a cell's relational operators 
/// (`index` and `mom`) should take `&C` directly.
///
/// Implements `Deref<Cell>`.
#[derive(Debug, Clone)]
pub struct RelCell<C> {
    pub index: CellIndex,
    pub cell: C
}

impl<C> RelCell<C> {
    /// Creates a mock cell with index and mom = 0 for testing.
    pub fn mock(cell: C) -> Self {
        RelCell {
            index: 0,
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

pub trait Alive: Cellular {
    fn is_alive(&self) -> bool;
    fn apoptosis(&mut self);
    fn birth(&self) -> Self;
}

/// Shifts a center of mass (`com`) with associated `mass` by `pos`.
pub fn shifted_com<B: Boundary<Coord = f32>>(
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
    let new_com = Pos::new(
        com.x + dx * added_mass / new_mass,
        com.y + dy * added_mass / new_mass,
    );
    // We call this to rewrap the position if necessary
    bound.valid_pos(new_com).expect("shifted the center to outside the available space").into()
}