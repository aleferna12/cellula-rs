//! Contains logic related to cells.

use crate::constants::CellIndex;
use crate::positional::boundaries::Boundary;
use crate::positional::pos::Pos;
use std::ops::{Deref, DerefMut};

/// Minimum components required to simulate a cell.
#[derive(Clone, Debug)]
pub struct BasicCell {
    /// Cell's current target area.
    pub target_area: u32,
    /// Cell's area.
    area: u32,
    /// Center of mass of the cell.
    center: Pos<f32>,
}

impl BasicCell {
    /// Returns an empty cell to be filled by methods like 
    /// [Habitable::spawn_cell()](crate::habitable::Habitable::spawn_cell())
    pub fn new_empty(target_area: u32) -> Self {
        Self {
            target_area,
            area: 0,
            center: Pos::new(0., 0.,)
        }
    }

    /// Makes a new, ready-to-go cell from a pre-existing state.
    ///
    /// Useful to initialise a cell from a backup.
    /// For most use cases, use [BasicCell::new_empty()] instead.
    pub fn new_ready(
        area: u32,
        center: Pos<f32>,
        target_area: u32
    ) -> Self {
        Self {
            area,
            center,
            target_area
        }
    }

    /// Returns the cell's area.
    pub fn area(&self) -> u32 {
        self.area
    }

    /// Returns the center of mass of the cell.
    pub fn center(&self) -> Pos<f32> {
        self.center
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
        self.is_valid() && self.target_area() > 0
    }

    fn apoptosis(&mut self) {
        self.target_area = 0
    }

    fn birth(&self) -> Self {
        let mut newborn = self.clone();
        newborn.area = 0;
        newborn
    }
}

/// Types that can be used a cell in a simulation.
pub trait Cellular {
    /// Returns the target area of the cell.
    fn target_area(&self) -> u32;
    /// Returns the area of the cell.
    fn area(&self) -> u32;
    /// Returns the center of mass of the cell.
    fn center(&self) -> Pos<f32>;
    // TODO: we should code this property into an enum type
    /// Returns whether the cell is still valid or not.
    /// 
    /// Invalid cells cannot recover from this state, and can effectively be ignored by the simulation algorithm. 
    fn is_valid(&self) -> bool;
    /// Shifts the center and area of the cell by granting (`add == true`) 
    /// or stealing (`add == false`) a position from it.
    fn shift_position(
        &mut self,
        pos: Pos<usize>,
        add: bool,
        bound: &impl Boundary<Coord = f32>
    );
}

/// Represents a cell that is bound to an [Environment](crate::environment::Environment).
///
/// Functions that do not need information about a cell's relational operators 
/// (`index` and `mom`) should take `&C` directly.
///
/// Implements [`Deref<Target = Cell>`].
#[derive(Debug, Clone)]
pub struct RelCell<C> {
    /// Relational cell index that unique to this cell in its 
    /// [Environment](crate::environment::Environment).
    pub index: CellIndex,
    /// Inner cell instance.
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

/// This trait indicates that a [Cellular] can be killed.
pub trait Alive: Cellular {
    /// Returns whether the cell is alive or not.
    fn is_alive(&self) -> bool;
    /// Kills the cell.
    fn apoptosis(&mut self);
    /// Returns a new cell that inherits properties from `self` but is empty and can be filled with 
    /// [Habitable::grant_position()](crate::habitable::Habitable::grant_position).
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