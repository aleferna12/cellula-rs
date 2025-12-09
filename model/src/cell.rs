//! Contains logic associated with [Cell].

use cellulars_lib::basic_cell::{shifted_com, Alive, BasicCell, Cellular};
use cellulars_lib::positional::boundaries::Boundary;
use cellulars_lib::positional::pos::Pos;
use std::ops::{Deref, DerefMut};
use bon::Builder;
use strum_macros::{Display, EnumString};

/// A cell that can track a chemical concentration and migrate towards its source.
#[derive(Clone, Debug, Builder)]
pub struct Cell {
    /// Area at which the cell divides.
    pub divide_area: u32,
    /// Target area for newborns of this cell (see [Alive::birth()]).
    pub newborn_target_area: u32,
    /// Current type of the cell.
    pub cell_type: CellType,
    /// Underlying basic cell.
    basic_cell: BasicCell,
    /// Center of the cell weighted by the chemical concentration at each cell position.
    chem_center: Pos<f32>,
    /// Total concentration of the chemical perceived by the cell.
    chem_mass: u32
}

impl Cell {
    /// Initialises an empty [Cell] to be filled progressively with [Cell::shift_position()].
    pub fn new_empty(target_area: u32, divide_area: u32, cell_type: CellType) -> Self {
        Self {
            basic_cell: BasicCell::new_empty(target_area),
            chem_center: Pos::new(0., 0.),
            chem_mass: 0,
            newborn_target_area: target_area,
            divide_area,
            cell_type,
        }
    }
    
    /// Returns the total concentration of the chemical perceived by the cell.
    pub fn chem_mass(&self) -> u32 {
        self.chem_mass
    }

    /// Returns the center of the cell weighted by the chemical concentration at each cell position.
    pub fn chem_center(&self) -> Pos<f32> {
        self.chem_center
    }

    /// Sets the area at which the cell divides when
    /// [MyEnvironment::reproduce()](crate::my_environment::MyEnvironment::reproduce()) is called.
    pub fn set_divide_area(&mut self, value: u32) {
        self.divide_area = value;
    }

    /// Adds or removes the chemical concentration `chem_at` at position `pos` from the cell.
    pub fn shift_chem<B: Boundary<Coord=f32>>(&mut self, pos: Pos<usize>, chem_at: u32, add: bool, bound: &B) {
        let shift = if add { 1 } else { -1 };
        if let Some(new_chem) = shifted_com(
            self.chem_center,
            pos,
            self.chem_mass as f32,
            chem_at as f32,
            shift,
            bound
        ) {
            self.chem_center = new_chem;
        } else {
            self.chem_center = self.center();
        }
        self.chem_mass = self.chem_mass
            .checked_add_signed(shift * chem_at as i32)
            .expect("overflow in `shift_chem`");
    }

    /// Updates parameters of the cell (called by [Pond::step()](cellulars_lib::step::Step::step())).
    pub fn update(&mut self) {
        if let CellType::Dividing = self.cell_type && self.target_area() < self.divide_area {
            let new_target_area = self.target_area() + 1;
            self.target_area = new_target_area;
        }
    }
}

impl Deref for Cell {
    type Target = BasicCell;

    fn deref(&self) -> &Self::Target {
        &self.basic_cell
    }
}

impl DerefMut for Cell {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.basic_cell
    }
}

impl Cellular for Cell {
    fn target_area(&self) -> u32 {
        self.basic_cell.target_area()
    }

    fn area(&self) -> u32 {
        self.basic_cell.area()
    }

    fn center(&self) -> Pos<f32> {
        self.basic_cell.center()
    }

    fn is_valid(&self) -> bool {
        self.basic_cell.is_valid()
    }

    fn shift_position(&mut self, pos: Pos<usize>, add: bool, bound: &impl Boundary<Coord=f32>) {
        self.basic_cell.shift_position(pos, add, bound)
    }
}

impl Alive for Cell {
    fn is_alive(&self) -> bool {
        self.basic_cell.is_alive()
    }

    fn apoptosis(&mut self) {
        self.basic_cell.apoptosis()
    }

    fn birth(&self) -> Self {
        let mut basic_cell = self.basic_cell.birth();
        basic_cell.target_area = self.newborn_target_area;
        Self { 
            basic_cell,
            chem_mass: 0,
            ..self.clone()
        }
    }
}

/// A cell is either migrating or dividing.
#[derive(Clone, Debug, EnumString, Display)]
#[strum(serialize_all = "kebab-case")]
pub enum CellType {
    /// A cell that is migrating.
    Migrating,
    /// A cell that is dividing.
    Dividing
}

#[cfg(test)]
mod tests {
    use super::*;
    use cellulars_lib::positional::boundaries::UnsafePeriodicBoundary;
    use cellulars_lib::positional::rect::Rect;

    fn make_unsafe_boundary() -> UnsafePeriodicBoundary<f32> {
        UnsafePeriodicBoundary::new(Rect::new((0., 0.).into(), (100., 100.).into()))
    }
    
    fn make_test_cell() -> Cell {
        Cell::new_empty(
            100,
            200,
            CellType::Migrating,
        )
    }

    #[test]
    fn test_shift_position_area_and_center() {
        let mut cell = make_test_cell();
        let bound = make_unsafe_boundary();

        cell.shift_position(Pos::new(10, 10), true, &bound);
        assert_eq!(cell.area(), 1);
        assert_eq!(cell.center(), Pos::new(10.0, 10.0));

        cell.shift_position(Pos::new(20, 20), true, &bound);
        assert_eq!(cell.area(), 2);
        assert_eq!(cell.center(), Pos::new(15.0, 15.0));

        cell.shift_position(Pos::new(10, 10), false, &bound);
        assert_eq!(cell.area(), 1);
        assert_eq!(cell.center(), Pos::new(20.0, 20.0));
    }

    #[test]
    fn test_shift_position_chem_center_and_mass() {
        let bound = make_unsafe_boundary();
        let mut cell = make_test_cell();

        // Add chem at (2, 3) with value 10
        cell.shift_chem(Pos::new(2, 3), 10, true, &bound);
        assert_eq!(cell.chem_mass, 10);
        assert_eq!(cell.chem_center, Pos::new(2., 3.));

        // Add chem at (4, 5) with value 10
        cell.shift_chem(Pos::new(4, 5), 10, true, &bound);
        assert_eq!(cell.chem_mass, 20);
        assert_eq!(cell.chem_center, Pos::new(3., 4.));

        // Remove chem from (2, 3)
        cell.shift_chem(Pos::new(2, 3), 10, false, &bound);
        assert_eq!(cell.chem_mass, 10);
        assert_eq!(cell.chem_center, Pos::new(4., 5.));
    }
}
