//! Contains logic associated with [`MyCell`].

use bon::Builder;
use cellulars_lib::base::cell::Cell;
use cellulars_lib::constants::FloatType;
use cellulars_lib::positional::boundaries::Boundary;
use cellulars_lib::positional::com::{Com, ShiftError};
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::traits::cellular::{Alive, Cellular, EmptyCell, HasCenter};
use strum_macros::{Display, EnumString};

/// A cell that can track a chemical concentration and migrate towards its source.
#[derive(Clone, Debug, Builder)]
pub struct MyCell {
    /// Area at which the cell divides.
    pub divide_area: u32,
    /// Target area for newborns of this cell (see [`Alive::birth()`]).
    pub newborn_target_area: u32,
    /// Current type of the cell.
    pub cell_type: CellType,
    /// Inner base cell.
    pub cell: Cell,
    /// Center of mass of the cell's perceived chemical concentration.
    chem_com: Com

}

impl MyCell {
    /// Initialises an empty [`MyCell`] to be filled progressively with [`MyCell::shift_position()`].
    pub fn new_empty(target_area: u32, divide_area: u32, cell_type: CellType) -> EmptyCell<Self> {
        EmptyCell::new(Self {
            cell: Cell::new_empty(target_area).into_cell(),
            chem_com: Com { pos: Pos::new(0., 0.), mass: 0 },
            newborn_target_area: target_area,
            divide_area,
            cell_type,
        }).expect("cell was not empty")
    }
    
    /// Returns the total concentration of the chemical perceived by the cell.
    pub fn chem_mass(&self) -> u32 {
        self.chem_com.mass
    }

    /// Returns the center of the cell weighted by the chemical concentration at each cell position.
    pub fn chem_center(&self) -> Pos<FloatType> {
        self.chem_com.pos
    }

    /// Sets the area at which the cell divides when
    /// [`Environment::reproduce()`](crate::my_environment::MyEnvironment::reproduce()) is called.
    pub fn set_divide_area(&mut self, value: u32) {
        self.divide_area = value;
    }

    /// Adds or removes the chemical concentration `chem_at` at position `pos` from the cell.
    pub fn shift_chem<B: Boundary<Coord = FloatType>>(&mut self, pos: Pos<usize>, chem_at: u32, adding: bool, boundary: &B) {
        let shifted = self.chem_com.shift(
            Com { pos: pos.cast_as(), mass: chem_at },
            adding,
            boundary
        );
        match shifted {
            Ok(new_com) => self.chem_com = new_com,
            Err(e) => log::warn!("Failed to shift chem center: {e}")
        }
    }

    /// Updates parameters of the cell (called by [`Pond::step()`](cellulars_lib::traits::step::Step::step())).
    pub fn update(&mut self) {
        if let CellType::Dividing = self.cell_type && self.target_area() < self.divide_area {
            let new_target_area = self.target_area() + 1;
            self.cell.target_area = new_target_area;
        }
    }
}

impl Cellular for MyCell {
    fn target_area(&self) -> u32 {
        self.cell.target_area()
    }

    fn area(&self) -> u32 {
        self.cell.area()
    }

    fn is_empty(&self) -> bool {
        self.cell.is_empty()
    }

    fn shift_position(
        &mut self,
        pos: Pos<usize>,
        adding: bool,
        bound: &impl Boundary<Coord = FloatType>
    ) -> Result<(), ShiftError> {
        self.cell.shift_position(pos, adding, bound)
    }
}

impl HasCenter for MyCell {
    fn center(&self) -> Pos<FloatType> {
        self.cell.center()
    }
}

impl Alive for MyCell {
    fn is_alive(&self) -> bool {
        self.cell.is_alive()
    }

    fn apoptosis(&mut self) {
        self.cell.apoptosis()
    }

    fn birth(&self) -> EmptyCell<Self> {
        let mut basic_cell = self.cell.birth().into_cell();
        basic_cell.target_area = self.newborn_target_area;
        EmptyCell::new(Self {
            chem_com: Com { pos: basic_cell.center(), mass: 0 },
            cell: basic_cell,
            ..self.clone()
        }).expect("failed to create empty cell")
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

    fn make_unsafe_boundary() -> UnsafePeriodicBoundary<FloatType> {
        UnsafePeriodicBoundary::new(Rect::new((0., 0.).into(), (100., 100.).into()))
    }
    
    fn make_test_cell() -> MyCell {
        MyCell::new_empty(
            100,
            200,
            CellType::Migrating,
        ).into_cell()
    }

    #[test]
    fn test_shift_position_area_and_center() {
        let mut cell = make_test_cell();
        let bound = make_unsafe_boundary();

        cell.shift_position(Pos::new(10, 10), true, &bound).unwrap();
        assert_eq!(cell.area(), 1);
        assert_eq!(cell.center(), Pos::new(10.0, 10.0));

        cell.shift_position(Pos::new(20, 20), true, &bound).unwrap();
        assert_eq!(cell.area(), 2);
        assert_eq!(cell.center(), Pos::new(15.0, 15.0));

        cell.shift_position(Pos::new(10, 10), false, &bound).unwrap();
        assert_eq!(cell.area(), 1);
        assert_eq!(cell.center(), Pos::new(20.0, 20.0));
    }

    #[test]
    fn test_shift_position_chem_center_and_mass() {
        let bound = make_unsafe_boundary();
        let mut cell = make_test_cell();

        // Add chem at (2, 3) with value 10
        cell.shift_chem(Pos::new(2, 3), 10, true, &bound);
        assert_eq!(cell.chem_com.mass, 10);
        assert_eq!(cell.chem_com.pos, Pos::new(2., 3.));

        // Add chem at (4, 5) with value 10
        cell.shift_chem(Pos::new(4, 5), 10, true, &bound);
        assert_eq!(cell.chem_com.mass, 20);
        assert_eq!(cell.chem_com.pos, Pos::new(3., 4.));

        // Remove chem from (2, 3)
        cell.shift_chem(Pos::new(2, 3), 10, false, &bound);
        assert_eq!(cell.chem_com.mass, 10);
        assert_eq!(cell.chem_com.pos, Pos::new(4., 5.));
    }
}
