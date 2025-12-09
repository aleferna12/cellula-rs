use crate::evolution::selector::Fit;
use cellulars_lib::basic_cell::{shifted_com, Alive, BasicCell, Cellular, RelCell};
use cellulars_lib::constants::CellIndex;
use cellulars_lib::positional::boundaries::Boundary;
use cellulars_lib::positional::pos::Pos;
use std::ops::{Deref, DerefMut};
use crate::evolution::bit_genome::BitGenome;

#[derive(Clone, Debug)]
pub struct Cell {
    pub basic_cell: BasicCell,
    pub chem_center: Pos<f32>,
    pub chem_mass: u32,
    pub delta_perimeter: Option<i32>,
    pub perimeter: u32,
    pub target_perimeter: u32,
    pub genome: BitGenome,
    pub ancestor: Option<CellIndex>
}

impl Cell {
    /// Initialises an empty migrating `Cell` to be filled progressively with `shift_position()`.
    pub fn new_empty(target_area: u32, target_perimeter: u32, genome: BitGenome) -> Self {
        Self {
            basic_cell: BasicCell::new_empty(target_area),
            chem_center: Pos::new(0., 0.),
            chem_mass: 0,
            delta_perimeter: None,
            perimeter: 0,
            ancestor: None,
            target_perimeter,
            genome
        }
    }

    pub fn chem_center(&self) -> Pos<f32> {
        self.chem_center
    }

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
        self.basic_cell.shift_position(pos, add, bound);
        self.perimeter = self.perimeter
            .checked_add_signed(self.delta_perimeter.expect("`delta_perimeter` not set"))
            .expect("overflow when shifting perimeter");
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
        Self { 
            basic_cell: self.basic_cell.birth(),
            chem_mass: 0,
            perimeter: 0,
            ..self.clone()
        }
    }
}

pub struct FitCell<'c> {
    pub cell: &'c RelCell<Cell>,
    pub half_fit: f32
}

impl Fit for FitCell<'_> {
    fn fitness(&self) -> f32 {
        let ratio = self.half_fit / self.cell.chem_mass as f32;
        1. / (1. + ratio * ratio)
    }
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
            250,
            BitGenome::new(0, 0, 0., 1).unwrap(),
        )
    }

    #[test]
    fn test_shift_position_area_and_center() {
        let mut cell = make_test_cell();
        cell.delta_perimeter = Some(1);
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
