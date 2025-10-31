use crate::evolution::genome::Genome;
use crate::evolution::grn::Grn;
use crate::evolution::selector::Fit;
use cellulars_lib::basic_cell::{Alive, BasicCell, Cellular};
use cellulars_lib::constants::CellIndex;
use cellulars_lib::positional::boundaries::Boundary;
use cellulars_lib::positional::pos::Pos;
use std::ops::{Deref, DerefMut};

#[derive(Clone, Debug)]
pub struct Cell {
    pub basic_cell: BasicCell,
    pub divide_area: u32,
    pub delta_perimeter: Option<i32>,
    pub perimeter: u32,
    pub target_perimeter: u32,
    pub genome: Grn<1, 1>,
    pub ancestor: Option<CellIndex>
}

impl Cell {
    /// Initialises an empty migrating `Cell` to be filled progressively with `shift_position()`.
    pub fn new_empty(target_area: u32, target_perimeter: u32, divide_area: u32, genome: Grn<1, 1>) -> Self {
        Self {
            basic_cell: BasicCell::new_empty(target_area),
            delta_perimeter: None,
            perimeter: 0,
            ancestor: None,
            target_perimeter,
            divide_area,
            genome
        }
    }

    pub fn divide_area(&self) -> u32 {
        self.divide_area
    }

    pub fn set_divide_area(&mut self, value: u32) {
        self.divide_area = value;
    }
    
    pub fn is_migrating(&self) -> bool {
        self.genome.nth_output_gene(0).active
    }
    
    pub fn is_dividing(&self) -> bool {
        !self.is_migrating()
    }

    pub fn update(&mut self) {
        if self.is_dividing() && self.target_area() < self.divide_area {
            let new_target_area = self.target_area() + 1;
            self.set_target_area(new_target_area);
        }
        self.genome.update_expression();
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
            perimeter: 0,
            ..self.clone()
        }
    }
}

impl Fit for Cell {
    fn fitness(&self) -> f32 {
        todo!()
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
            200,
            Grn::empty(),
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
}
