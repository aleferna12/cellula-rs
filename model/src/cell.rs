use crate::genetics::genome::Genome;
use crate::genetics::grn::Grn;
use cellulars_lib::cellular::Cellular;
use cellulars_lib::positional::boundary::Boundary;
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::selector::Fit;

pub trait CanMigrate: Cellular {
    fn is_migrating(&self) -> bool;
}

#[derive(Clone, Debug)]
pub struct Cell {
    pub area: u32,
    pub target_area: u32,
    pub divide_area: u32,
    pub center: Pos<f32>,
    pub chem_center: Pos<f32>,
    pub chem_mass: f32,
    pub genome: Grn<1, 1>
}

impl Cell {
    /// Initialises an empty migrating `Cell` to be filled progressively with `shift_position()`.
    pub fn new_empty(target_area: u32, divide_area: u32, genome: Grn<1, 1>) -> Self {
        Self {
            area: 0,
            target_area,
            divide_area,
            center: Pos::new(0., 0.),
            chem_center: Pos::new(0., 0.),
            chem_mass: 0.,
            genome
        }
    }

    pub fn chem_center(&self) -> Pos<f32> {
        self.chem_center
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

    // TODO!: Pretty sure i can merge this and shift_position
    pub(crate) fn shift_chem<B: Boundary<Coord=f32>>(&mut self, pos: Pos<usize>, chem_at: f32, add: bool, bound: &B) {
        let shift = if add { 1 } else { -1 };
        if let Some(new_chem) = shifted_com(
            self.chem_center,
            pos,
            self.chem_mass,
            chem_at,
            shift,
            bound
        ) {
            self.chem_center = new_chem;
        } else {
            self.chem_center = self.center;
        }
        self.chem_mass += shift as f32 * chem_at;
        self.genome.input_signals[0] = self.chem_mass;
    }
}

impl Cellular for Cell {
    fn target_area(&self) -> u32 {
        self.target_area
    }

    fn set_target_area(&mut self, value: u32) {
        self.target_area = value;
    }

    fn area(&self) -> u32 {
        self.area
    }

    fn center(&self) -> Pos<f32> {
        self.center
    }

    fn shift_position<B: Boundary<Coord = f32>>(
        &mut self,
        pos: Pos<usize>,
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
        }
        self.area = self.area.saturating_add_signed(shift);
    }

    fn update(&mut self) {
        if self.is_dividing() && self.target_area < self.divide_area {
            self.target_area += 1;
        }
        self.genome.update_expression();
    }

    fn birth(&self) -> Self {
        Self::new_empty(
            self.target_area,
            self.divide_area,
            self.genome.clone()
        )
    }

    fn die(&mut self) {
        self.target_area = 0;
    }

    fn is_alive(&self) -> bool {
        self.target_area > 0
    }

    fn is_valid(&self) -> bool {
        self.area > 0
    }
}

impl Fit for Cell {
    fn fitness(&self) -> f32 {
        self.chem_mass
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
    let new_com = Pos::new(
        com.x + dx * added_mass / new_mass,
        com.y + dy * added_mass / new_mass,
    );
    // We call this to rewrap the position if necessary
    bound.valid_pos(new_com).expect("shifted the center to outside the available space").into()
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::genetics::mock_genome::MockGenome;
//     use crate::positional::boundary::UnsafePeriodicBoundary;
//     use crate::positional::pos::Pos;
//     use crate::positional::rect::Rect;
// 
//     fn make_unsafe_boundary() -> UnsafePeriodicBoundary<f32> {
//         UnsafePeriodicBoundary::new(Rect::new((0., 0.).into(), (100., 100.).into()))
//     }
//     
//     fn make_test_cell() -> Cell<MockGenome> {
//         Cell::new_empty(
//             100,
//             200,
//             MockGenome::new(0)
//         )
//     }
// 
//     #[test]
//     fn test_shift_position_area_and_center() {
//         let mut cell = make_test_cell();
//         let bound = make_unsafe_boundary();
// 
//         cell.shift_position(Pos::new(10, 10), true, &bound);
//         assert_eq!(cell.area, 1);
//         assert_eq!(cell.center, Pos::new(10.0, 10.0));
// 
//         cell.shift_position(Pos::new(20, 20), true, &bound);
//         assert_eq!(cell.area, 2);
//         assert_eq!(cell.center, Pos::new(15.0, 15.0));
// 
//         cell.shift_position(Pos::new(10, 10), false, &bound);
//         assert_eq!(cell.area, 1);
//         assert_eq!(cell.center, Pos::new(20.0, 20.0));
//     }
// 
//     #[test]
//     fn test_shift_position_chem_center_and_mass() {
//         let bound = make_unsafe_boundary();
//         let mut cell = make_test_cell();
// 
//         // Add chem at (2, 3) with value 10
//         cell.shift_chem(Pos::new(2, 3), 10., true, &bound);
//         assert_eq!(cell.chem_mass, 10.);
//         assert_eq!(cell.chem_center, Pos::new(2., 3.));
// 
//         // Add chem at (4, 5) with value 10
//         cell.shift_chem(Pos::new(4, 5), 10., true, &bound);
//         assert_eq!(cell.chem_mass, 20.);
//         assert_eq!(cell.chem_center, Pos::new(3., 4.));
// 
//         // Remove chem from (2, 3)
//         cell.shift_chem(Pos::new(2, 3), 10., false, &bound);
//         assert_eq!(cell.chem_mass, 10.);
//         assert_eq!(cell.chem_center, Pos::new(4., 5.));
//     }
// }
