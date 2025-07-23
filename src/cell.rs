use crate::constants::Spin;
use crate::environment::LatticeEntity;
use crate::genome::CellType;
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
pub struct RelCell<G> {
    pub spin: Spin,
    pub mom: Spin,
    pub(crate) cell: Cell<G>
}

impl<G> RelCell<G> {
    /// Creates a mock cell with spin and mom = `LatticeEntity<()>::first_cell_spin()` for testing.
    pub fn mock(cell: Cell<G>) -> Self {
        RelCell {
            spin: LatticeEntity::first_cell_spin(),
            mom: LatticeEntity::first_cell_spin(),
            cell
        }
    }
}

impl<G> Deref for RelCell<G> {
    type Target = Cell<G>;

    fn deref(&self) -> &Self::Target {
        &self.cell
    }
}

impl<G> DerefMut for RelCell<G> {
    fn deref_mut(&mut self) -> &mut <Self as Deref>::Target {
        &mut self.cell
    }
}

#[derive(Clone, Debug)]
pub struct Cell<G> {
    pub area: u32,
    pub target_area: u32,
    pub center: Pos<f32>,
    pub light_center: Pos<f32>,
    pub light_mass: u32,
    pub cell_type: CellType,
    pub genome: G
}

impl<G> Cell<G> {
    /// Initialises an empty migrating `Cell` to be filled progressively with `shift_position()`.
    pub fn new(target_area: u32, genome: G) -> Self {
        Self {
            area: 0,
            target_area,
            center: Pos::new(0., 0.),
            light_center: Pos::new(0., 0.),
            light_mass: 0,
            cell_type: CellType::Migrate,
            genome
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
        } else { 
            self.light_center = self.center;
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
    let new_com = Pos::new(
        com.x + dx * added_mass / new_mass,
        com.y + dy * added_mass / new_mass,
    );
    // We call this to rewrap the position if necessary
    bound.valid_pos(new_com).expect("Shifted the center to outside the available space").into()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::genome::MockGenome;
    use crate::positional::boundary::UnsafePeriodicBoundary;
    use crate::positional::pos::Pos;
    use crate::positional::rect::Rect;

    fn make_unsafe_boundary() -> UnsafePeriodicBoundary<f32> {
        UnsafePeriodicBoundary::new(Rect::new((0., 0.).into(), (100., 100.).into()))
    }
    
    fn make_test_cell() -> Cell<MockGenome> {
        Cell::new(100, MockGenome::new(0))
    }

    #[test]
    fn test_shift_position_area_and_center() {
        let mut cell = make_test_cell();
        let bound = make_unsafe_boundary();

        cell.shift_position(Pos::new(10, 10), 0, true, &bound);
        assert_eq!(cell.area, 1);
        assert_eq!(cell.center, Pos::new(10.0, 10.0));

        cell.shift_position(Pos::new(20, 20), 0, true, &bound);
        assert_eq!(cell.area, 2);
        assert_eq!(cell.center, Pos::new(15.0, 15.0));

        cell.shift_position(Pos::new(10, 10), 0, false, &bound);
        assert_eq!(cell.area, 1);
        assert_eq!(cell.center, Pos::new(20.0, 20.0));
    }

    #[test]
    fn test_shift_position_light_center_and_mass() {
        let bound = make_unsafe_boundary();
        let mut cell = make_test_cell();

        // Add light at (2, 3) with value 10
        cell.shift_position(Pos::new(2, 3), 10, true, &bound);
        assert_eq!(cell.light_mass, 10);
        assert_eq!(cell.light_center, Pos::new(2., 3.));

        // Add light at (4, 5) with value 10
        cell.shift_position(Pos::new(4, 5), 10, true, &bound);
        assert_eq!(cell.light_mass, 20);
        assert_eq!(cell.light_center, Pos::new(3., 4.));

        // Remove light from (2, 3)
        cell.shift_position(Pos::new(2, 3), 10, false, &bound);
        assert_eq!(cell.light_mass, 10);
        assert_eq!(cell.light_center, Pos::new(4., 5.));
    }
}
