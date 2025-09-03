use std::ops::{Deref, DerefMut};
use crate::cell::Cell;
use crate::constants::{BoundaryType, NeighbourhoodType};
use cellulars_lib::constants::Spin;
use cellulars_lib::environment::Environment;
use cellulars_lib::lattice::Lattice;
use cellulars_lib::positional::boundary::AsLatticeBoundary;
use cellulars_lib::space::Space;
use cellulars_lib::spatial::Spatial;

pub struct ChemSpace {
    space: Space<BoundaryType>,
    pub chem_lattice: Lattice<u32>,
}

impl ChemSpace {
    pub fn new(bound: BoundaryType) -> Result<Self, <BoundaryType as AsLatticeBoundary>::Error> {
        let space = Space::new(bound)?;
        Ok(Self {
            chem_lattice: space.cell_lattice.clone(),
            space,
        })
    }
}

impl Deref for ChemSpace {
    type Target = Space<BoundaryType>;

    fn deref(&self) -> &Self::Target {
        &self.space
    }
}

impl DerefMut for ChemSpace {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.space
    }
}

impl Spatial for ChemSpace {
    type Boundary = BoundaryType;
    fn cell_lattice(&self) -> &Lattice<Spin> {
        &self.cell_lattice
    }

    fn cell_lattice_mut(&mut self) -> &mut Lattice<Spin> {
        &mut self.cell_lattice
    }

    fn boundary(&self) -> &BoundaryType {
        &self.bound
    }

    fn lattice_boundary(&self) -> &<BoundaryType as AsLatticeBoundary>::LatticeBoundary {
        &self.lat_bound
    }
}

pub type ChemEnvironment = Environment<Cell, NeighbourhoodType, ChemSpace>;

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::cell::{Cell, RelCell};
//     use crate::constants::BoundaryType;
//     use crate::genetics::mock_genome::MockGenome;
//     use crate::positional::neighbourhood::MooreNeighbourhood;
//     use crate::positional::pos::Pos;
//
//     fn space_cell_pair(positions: &[Pos<usize>]) -> (Space<BoundaryType>, RelCell<impl Cellular>) {
//         let mut cell = RelCell::mock(Cell::new_empty(
//             10,
//             20,
//             MockGenome::new(0)
//         ));
//         let mut space = Space::new(BoundaryType::new(Rect::new(
//             (0., 0.).into(),
//             (10., 10.).into()
//         ))).unwrap();
//         for pos in positions {
//             space.cell_lattice[*pos] = cell.spin;
//             cell.shift_position(*pos, true, &space.bound)
//         }
//         (space, cell)
//     }
//
//     #[test]
//     fn test_box_cell_positions() {
//         let positions = [
//             Pos::new(5, 5),
//             Pos::new(5, 6),
//             Pos::new(6, 5),
//             Pos::new(6, 6),
//         ];
//         let (space, cell) = space_cell_pair(&positions);
//         let boxed_positions = space.search_cell_box(&cell, 2.0);
//         assert_eq!(boxed_positions.len(), positions.len());
//         for pos in &positions {
//             assert!(boxed_positions.contains(pos));
//         }
//     }
//
//     #[test]
//     fn test_contiguous_cell_positions() {
//         let positions = [
//             Pos::new(5, 5),
//             Pos::new(5, 6),
//             Pos::new(6, 5),
//             Pos::new(6, 6),
//         ];
//         let (space, cell) = space_cell_pair(&positions);
//         let neighbourhood = MooreNeighbourhood::new(1);
//         let contiguous_positions = space.search_cell_contiguous(&cell, &neighbourhood);
//
//         // Should find all 4 contiguous positions
//         assert_eq!(contiguous_positions.len(), positions.len());
//         for pos in &positions {
//             assert!(contiguous_positions.contains(pos));
//         }
//     }
//
//     #[test]
//     fn test_contiguous_cell_positions_discontiguous() {
//         let positions = [
//             Pos::new(5, 5),
//             Pos::new(5, 6),
//             Pos::new(7, 7), // discontiguous point
//         ];
//         let (space, cell) = space_cell_pair(&positions);
//         let neighbourhood = MooreNeighbourhood::new(1);
//         let result = space.search_cell_contiguous(&cell, &neighbourhood);
//
//         assert_eq!(result.len(), 2);
//     }
//
//     #[test]
//     fn test_cell_neighbours() {
//         let positions = [
//             Pos::new(2, 2),
//             Pos::new(2, 1),
//         ];
//         let (mut space, cell) = space_cell_pair(&positions);
//         let neighbourhood = MooreNeighbourhood::new(1);
//
//         let neighbour_spins = [cell.spin + 1, cell.spin + 2];
//         space.cell_lattice[Pos::new(1, 2)] = neighbour_spins[0];
//         space.cell_lattice[Pos::new(2, 0)] = neighbour_spins[1];
//
//         let neighs = space.cell_neighbours(&cell, 1.0, &neighbourhood);
//
//         assert!(neighs.contains(&neighbour_spins[0]));
//         assert!(neighs.contains(&neighbour_spins[1]));
//         assert!(!neighs.contains(&cell.spin));
//     }
// }