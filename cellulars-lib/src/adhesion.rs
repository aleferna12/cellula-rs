use crate::basic_cell::RelCell;
use crate::lattice_entity::LatticeEntity;
use crate::lattice_entity::LatticeEntity::*;
use std::fmt::Debug;

pub trait AdhesionSystem {
    fn adhesion_energy<C>(&self, entity1: LatticeEntity<&RelCell<C>>, entity2: LatticeEntity<&RelCell<C>>) -> f32;
}

pub struct StaticAdhesion {
    pub cell_energy: f32,
    pub medium_energy: f32,
    pub solid_energy: f32
}

impl AdhesionSystem for StaticAdhesion {
    fn adhesion_energy<C>(
        &self, 
        entity1: LatticeEntity<&RelCell<C>>,
        entity2: LatticeEntity<&RelCell<C>>
    ) -> f32 {
        match (entity1, entity2) {
            (SomeCell(c1), SomeCell(c2)) => {
                if c1.spin == c2.spin {
                    0.
                } else {
                    2. * self.cell_energy
                }
            }
            (SomeCell(_), Medium) | (Medium, SomeCell(_)) => self.medium_energy,
            (SomeCell(_), Solid) | (Solid, SomeCell(_)) => self.solid_energy,
            _ => 0.
        }
    }
}

// #[cfg(test)]
// mod tests {
//     use super::*;
//     use crate::constants::Spin;
//     use crate::genetics::mock_genome::MockGenome;
// 
//     // Helper to create a RelCell with given spin and mom by mocking and overriding
//     fn make_rel_with_spin(spin: Spin, mom: Spin) -> RelCell<Cell<MockGenome>> {
//         RelCell {
//             spin,
//             mom,
//             cell: Cell::new_empty(
//                 10,
//                 200,
//                 MockGenome::new(0)
//             )
//         }
//     }
//     
//     fn make_static_adhesion() -> StaticAdhesion {
//         StaticAdhesion {
//             cell_energy: 3.,
//             medium_energy: 1.5,
//             solid_energy: 2.,
//         }
//     }
// 
//     #[test]
//     fn test_static_adhesion() {
//         let static_adhesion = make_static_adhesion();
// 
//         let cell1 = make_rel_with_spin(1, 1);
//         let cell2 = make_rel_with_spin(2, 1);
// 
//         assert_eq!(
//             static_adhesion.adhesion_energy(SomeCell(&cell1), SomeCell(&cell1)),
//             0.
//         );
//         assert_eq!(
//             static_adhesion.adhesion_energy(SomeCell(&cell1), SomeCell(&cell2)),
//             2. * static_adhesion.cell_energy
//         );
//         assert_eq!(
//             static_adhesion.adhesion_energy(SomeCell(&cell1), Medium),
//             static_adhesion.medium_energy
//         );
//         assert_eq!(
//             static_adhesion.adhesion_energy(Solid, SomeCell(&cell1)),
//             static_adhesion.solid_energy
//         );
//     }
// 
//     #[test]
//     fn test_clonal_adhesion() {
//         let max_spin = 5;
//         let mut clonal_adhesion = ClonalAdhesion::new(max_spin, make_static_adhesion());
// 
//         let cell1 = make_rel_with_spin(1, 1);
//         let cell2 = make_rel_with_spin(2, 1);
//         // Initially clone_pairs empty
//         assert_eq!(
//             clonal_adhesion.adhesion_energy(SomeCell(&cell1), SomeCell(&cell1)),
//             0.
//         );
//         assert_eq!(
//             clonal_adhesion.adhesion_energy(SomeCell(&cell1), SomeCell(&cell2)),
//             2. * clonal_adhesion.static_adhesion.medium_energy
//         );
// 
//         // Manually set clone pair between spin 1 and 2
//         clonal_adhesion.clone_pairs[(1, 2)] = true;
//         assert_eq!(
//             clonal_adhesion.adhesion_energy(SomeCell(&cell1), SomeCell(&cell2)),
//             2. * clonal_adhesion.static_adhesion.cell_energy
//         );
//         // ClonalAdhesion falls back to StaticAdhesion for Medium and Solid
//         assert_eq!(
//             clonal_adhesion.adhesion_energy(SomeCell(&cell1), Medium),
//             clonal_adhesion.static_adhesion.medium_energy
//         );
//         assert_eq!(
//             clonal_adhesion.adhesion_energy(Solid, SomeCell(&cell1)),
//             clonal_adhesion.static_adhesion.solid_energy
//         );
//     }
// }
