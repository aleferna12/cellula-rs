use crate::cellular::{Cellular, RelCell};
use crate::constants::Spin;
use crate::environment::Environment;
use crate::lattice_entity::LatticeEntity;
use crate::lattice_entity::LatticeEntity::*;
use crate::positional::neighbourhood::Neighbourhood;
use crate::spatial::Spatial;
use crate::symmetric_table::SymmetricTable;
use std::collections::HashSet;
use std::fmt::Debug;

pub trait AdhesionSystem {
    fn adhesion_energy<C>(&self, entity1: LatticeEntity<&RelCell<C>>, entity2: LatticeEntity<&RelCell<C>>) -> f32;
}

pub struct ClonalAdhesion {
    pub static_adhesion: StaticAdhesion,
    // TODO!: try red-black trees in each cell
    pub clone_pairs: SymmetricTable<bool>
}

impl ClonalAdhesion {
    pub fn new(max_spin: Spin, static_adhesion: StaticAdhesion) -> Self {
        Self {
            static_adhesion,
            clone_pairs: SymmetricTable::new(max_spin as usize)
        }
    }
    
    // TODO!: remove and implement neighbour tracking
    pub fn update_clones(
        &mut self,
        cell_spin: Spin,
        cell_search_radius: f32,
        env: &Environment<
            impl Cellular + Debug, 
            impl Neighbourhood,
            impl Spatial
        >
    ) -> Option<Vec<Spin>> {
        let entity = env.cells.get_entity(cell_spin);
        if entity.spin() < LatticeEntity::first_cell_spin() {
            return None;
        }
        
        let cell = entity.unwrap_cell();
        let cell_neighs = env.cell_neighbours(
            cell,
            cell_search_radius,
        );
        
        let mom_cell = env.cells.get_entity(cell.mom).expect_cell("cell's mom is not a cell");
        let mom_neighs = env.cell_neighbours(
            mom_cell,
            cell_search_radius,
        );
        
        let mom_clones = HashSet::from_iter(
            (LatticeEntity::first_cell_spin()..(env.cells.n_cells() + LatticeEntity::first_cell_spin()))
                .filter(|spin| {
                    self.clone_pairs[(mom_cell.spin as usize, *spin as usize)]
                })
        );
        for spin in mom_clones.difference(&mom_neighs) {
            self.clone_pairs[(mom_cell.spin as usize, *spin as usize)] = false;
        }
        let clones: Vec<_> = cell_neighs.intersection(&mom_clones).copied().collect();
        for spin in &clones {
            self.clone_pairs[(cell.spin as usize, *spin as usize)] = true;
        }
        self.clone_pairs[(cell.spin as usize, mom_cell.spin as usize)] = true;
        Some(clones)
    }
}

impl AdhesionSystem for ClonalAdhesion {
    fn adhesion_energy<C>(&self, entity1: LatticeEntity<&RelCell<C>>, entity2: LatticeEntity<&RelCell<C>>) -> f32 {
        if let (SomeCell(c1), SomeCell(c2)) = (entity1, entity2) {
            if c1.spin == c2.spin {
                return 0.
            }
            if self.clone_pairs[(c1.spin as usize, c2.spin as usize)] {
                return 2. * self.static_adhesion.cell_energy;
            }
            return 2. * self.static_adhesion.medium_energy;
        }
        // TODO: Its potentially more efficient to handle all cases here since we skip checking cell-cell again
        // Handle all other cases
        self.static_adhesion.adhesion_energy(entity1, entity2)
    }
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
