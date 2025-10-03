use crate::chem_environment::ChemEnvironment;
use cellulars_lib::adhesion::{AdhesionSystem, StaticAdhesion};
use cellulars_lib::basic_cell::RelCell;
use cellulars_lib::constants::CellIndex;
use cellulars_lib::entity::Entity;
use cellulars_lib::environment::Habitable;
use cellulars_lib::symmetric_table::SymmetricTable;
use std::collections::HashSet;

#[derive(Clone)]
pub struct ClonalAdhesion {
    pub clone_energy: f32,
    pub static_adhesion: StaticAdhesion,
    pub clones_table: SymmetricTable<bool>
}

impl ClonalAdhesion {
    pub fn new(clone_energy: f32, static_adhesion: StaticAdhesion, clones_table: SymmetricTable<bool>) -> Self {
        Self {
            clone_energy,
            static_adhesion,
            clones_table
        }
    }

    pub fn update_clones(
        &mut self,
        cell_index: CellIndex,
        env: &ChemEnvironment
    ) {
        let cell = env.cells().get_cell(cell_index);
        let cell_neighs = env.cell_neighbours(cell, 2.0);

        let mom_cell = env.cells().get_cell(cell.mom);
        let mom_neighs = env.cell_neighbours(
            mom_cell,
            2.0
        );

        let mom_clones = HashSet::from_iter(
            (0..env.cells().n_cells())
                .filter_map(|index| {
                    if self.clones_table[(mom_cell.index as usize, index as usize)] {
                        Some(Entity::Some(index))
                    } else {
                        None
                    }
                })
        );
        for spin in mom_clones.difference(&mom_neighs) {
            if let Entity::Some(index) = spin {
                self.clones_table[(mom_cell.index as usize, *index as usize)] = false;
            }
        }
        let clones: Vec<_> = cell_neighs.intersection(&mom_clones).collect();
        for spin in &clones {
            if let Entity::Some(index) = spin {
                self.clones_table[(cell.index as usize, *index as usize)] = true;
            }
        }
        self.clones_table[(cell.index as usize, mom_cell.index as usize)] = true;
    }
}

impl AdhesionSystem for ClonalAdhesion {
    fn adhesion_energy<C>(&self, entity1: Entity<&RelCell<C>>, entity2: Entity<&RelCell<C>>) -> f32 {
        if let (Entity::Some(c1), Entity::Some(c2)) = (entity1, entity2) {
            if c1.index == c2.index {
                return 0.
            }
            if self.clones_table[(c1.index as usize, c2.index as usize)] {
                return 2. * self.clone_energy;
            }
        }
        // Handle all other cases
        self.static_adhesion.adhesion_energy(entity1, entity2)
    }
}

#[cfg(test)]
mod tests {
    use crate::clonal_adhesion::ClonalAdhesion;
    use cellulars_lib::adhesion::{AdhesionSystem, StaticAdhesion};
    use cellulars_lib::basic_cell::{BasicCell, RelCell};
    use cellulars_lib::constants::CellIndex;
    use cellulars_lib::entity::Entity;
    use cellulars_lib::symmetric_table::SymmetricTable;

    // Helper to create a RelCell with given index and mom by mocking and overriding
    fn make_rel_with_index(index: CellIndex, mom: CellIndex) -> RelCell<BasicCell> {
        RelCell {
            index,
            mom,
            cell: BasicCell::new_empty(10)
        }
    }

    fn make_static_adhesion() -> StaticAdhesion {
        StaticAdhesion {
            cell_energy: 2.,
            medium_energy: 3.,
            solid_energy: 4.,
        }
    }

    #[test]
    fn test_clonal_adhesion() {
        let max_index = 5;
        let mut clonal_adhesion = ClonalAdhesion::new(
            1.,
            make_static_adhesion(),
            SymmetricTable::new(max_index)
        );

        let cell1 = make_rel_with_index(1, 1);
        let cell2 = make_rel_with_index(2, 1);
        // Initially clone_pairs empty
        assert_eq!(
            clonal_adhesion.adhesion_energy(Entity::Some(&cell1), Entity::Some(&cell1)),
            0.
        );
        assert_eq!(
            clonal_adhesion.adhesion_energy(Entity::Some(&cell1), Entity::Some(&cell2)),
            2. * clonal_adhesion.static_adhesion.cell_energy
        );

        // Manually set clone pair between indexes 1 and 2
        clonal_adhesion.clones_table[(1, 2)] = true;
        assert_eq!(
            clonal_adhesion.adhesion_energy(Entity::Some(&cell1), Entity::Some(&cell2)),
            2. * clonal_adhesion.clone_energy
        );
        // ClonalAdhesion falls back to StaticAdhesion for Medium and Solid
        assert_eq!(
            clonal_adhesion.adhesion_energy(Entity::Some(&cell1), Entity::Medium),
            clonal_adhesion.static_adhesion.medium_energy
        );
        assert_eq!(
            clonal_adhesion.adhesion_energy(Entity::Solid, Entity::Some(&cell1)),
            clonal_adhesion.static_adhesion.solid_energy
        );
    }
}