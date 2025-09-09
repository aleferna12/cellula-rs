use crate::chem_environment::ChemEnvironment;
use cellulars_lib::adhesion::{AdhesionSystem, StaticAdhesion};
use cellulars_lib::basic_cell::RelCell;
use cellulars_lib::constants::Spin;
use cellulars_lib::environment::Habitable;
use cellulars_lib::lattice_entity::LatticeEntity;
use cellulars_lib::lattice_entity::LatticeEntity::SomeCell;
use cellulars_lib::symmetric_table::SymmetricTable;
use std::collections::HashSet;

#[derive(Clone)]
pub struct ClonalAdhesion {
    pub clone_energy: f32,
    pub static_adhesion: StaticAdhesion,
    pub clone_pairs: SymmetricTable<bool>
}

impl ClonalAdhesion {
    pub fn new(max_spin: Spin, clone_energy: f32, static_adhesion: StaticAdhesion) -> Self {
        Self {
            clone_energy,
            static_adhesion,
            clone_pairs: SymmetricTable::new(max_spin as usize)
        }
    }

    // TODO!: remove and implement neighbour tracking
    pub fn update_clones(
        &mut self,
        cell_spin: Spin,
        env: &ChemEnvironment
    ) -> Option<Vec<Spin>> {
        let entity = env.cells().get_entity(cell_spin);
        let cell = entity.unwrap_cell();
        let cell_neighs = env.cell_neighbours(cell, 2.0);

        let mom_cell = env.cells().get_entity(cell.mom).expect_cell("cell's mom is not a cell");
        let mom_neighs = env.cell_neighbours(
            mom_cell,
            2.0
        );

        let mom_clones = HashSet::from_iter(
            (LatticeEntity::first_cell_spin()..(env.cells().n_cells() + LatticeEntity::first_cell_spin()))
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
    use cellulars_lib::constants::Spin;
    use cellulars_lib::lattice_entity::LatticeEntity::{Medium, Solid, SomeCell};

    // Helper to create a RelCell with given spin and mom by mocking and overriding
    fn make_rel_with_spin(spin: Spin, mom: Spin) -> RelCell<BasicCell> {
        RelCell {
            spin,
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
        let max_spin = 5;
        let mut clonal_adhesion = ClonalAdhesion::new(max_spin, 1., make_static_adhesion());

        let cell1 = make_rel_with_spin(1, 1);
        let cell2 = make_rel_with_spin(2, 1);
        // Initially clone_pairs empty
        assert_eq!(
            clonal_adhesion.adhesion_energy(SomeCell(&cell1), SomeCell(&cell1)),
            0.
        );
        assert_eq!(
            clonal_adhesion.adhesion_energy(SomeCell(&cell1), SomeCell(&cell2)),
            2. * clonal_adhesion.static_adhesion.cell_energy
        );

        // Manually set clone pair between spin 1 and 2
        clonal_adhesion.clone_pairs[(1, 2)] = true;
        assert_eq!(
            clonal_adhesion.adhesion_energy(SomeCell(&cell1), SomeCell(&cell2)),
            2. * clonal_adhesion.clone_energy
        );
        // ClonalAdhesion falls back to StaticAdhesion for Medium and Solid
        assert_eq!(
            clonal_adhesion.adhesion_energy(SomeCell(&cell1), Medium),
            clonal_adhesion.static_adhesion.medium_energy
        );
        assert_eq!(
            clonal_adhesion.adhesion_energy(Solid, SomeCell(&cell1)),
            clonal_adhesion.static_adhesion.solid_energy
        );
    }
}