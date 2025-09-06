use crate::chem_environment::ChemEnvironment;
use cellulars_lib::adhesion::{AdhesionSystem, StaticAdhesion};
use cellulars_lib::basic_cell::RelCell;
use cellulars_lib::constants::Spin;
use cellulars_lib::environment::Habitable;
use cellulars_lib::lattice_entity::LatticeEntity;
use cellulars_lib::lattice_entity::LatticeEntity::SomeCell;
use cellulars_lib::symmetric_table::SymmetricTable;
use std::collections::HashSet;

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
        env: &ChemEnvironment
    ) -> Option<Vec<Spin>> {
        let entity = env.cells().get_entity(cell_spin);
        let cell = entity.unwrap_cell();
        let cell_neighs = env.cell_neighbours(cell);

        let mom_cell = env.cells().get_entity(cell.mom).expect_cell("cell's mom is not a cell");
        let mom_neighs = env.cell_neighbours(
            mom_cell,
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
                return 2. * self.static_adhesion.cell_energy;
            }
            return 2. * self.static_adhesion.medium_energy;
        }
        // TODO: Its potentially more efficient to handle all cases here since we skip checking cell-cell again
        // Handle all other cases
        self.static_adhesion.adhesion_energy(entity1, entity2)
    }
}