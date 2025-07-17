use crate::cell::RelCell;
use crate::constants::Spin;
use crate::environment::LatticeEntity::*;
use crate::environment::{Environment, LatticeEntity};
use crate::io::parameters::StaticAdhesionParameters;
use crate::spin_table::SpinTable;
use std::collections::HashSet;

pub trait AdhesionSystem {
    // This arguably should receive info about which specific site is being copied 
    // It would be useful for hybrid models where adhesion properties depend on concentration of proteins etc
    // See https://compucell3dreferencemanual.readthedocs.io/en/latest/adhesion_flex_plugin.html
    // Although, in a puritan interpretation of CPM, the Hamiltonian is a property of the system and anything that is 
    // copy-attempt-dependent should be a bias...
    fn adhesion_energy(&self, entity1: LatticeEntity<&RelCell>, entity2: LatticeEntity<&RelCell>) -> f32;
}

// TODO!: Start by modeling adhesion based on whether two cells shared a boundary when they were born
//  To better maintain cluster shape, this can be extended to both determine normal adhesion and also have springs
//  connecting the two clonal cells
pub struct ClonalAdhesion {
    pub static_adhesion: StaticAdhesion,
    // TODO!: try red-black trees in each cell
    pub clone_pairs: SpinTable<bool>
}

impl ClonalAdhesion {
    pub fn new(params: StaticAdhesionParameters, max_spin: Spin) -> Self {
        Self {
            static_adhesion: StaticAdhesion::from(params),
            clone_pairs: SpinTable::new(max_spin)
        }
    }
    
    pub fn update_clones(
        &mut self,
        cell_spin: Spin,
        env: &Environment
    ) -> Option<Vec<Spin>> {
        let entity = env.cells.get_entity(cell_spin);
        if entity.spin() < LatticeEntity::first_cell_spin() {
            return None;
        }
        
        let cell = entity.unwrap_cell();
        let cell_neighs = env.space.cell_neighbours(
            cell,
            env.cell_search_radius,
            &env.neighbourhood
        );
        
        let mom_cell = env.cells.get_entity(cell.mom).expect_cell("cell's mom is not a cell");
        let mom_neighs = env.space.cell_neighbours(
            mom_cell,
            env.cell_search_radius,
            &env.neighbourhood
        );
        
        let mom_clones = HashSet::from_iter(
            (LatticeEntity::first_cell_spin()..=env.cells.n_cells())
                .filter(|spin| {
                    self.clone_pairs[(mom_cell.spin, *spin)]
                })
        );
        for spin in mom_clones.difference(&mom_neighs) {
            self.clone_pairs[(mom_cell.spin, *spin)] = false;
        }
        let clones: Vec<_> = cell_neighs.intersection(&mom_clones).copied().collect();
        for spin in &clones {
            self.clone_pairs[(cell.spin, *spin)] = true;
        }
        self.clone_pairs[(cell.spin, mom_cell.spin)] = true;
        Some(clones)
    }
}

impl AdhesionSystem for ClonalAdhesion {
    fn adhesion_energy(&self, entity1: LatticeEntity<&RelCell>, entity2: LatticeEntity<&RelCell>) -> f32 {
        if let (SomeCell(c1), SomeCell(c2)) = (entity1, entity2) {
            if c1.spin == c2.spin {
                return 0.
            }
            if self.clone_pairs[(c1.spin, c2.spin)] {
                return 2. * self.static_adhesion.cell_energy;
            }
            return 2. * self.static_adhesion.medium_energy;
        }
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
    fn adhesion_energy(&self, entity1: LatticeEntity<&RelCell>, entity2: LatticeEntity<&RelCell>) -> f32 {
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

impl From<StaticAdhesionParameters> for StaticAdhesion {
    fn from(params: StaticAdhesionParameters) -> Self {
        Self {
            cell_energy: params.cell_energy,
            medium_energy: params.medium_energy,
            solid_energy: params.solid_energy,
        }
    }
}