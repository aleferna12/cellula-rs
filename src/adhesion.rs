use std::collections::HashSet;
use std::ptr;
use crate::cell::RelCell;
use crate::constants::Spin;
use crate::environment::LatticeEntity;
use crate::environment::LatticeEntity::*;
use crate::parameters::StaticAdhesionParameters;

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
    // TODO!: should this be stored as an array in each cell (replacing spin as a cell property)? Benchmark
    //  the current implementation costs almost 25% of performance compared to StaticAdhesion
    //  best solution is probably a big table in the heap that we can access with spins
    pub clone_pairs: HashSet<(Spin, Spin)>
}

impl ClonalAdhesion {
    fn canonicalize(pair: (Spin, Spin)) -> (Spin, Spin) {
        if pair.0 > pair.1 {
            return (pair.1, pair.0);
        }
        pair
    }
    
    // TODO!: This is horrible and doesnt work.
    //  we need to check that mom was a clone with neighbour before inserting the clone
    //  otherwise cells can attach to new groups by being neighbours
    pub fn update_clones(&mut self, spin: Spin, neigh_spins: impl Iterator<Item = Spin>) {
        for pair in self.clone_pairs.iter().copied().collect::<Vec<_>>() {
            if pair.0 == spin || pair.1 == spin {
                self.clone_pairs.remove(&pair);
            }
        }
        for neigh_spin in neigh_spins {
            if neigh_spin < LatticeEntity::first_cell_spin() {
                continue
            }
            self.clone_pairs.insert(Self::canonicalize((spin, neigh_spin)));
        }
    }
}

impl AdhesionSystem for ClonalAdhesion {
    fn adhesion_energy(&self, entity1: LatticeEntity<&RelCell>, entity2: LatticeEntity<&RelCell>) -> f32 {
        if let (SomeCell(c1), SomeCell(c2)) = (entity1, entity2) {
            let canonical = Self::canonicalize((c1.spin, c2.spin));
            if self.clone_pairs.contains(&canonical) {
                return self.static_adhesion.cell_energy;
            }
        }
        // Handle all other cases
        self.static_adhesion.adhesion_energy(entity1, entity2)
    }
}

impl From<StaticAdhesionParameters> for ClonalAdhesion {
    fn from(params: StaticAdhesionParameters) -> Self {
        Self {
            static_adhesion: params.into(),
            clone_pairs: HashSet::default()
        }
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
                if ptr::eq(c1, c2) {
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