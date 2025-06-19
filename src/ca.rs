use std::f32::consts::E;
use std::ptr;
use rand::{Rng, RngCore};
use crate::cell::Cell;
use crate::environment::Environment;
use crate::lattice::LatticeEntity;
use crate::lattice::LatticeEntity::SomeCell;
use crate::pos::Pos2D;

// This could be a module but it's convenient to be able to access the relevant parameters 
// Also we might eventually want to implement multiple CA choices, in which case I can "easily" make CA a trait
pub struct CA {
    pub boltz_t: f32,
    pub size_lambda: f32,
}
impl CA {
    pub fn new(boltz_t: f32, size_lambda: f32) -> CA {
        CA {
            boltz_t,
            size_lambda
        }
    }

    pub fn step(&self, env: &mut Environment, rng: &mut impl RngCore) {
        // TODO: ensure this makes sense for neigh_r > 1
        let edge_per_pos = env.neigh_r as f32 / 2f32;
        let mut to_visit = env.edge_bk.len() as f32 / edge_per_pos;
        while 0f32 < to_visit {
            let edge_i = env.edge_bk.random_index(rng);
            let edge = env.edge_bk.at(edge_i);
            // TODO: is this really faster than just keeping both edges in the IndexSet? Benchmark
            let (pos_from, pos_to) = if rng.random::<f32>() < 0.5 {
                (edge.p1, edge.p2)
            } else {
                (edge.p2, edge.p1)
            };
            let delta_h = self.delta_hamiltonian(env, pos_from, pos_to);
            if self.accept_copy(rng, delta_h) {
                let sigma_from = env.cell_lattice[pos_from];
                let sigma_to = env.cell_lattice[pos_to];
                env.cell_lattice[pos_to] = sigma_from;
                if let SomeCell(cell) = env.get_entity_mut(sigma_from) {
                    cell.area += 1;
                }
                if let SomeCell(cell) = env.get_entity_mut(sigma_to) {
                    cell.area -= 1;
                }
                let (removed, added) = env.update_edges(pos_to);
                // TODO: ensure this makes sense for neigh_r > 1
                to_visit += (added as f32 - removed as f32) / edge_per_pos;
            }
            to_visit -= 1f32;
        }
    }

    pub fn accept_copy(&self, rng: &mut impl Rng, delta_h: f32) -> bool {
        delta_h < 0f32 || rng.random::<f32>() < E.powf(-delta_h / self.boltz_t)
    }
    
    // TODO: should this just take the entities? think about what information deltaH should have access to
    pub fn delta_hamiltonian(&self, env: &Environment, pos_from: Pos2D<usize>, pos_to: Pos2D<usize>) -> f32 {
        let mut delta_h = 0f32;
        let entity_from = env.get_entity(env.cell_lattice[pos_from]);
        let entity_to = env.get_entity(env.cell_lattice[pos_to]);
        delta_h += self.delta_hamiltonian_size(entity_from, entity_to);
        let neighs = env.cell_lattice
            .validate(pos_to.moore_neighs(env.neigh_r))
            .map(|neigh_pos| env.get_entity(env.cell_lattice[neigh_pos]));
        delta_h += self.delta_hamiltonian_adhesion(entity_from, entity_to, neighs);
        delta_h
    }
    
    pub fn delta_hamiltonian_size(&self, entity_from: LatticeEntity<&Cell>, entity_to: LatticeEntity<&Cell>) -> f32 {
        let mut delta_h = 0f32;
        if let SomeCell(cell) = entity_from {
            delta_h += self.size_energy_diff(true, cell.area, cell.target_area);
        }
        if let SomeCell(cell) = entity_to {
            delta_h += self.size_energy_diff(false, cell.area, cell.target_area);
        }
        delta_h
    }

    pub fn delta_hamiltonian_adhesion<'a>(
        &self,
        entity_from: LatticeEntity<&Cell>,
        entity_to: LatticeEntity<&Cell>,
        neigh_cells: impl Iterator<Item = LatticeEntity<&'a Cell>>
    ) -> f32 {
        let mut energy = 0f32;
        for neigh in neigh_cells {
            energy -= self.adhesion_energy(entity_to, neigh);
            energy += self.adhesion_energy(entity_from, neigh);
        }
        energy
    }

    pub fn size_energy_diff(&self, area_increased: bool, area: u32, target_area: u32) -> f32 {
        let delta_area = if area_increased { 1f32 } else { -1f32 };
        2f32 * self.size_lambda * delta_area * (area as f32 - target_area as f32) + self.size_lambda
    }
    
    pub fn adhesion_energy(&self, entity1: LatticeEntity<&Cell>, entity2: LatticeEntity<&Cell>) -> f32 {
        if let (SomeCell(cell1), SomeCell(cell2)) = (entity1, entity2) {
            if ptr::eq(cell1, cell2) { 
                return 0f32;
            }
            return 10f32;
        }
        20f32
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_delta_hamiltonian_size() {
        let ca = CA::new(12f32, 1f32);
        let cell1 = Cell::new(100, 100);
        let cell2 = Cell::new(100, 100);
        let dh = ca.delta_hamiltonian_size(SomeCell(&cell1), SomeCell(&cell2));
        assert_eq!(dh, 2f32);
    }
}