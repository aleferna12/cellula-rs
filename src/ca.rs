use std::f32::consts::E;
use rand::{Rng, RngCore};
use crate::cell::Cell;
use crate::environment::Environment;
use crate::lattice::LatticeEntity;
use crate::lattice::LatticeEntity::SomeCell;

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
            let (p1, p2) = if rng.random::<f32>() < 0.5 {
                (edge.p1, edge.p2)
            } else {
                (edge.p2, edge.p1)
            };
            let sigma_from = env.cell_lattice[p1];
            let sigma_to = env.cell_lattice[p2];
            let delta_h = self.delta_hamiltonian(
                env.get_cell(sigma_from),
                env.get_cell(sigma_to)
            );
            if self.accept_copy(rng, delta_h) {
                env.cell_lattice[p2] = sigma_from;
                if let SomeCell(cell) = env.get_cell_mut(sigma_from) {
                    cell.area += 1;
                }
                if let SomeCell(cell) = env.get_cell_mut(sigma_to) {
                    cell.area -= 1;
                }
                let (removed, added) = env.update_edges(p2);
                // TODO: ensure this makes sense for neigh_r > 1
                to_visit += (added as f32 - removed as f32) / edge_per_pos;
            }
            to_visit -= 1f32;
        }
    }

    pub fn accept_copy(&self, rng: &mut impl Rng, delta_h: f32) -> bool {
        delta_h < 0f32 || rng.random::<f32>() < E.powf(-delta_h / self.boltz_t)
    }

    // TODO: add adhesion
    pub fn delta_hamiltonian(&self, cell_from: LatticeEntity<&Cell>, cell_to: LatticeEntity<&Cell>) -> f32 {
        let mut delta_h = 0f32;
        if let SomeCell(cell) = cell_from {
            delta_h += self.delta_hamiltonian_size(
                1,
                cell.area,
                cell.target_area,
                self.size_lambda
            )
        }
        if let SomeCell(cell) = cell_to {
            delta_h += self.delta_hamiltonian_size(
                -1,
                cell.area,
                cell.target_area,
                self.size_lambda
            )
        }
        delta_h
    }

    pub fn delta_hamiltonian_size(&self, delta_area: i32, area: u32, target_area: u32, size_lambda: f32) -> f32 {
        2f32 * size_lambda * delta_area as f32 * (area as f32 - target_area as f32) + size_lambda
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_delta_hamiltonian() {
        let ca = CA::new(12f32, 1f32);
        let cell1 = Cell::new(100, 100);
        let cell2 = Cell::new(100, 100);
        let dh = ca.delta_hamiltonian(SomeCell(&cell1), SomeCell(&cell2));
        assert_eq!(dh, 2f32);
    }
}