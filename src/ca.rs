use std::f32::consts::E;
use std::ptr;
use rand::Rng;
use crate::boundary::Boundary;
use crate::cell::Cell;
use crate::environment::Environment;
use crate::environment::LatticeEntity;
use crate::environment::LatticeEntity::{Medium, SomeCell, Solid};
use crate::neighbourhood::Neighbourhood;
use crate::pos::Pos2D;

// This could be a module but it's convenient to be able to access the relevant parameters 
// Also we might eventually want to implement multiple CA choices, in which case I can "easily" make CA a trait 
// that just implements `step()`
pub struct CA {
    pub boltz_t: f32,
    pub size_lambda: f32,
    pub cell_energy: f32,
    pub med_energy: f32,
    pub solid_energy: f32
}
impl CA {
    pub fn new(boltz_t: f32, size_lambda: f32, cell_energy: f32, med_energy: f32, solid_energy: f32) -> CA {
        CA {
            boltz_t,
            size_lambda,
            cell_energy,
            med_energy,
            solid_energy
        }
    }

    pub fn step(&self, env: &mut Environment, rng: &mut impl Rng) {
        let mut to_visit = env.edge_book.len() as f32 / env.edge_per_pos() as f32;
        while 0. < to_visit {
            let edge_i = env.edge_book.random_index(rng);
            let edge = env.edge_book.at(edge_i);
            // This is WAY faster than keeping the symmetric edge in EdgeBook (like 2x as fast!)
            // or at least, this is the case when using IndexSet, I would assume its somewhat implementation-dependent
            let (pos_from, pos_to) = if rng.random::<f32>() < 0.5 {
                (edge.p1, edge.p2)
            } else {
                (edge.p2, edge.p1)
            };
            to_visit += self.attempt_site_copy(env, rng, pos_from, pos_to);
            to_visit -= 1.;
        }
    }

    /// Attempts to execute the selected site copy.
    /// 
    /// # Returns:
    /// 
    /// The number of extra updates that the copy attempt incurred (not whether it was successful or not!).
    pub fn attempt_site_copy(
        &self,
        env: &mut Environment,
        rng: &mut impl Rng,
        pos_from: Pos2D<usize>,
        pos_to: Pos2D<usize>
    ) -> f32 {
        let sigma_to = env.cell_lattice[pos_to];
        if sigma_to == Solid.discriminant() {
            return 0.;
        }
        // If was going to copy from a Solid, create a Medium cell instead 
        let sigma_from = {
            let sigma = env.cell_lattice[pos_from];
            if sigma == Solid.discriminant() { Medium.discriminant() } else { sigma }
        };

        let entity_from = env.get_entity(sigma_from);
        let entity_to = env.get_entity(sigma_to);
        let neigh_entities = env.cell_lattice.bound.valid_positions(
            env.neighbourhood.neighbours(pos_to.into())
        ).map(|neigh| {
            env.get_entity(env.cell_lattice[Pos2D::<usize>::from(neigh)])
        });
        
        let delta_h = self.delta_hamiltonian(entity_from, entity_to, neigh_entities);
        if !self.accept_site_copy(rng, delta_h) {
            return 0.;
        }
        
        // Executes the copy
        env.cell_lattice[pos_to] = sigma_from;
        if let SomeCell(cell) = env.get_entity_mut(sigma_from) {
            cell.area += 1;
        }
        if let SomeCell(cell) = env.get_entity_mut(sigma_to) {
            cell.area -= 1;
        }
        let (removed, added) = env.update_edges(pos_to);
        (added as f32 - removed as f32) / env.edge_per_pos() as f32
    }

    pub fn accept_site_copy(&self, rng: &mut impl Rng, delta_h: f32) -> bool {
        delta_h < 0. || rng.random::<f32>() < E.powf(-delta_h / self.boltz_t)
    }

    pub fn delta_hamiltonian<'a>(
        &self,
        entity_from: LatticeEntity<&Cell>,
        entity_to: LatticeEntity<&Cell>,
        neigh_entities: impl Iterator<Item = LatticeEntity<&'a Cell>>
    ) -> f32 {
        let mut delta_h = 0.;
        delta_h += self.delta_hamiltonian_size(entity_from, entity_to);
        delta_h += self.delta_hamiltonian_adhesion(entity_from, entity_to, neigh_entities);
        delta_h
    }
    
    pub fn delta_hamiltonian_size(&self, entity_from: LatticeEntity<&Cell>, entity_to: LatticeEntity<&Cell>) -> f32 {
        let mut delta_h = 0.;
        if let SomeCell(cell) = entity_from {
            delta_h += self.size_energy_diff(true, cell.area, cell.target_area);
        }
        if let SomeCell(cell) = entity_to {
            delta_h += self.size_energy_diff(false, cell.area, cell.target_area);
        }
        delta_h
    }

    // TODO!: test
    pub fn delta_hamiltonian_adhesion<'a>(
        &self,
        entity_from: LatticeEntity<&Cell>,
        entity_to: LatticeEntity<&Cell>,
        neigh_entities: impl Iterator<Item = LatticeEntity<&'a Cell>>
    ) -> f32 {
        let mut energy = 0.;
        for neigh in neigh_entities {
            energy -= self.adhesion_energy(entity_to, neigh);
            energy += self.adhesion_energy(entity_from, neigh);
        }
        energy
    }

    pub fn size_energy_diff(&self, area_increased: bool, area: u32, target_area: u32) -> f32 {
        let delta_area = if area_increased { 1. } else { -1. };
        2. * self.size_lambda * delta_area * (area as f32 - target_area as f32) + self.size_lambda
    }

    pub fn adhesion_energy(&self, entity1: LatticeEntity<&Cell>, entity2: LatticeEntity<&Cell>) -> f32 {
        match (entity1, entity2) {
            (SomeCell(c1), SomeCell(c2)) => {
                if ptr::eq(c1, c2) {
                    0.
                } else {
                    2. * self.cell_energy
                }
            }
            (SomeCell(_), Medium) | (Medium, SomeCell(_)) => self.med_energy,
            (SomeCell(_), Solid) | (Solid, SomeCell(_)) => self.solid_energy,
            _ => 0.
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_delta_hamiltonian_size() {
        let ca = CA::new(12., 1., 10., 20., 20.);
        let cell1 = Cell::new(100, 100);
        let cell2 = Cell::new(100, 100);
        let dh = ca.delta_hamiltonian_size(SomeCell(&cell1), SomeCell(&cell2));
        assert_eq!(dh, 2.);
    }
}