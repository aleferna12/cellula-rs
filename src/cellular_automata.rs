use crate::genome::CellType;
use crate::adhesion::AdhesionSystem;
use crate::cell::{Cell, RelCell};
use crate::environment::Environment;
use crate::environment::LatticeEntity;
use crate::environment::LatticeEntity::{Medium, Solid, SomeCell};
use crate::positional::boundary::Boundary;
use crate::positional::neighbourhood::Neighbourhood;
use crate::positional::pos::Pos;
use rand::Rng;
use std::f32::consts::E;

// This could be a module but it's convenient to be able to access the relevant parameters 
// Also we might eventually want to implement multiple CA choices, in which case I can "easily" make CA a trait 
// that just implements `step()`
pub struct CellularAutomata<A> {
    pub boltz_t: f32,
    pub size_lambda: f32,
    pub chemotaxis_mu: f32,
    pub adhesion: A
}

impl<A: AdhesionSystem> CellularAutomata<A> {
    pub fn new(boltz_t: f32, size_lambda: f32, chemotaxis_mu: f32, adhesion: A) -> Self {
        Self {
            boltz_t,
            size_lambda,
            chemotaxis_mu,
            adhesion
        }
    }
    
    pub fn step(&self, env: &mut Environment, rng: &mut impl Rng) {
        let mut to_visit = env.edge_book.len() as f32 / env.neighbourhood.n_neighs() as f32;
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
        pos_from: Pos<usize>,
        pos_to: Pos<usize>
    ) -> f32 {
        let spin_to = env.space.cell_lattice[pos_to];
        if spin_to == Solid.discriminant() {
            return 0.;
        }
        // If was going to copy from a Solid, create a Medium cell instead 
        let spin_from = {
            let spin = env.space.cell_lattice[pos_from];
            if spin == Solid.discriminant() { Medium.discriminant() } else { spin }
        };

        let entity_from = env.cells.get_entity(spin_from);
        let entity_to = env.cells.get_entity(spin_to);
        let neigh_entities = env.space.lat_bound.valid_positions(
            env.neighbourhood.neighbours(pos_to.into())
        ).map(|neigh| {
            env.cells.get_entity(env.space.cell_lattice[Pos::<usize>::from(neigh)])
        });
        
        let mut delta_h = self.delta_hamiltonian(entity_from, entity_to, neigh_entities);
        if let SomeCell(cell) = entity_from {
            if env.cells.migrate && let CellType::Migrate = cell.cell_type {
                delta_h += self.chemotaxis_bias(&cell, pos_to, self.chemotaxis_mu, &env.space.bound);
            }
        }
        if !self.accept_site_copy(rng, delta_h) {
            return 0.;
        }
        
        // Executes the copy
        env.space.cell_lattice[pos_to] = spin_from;
        if let SomeCell(cell) = env.cells.get_entity_mut(spin_from) {
            cell.shift_position(pos_to,env.space.light_lattice[pos_to], true, &env.space.bound);
        }
        if let SomeCell(cell) = env.cells.get_entity_mut(spin_to) {
            cell.shift_position(pos_to,env.space.light_lattice[pos_to], false, &env.space.bound);
        }
        let (removed, added) = env.update_edges(pos_to);
        // Times 2 to represent the symmetric edge
        2. * (added as f32 - removed as f32) / env.neighbourhood.n_neighs() as f32
    }
    
    pub fn chemotaxis_bias<G, B: Boundary<Coord = f32>>(
        &self,
        cell: &Cell<G>,
        pos_to: Pos<usize>,
        chemotaxis_mu: f32,
        bound: &B
    ) -> f32 {
        let (dx1, dy1) = bound.displacement(
            cell.center,
            Pos::new(pos_to.x as f32, pos_to.y as f32)
        );
        let (dx2, dy2) = bound.displacement(
            cell.center,
            cell.light_center
        );

        let dot = dx1 * dx2 + dy1 * dy2;
        let norm1_sq = dx1 * dx1 + dy1 * dy1;
        let norm2_sq = dx2 * dx2 + dy2 * dy2;
        let denom = (norm1_sq * norm2_sq).sqrt();

        if denom <= 0. {
            0.
        } else {
            -chemotaxis_mu * (dot / denom)
        }
    }

    pub fn accept_site_copy(&self, rng: &mut impl Rng, delta_h: f32) -> bool {
        delta_h < 0. || rng.random::<f32>() < E.powf(-delta_h / self.boltz_t)
    }

    pub fn delta_hamiltonian<'a, G: 'a>(
        &self,
        entity_from: LatticeEntity<&RelCell<G>>,
        entity_to: LatticeEntity<&RelCell<G>>,
        neigh_entities: impl Iterator<Item = LatticeEntity<&'a RelCell<G>>>
    ) -> f32 {
        let mut delta_h = 0.;
        delta_h += self.delta_hamiltonian_size(entity_from, entity_to);
        delta_h += self.delta_hamiltonian_adhesion(entity_from, entity_to, neigh_entities);
        delta_h
    }
    
    pub fn delta_hamiltonian_size<G>(
        &self, 
        entity_from: LatticeEntity<&RelCell<G>>, 
        entity_to: LatticeEntity<&RelCell<G>>
    ) -> f32 {
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
    pub fn delta_hamiltonian_adhesion<'a, G: 'a>(
        &self,
        entity_from: LatticeEntity<&RelCell<G>>,
        entity_to: LatticeEntity<&RelCell<G>>,
        neigh_entities: impl Iterator<Item = LatticeEntity<&'a RelCell<G>>>
    ) -> f32 {
        let mut energy = 0.;
        for neigh in neigh_entities {
            energy -= self.adhesion.adhesion_energy(entity_to, neigh);
            energy += self.adhesion.adhesion_energy(entity_from, neigh);
        }
        energy
    }

    pub fn size_energy_diff(&self, area_increased: bool, area: u32, target_area: u32) -> f32 {
        let delta_area = if area_increased { 1. } else { -1. };
        2. * self.size_lambda * delta_area * (area as f32 - target_area as f32) + self.size_lambda
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adhesion::ClonalAdhesion;
    use crate::cell::Cell;
    use crate::genome::MockGenome;
    use crate::io::parameters::StaticAdhesionParameters;

    #[test]
    fn test_delta_hamiltonian_size() {
        let adh_params = StaticAdhesionParameters {
            cell_energy: 10.,
            medium_energy: 20.,
            solid_energy: 20.
        };
        let ca = CellularAutomata::new(
            16., 
            1.,
            1.,
            ClonalAdhesion::new(adh_params, 10)
        );
        let cell = RelCell::mock(Cell::new_empty(100, MockGenome::new(0)));
        let dh = ca.delta_hamiltonian_size(SomeCell(&cell), SomeCell(&cell.clone()));
        assert_eq!(dh, 2.);
    }
}