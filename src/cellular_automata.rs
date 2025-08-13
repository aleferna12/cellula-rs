use crate::adhesion::AdhesionSystem;
use crate::cell::{CanMigrate, Cellular, ChemSniffer, RelCell};
use crate::environment::{EdgesUpdate, Environment};
use crate::environment::LatticeEntity;
use crate::environment::LatticeEntity::{Medium, Solid, SomeCell};
use crate::positional::boundary::{AsLatticeBoundary, Boundary};
use crate::positional::neighbourhood::Neighbourhood;
use crate::positional::pos::Pos;
use rand::Rng;
use std::f32::consts::E;
use crate::constants::Spin;

// This could be a module but it's convenient to be able to access the relevant parameters 
// Also we might eventually want to implement multiple CA choices, in which case I can "easily" make CA a trait 
// that just implements `step()`
pub struct CellularAutomata<A> {
    pub boltz_t: f32,
    pub size_lambda: f32,
    pub chemotaxis_mu: f32,
    pub adhesion: A
}

impl<A> CellularAutomata<A> {
    pub fn new(boltz_t: f32, size_lambda: f32, chemotaxis_mu: f32, adhesion: A) -> Self {
        Self {
            boltz_t,
            size_lambda,
            chemotaxis_mu,
            adhesion
        }
    }

    pub fn chemotaxis_bias<B: Boundary<Coord = f32>>(
        &self,
        cell: &(impl CanMigrate + ChemSniffer),
        pos_to: Pos<usize>,
        chemotaxis_mu: f32,
        bound: &B
    ) -> f32 {
        let (dx1, dy1) = bound.displacement(
            cell.center(),
            Pos::new(pos_to.x as f32, pos_to.y as f32)
        );
        let (dx2, dy2) = bound.displacement(
            cell.center(),
            cell.chem_center()
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

    pub fn delta_hamiltonian_size<C: Cellular>(
        &self,
        entity_source: LatticeEntity<&RelCell<C>>,
        entity_target: LatticeEntity<&RelCell<C>>
    ) -> f32 {
        let mut delta_h = 0.;
        if let SomeCell(cell) = entity_source {
            delta_h += self.size_energy_diff(true, cell.area(), cell.target_area());
        }
        if let SomeCell(cell) = entity_target {
            delta_h += self.size_energy_diff(false, cell.area(), cell.target_area());
        }
        delta_h
    }

    pub fn accept_site_copy(&self, rng: &mut impl Rng, delta_h: f32) -> bool {
        delta_h < 0. || rng.random::<f32>() < E.powf(-delta_h / self.boltz_t)
    }

    pub fn size_energy_diff(&self, area_increased: bool, area: u32, target_area: u32) -> f32 {
        let delta_area = if area_increased { 1. } else { -1. };
        2. * self.size_lambda * delta_area * (area as f32 - target_area as f32) + self.size_lambda
    }
}

impl<A: AdhesionSystem> CellularAutomata<A> {
    pub fn step(
        &self, 
        env: &mut Environment<
            impl CanMigrate + ChemSniffer, 
            impl Neighbourhood, 
            impl AsLatticeBoundary<Coord = f32>
        >, 
        rng: &mut impl Rng
    ) {
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
    /// The number of extra updates that the copy attempt incurred.
    pub fn attempt_site_copy(
        &self,
        env: &mut Environment<
            impl CanMigrate + ChemSniffer, 
            impl Neighbourhood, 
            impl AsLatticeBoundary<Coord = f32>
        >,
        rng: &mut impl Rng,
        pos_source: Pos<usize>,
        pos_target: Pos<usize>
    ) -> f32 {
        let spin_target = env.space.cell_lattice[pos_target];
        if spin_target == Solid.discriminant() {
            return 0.;
        }
        // If was going to copy from a Solid, create a Medium cell instead 
        let spin_source = {
            let spin = env.space.cell_lattice[pos_source];
            if spin == Solid.discriminant() { Medium.discriminant() } else { spin }
        };

        let entity_source = env.cells.get_entity(spin_source);
        let entity_target = env.cells.get_entity(spin_target);
        let neigh_entities = env.space.lat_bound.valid_positions(
            env.neighbourhood.neighbours(pos_target.to_isize())
        ).map(|neigh| {
            env.cells.get_entity(env.space.cell_lattice[neigh.to_usize()])
        });

        let mut delta_h = self.delta_hamiltonian(entity_source, entity_target, neigh_entities);
        if let SomeCell(cell) = entity_source {
            if env.cells.migrate && cell.is_migrating() {
                delta_h += self.chemotaxis_bias(&cell.cell, pos_target, self.chemotaxis_mu, &env.space.bound);
            }
        }
        if !self.accept_site_copy(rng, delta_h) {
            return 0.;
        }
        let edges_update = self.shift_position(
            pos_target,
            spin_target,
            spin_source,
            env,
        );
        // Times 2 to represent the symmetric edge
        2. * (edges_update.added as f32 - edges_update.removed as f32) / env.neighbourhood.n_neighs() as f32
    }

    pub fn shift_position<C: ChemSniffer>(
        &self,
        pos: Pos<usize>,
        from: Spin,
        to: Spin,
        env: &mut Environment<
            C,
            impl Neighbourhood,
            impl AsLatticeBoundary<Coord = f32>
        >,
    ) -> EdgesUpdate {
        // Executes the copy
        env.space.cell_lattice[pos] = to;
        let chem_at = env.space.chem_lattice[pos] as f32;
        if let SomeCell(cell) = env.cells.get_entity_mut(to) {
            cell.shift_position(pos, true, &env.space.bound);
            cell.shift_chem(pos, chem_at, true, &env.space.bound);
        }
        if let SomeCell(cell) = env.cells.get_entity_mut(from) {
            cell.shift_position(pos, false, &env.space.bound);
            cell.shift_chem(pos, chem_at, false, &env.space.bound);
        }
        env.update_edges(pos)
    }

    pub fn delta_hamiltonian<'a, C: 'a + Cellular>(
        &self,
        entity_source: LatticeEntity<&RelCell<C>>,
        entity_target: LatticeEntity<&RelCell<C>>,
        neigh_entities: impl Iterator<Item = LatticeEntity<&'a RelCell<C>>>
    ) -> f32 {
        let mut delta_h = 0.;
        delta_h += self.delta_hamiltonian_size(entity_source, entity_target);
        delta_h += self.delta_hamiltonian_adhesion(entity_source, entity_target, neigh_entities);
        delta_h
    }
    
    // TODO!: test
    pub fn delta_hamiltonian_adhesion<'a, C: 'a>(
        &self,
        entity_source: LatticeEntity<&RelCell<C>>,
        entity_target: LatticeEntity<&RelCell<C>>,
        neigh_entities: impl Iterator<Item = LatticeEntity<&'a RelCell<C>>>
    ) -> f32 {
        let mut energy = 0.;
        for neigh in neigh_entities {
            energy -= self.adhesion.adhesion_energy(entity_target, neigh);
            energy += self.adhesion.adhesion_energy(entity_source, neigh);
        }
        energy
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::adhesion::{ClonalAdhesion, StaticAdhesion};
    use crate::cell::Cell;
    use crate::genome::MockGenome;

    #[test]
    fn test_delta_hamiltonian_size() {
        let adh = StaticAdhesion {
            cell_energy: 10.,
            medium_energy: 20.,
            solid_energy: 20.
        };
        let ca = CellularAutomata::new(
            16., 
            1.,
            1.,
            ClonalAdhesion::new(10, adh)
        );
        let cell = RelCell::mock(Cell::new_empty(
            100, 
            200, 
            MockGenome::new(0)
        ));
        let dh = ca.delta_hamiltonian_size(SomeCell(&cell), SomeCell(&cell.clone()));
        assert_eq!(dh, 2.);
    }
}