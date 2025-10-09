use crate::cell::Cell;
use crate::chem_environment::ChemEnvironment;
use crate::clonal_adhesion::ClonalAdhesion;
use bon::Builder;
use cellulars_lib::adhesion::AdhesionSystem;
use cellulars_lib::basic_cell::Cellular;
use cellulars_lib::entity::{Entity, Spin};
use cellulars_lib::environment::Habitable;
use cellulars_lib::positional::boundaries::Boundary;
use cellulars_lib::positional::neighbourhood::Neighbourhood;
use cellulars_lib::positional::pos::Pos;
use rand::Rng;
use std::f32::consts::E;

// This could be a module but it's convenient to be able to access the relevant parameters
// Also we might eventually want to implement multiple CA choices, in which case I can "easily" make CA a trait 
// that just implements `step()`
#[derive(Clone, Builder)]
pub struct CellularAutomata {
    pub boltz_t: f32,
    pub size_lambda: f32,
    pub chemotaxis_mu: f32,
    pub enable_migration: bool,
    pub adhesion: ClonalAdhesion
}

impl CellularAutomata {
    pub fn chemotaxis_bias<B: Boundary<Coord = f32>>(
        &self,
        cell: &Cell,
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

    pub fn delta_hamiltonian_size(
        &self,
        spin_source: Spin,
        spin_target: Spin,
        env: &ChemEnvironment
    ) -> f32 {
        let mut delta_h = 0.;
        if let Spin::Some(cell_index) = spin_source {
            let cell = env.cells.get_cell(cell_index);
            delta_h += self.size_energy_diff(true, cell.area(), cell.target_area());
        }
        if let Entity::Some(cell_index) = spin_target {
            let cell = env.cells.get_cell(cell_index);
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

impl CellularAutomata {
    pub fn step(
        &self, 
        env: &mut ChemEnvironment,
        rng: &mut impl Rng
    ) {
        let mut to_visit = 2. * env.edge_book.len() as f32 / env.neighbourhood.n_neighs() as f32;
        while 0. < to_visit {
            let edge_i = env.edge_book.random_index(rng);
            let edge = env.edge_book.at(edge_i);
            // This is WAY faster than keeping the symmetric edge in EdgeBook (like 2x as fast!)
            // or at least, this is the case when using IndexSet, I would assume its somewhat implementation-dependent
            let (pos_source, pos_target) = if rng.random::<f32>() < 0.5 {
                (edge.p1, edge.p2)
            } else {
                (edge.p2, edge.p1)
            };
            to_visit += self.attempt_site_copy(env, rng, pos_source, pos_target);
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
        env: &mut ChemEnvironment,
        rng: &mut impl Rng,
        pos_source: Pos<usize>,
        pos_target: Pos<usize>
    ) -> f32 {
        let spin_target = env.cell_lattice[pos_target];
        if spin_target == Entity::Solid {
            return 0.;
        }
        // If was going to copy from a Solid, create a Medium cell instead 
        let spin_source = {
            let spin = env.cell_lattice[pos_source];
            if spin == Entity::Solid { Entity::Medium } else { spin }
        };
        let neigh_spins = env.bounds.lattice_boundary.valid_positions(
            env.neighbourhood.neighbours(pos_target.to_isize())
        ).map(|neigh| {
            env.cell_lattice[neigh.to_usize()]
        });

        let mut delta_h = self.delta_hamiltonian(spin_source, spin_target, neigh_spins, env);
        if let Entity::Some(cell_index) = spin_source {
            let cell = env.cells.get_cell(cell_index);
            if self.enable_migration && cell.is_migrating() {
                delta_h += self.chemotaxis_bias(&cell.cell, pos_target, self.chemotaxis_mu, &env.bounds.boundary);
            }
        }
        if !self.accept_site_copy(rng, delta_h) {
            return 0.;
        }
        let edges_update = env.grant_position(
            pos_target,
            spin_source
        );
        // Times 2 to represent the symmetric edge
        2. * (edges_update.added as f32 - edges_update.removed as f32) / env.neighbourhood.n_neighs() as f32
    }

    pub fn delta_hamiltonian(
        &self,
        spin_source: Spin,
        spin_target: Spin,
        neigh_spins: impl Iterator<Item = Spin>,
        env: &ChemEnvironment
    ) -> f32 {
        let mut delta_h = 0.;
        delta_h += self.delta_hamiltonian_size(spin_source, spin_target, env);
        delta_h += self.delta_hamiltonian_adhesion(spin_source, spin_target, neigh_spins, env);
        delta_h
    }
    
    // TODO!: test
    pub fn delta_hamiltonian_adhesion(
        &self,
        spin_source: Spin,
        spin_target: Spin,
        neigh_spin: impl Iterator<Item = Spin>,
        env: &ChemEnvironment,
    ) -> f32 {
        let mut energy = 0.;
        for neigh in neigh_spin {
            energy -= self.adhesion.adhesion_energy(
                spin_target, 
                neigh, 
                &env.clones_table
            );
            energy += self.adhesion.adhesion_energy(
                spin_source, 
                neigh,
                &env.clones_table
            );
        }
        energy
    }
}