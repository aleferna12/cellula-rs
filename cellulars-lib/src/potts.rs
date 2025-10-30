use crate::basic_cell::Cellular;
use crate::habitable::Habitable;
use crate::positional::neighbourhood::Neighbourhood;
use crate::positional::pos::Pos;
use crate::spin::Spin;
use rand::Rng;
use std::f32::consts::E;

pub trait Potts {
    type Environment: Habitable;
    
    fn boltz_t(&self) -> f32;
    
    fn size_lambda(&self) -> f32;
    
    fn copy_biases(&self, _pos_source: Pos<usize>, _pos_target: Pos<usize>, _env: &Self::Environment) -> f32 {
        0.
    }

    fn delta_hamiltonian_size(
        &self,
        spin_source: Spin,
        spin_target: Spin,
        env: &Self::Environment,
    ) -> f32 {
        let mut delta_h = 0.;
        if let Spin::Some(cell_index) = spin_source {
            let cell = env.env().cells.get_cell(cell_index);
            delta_h += self.size_energy_diff(true, cell.area(), cell.target_area());
        }
        if let Spin::Some(cell_index) = spin_target {
            let cell = env.env().cells.get_cell(cell_index);
            delta_h += self.size_energy_diff(false, cell.area(), cell.target_area());
        }
        delta_h
    }

    fn accept_site_copy(&self, rng: &mut impl Rng, delta_h: f32) -> bool {
        delta_h < 0. || rng.random::<f32>() < E.powf(-delta_h / self.boltz_t())
    }

    fn size_energy_diff(&self, area_increased: bool, area: u32, target_area: u32) -> f32 {
        let delta_area = if area_increased { 1. } else { -1. };
        2. * self.size_lambda() * delta_area * (area as f32 - target_area as f32) + self.size_lambda()
    }

    fn step(
        &self,
        env: &mut Self::Environment,
        rng: &mut impl Rng
    ) {
        let mut to_visit = 2. * env.env().edge_book.len() as f32 / env.env().neighbourhood.n_neighs() as f32;
        while 0. < to_visit {
            let edge_i = env.env().edge_book.random_index(rng);
            let edge = env.env().edge_book.at(edge_i);
            // This is WAY faster than keeping the symmetric edge in EdgeBook (like 2x as fast!)
            // or at least, this is the case when using IndexSet, I would assume its somewhat implementation-dependent
            let (pos_source, pos_target) = if rng.random::<f32>() < 0.5 {
                (edge.p1, edge.p2)
            } else {
                (edge.p2, edge.p1)
            };
            to_visit += self.attempt_site_copy(pos_source, pos_target, env, rng);
            to_visit -= 1.;
        }
    }

    /// Attempts to execute the selected site copy.
    ///
    /// # Returns:
    ///
    /// The number of extra updates that the copy attempt incurred.
    fn attempt_site_copy(
        &self,
        pos_source: Pos<usize>,
        pos_target: Pos<usize>,
        env: &mut Self::Environment,
        rng: &mut impl Rng
    ) -> f32 {
        let spin_target = env.env().cell_lattice[pos_target];
        if spin_target == Spin::Solid {
            return 0.;
        }
        let spin_source = {
            let spin = env.env().cell_lattice[pos_source];
            // If was going to copy from a Solid, treat it as a Medium cell instead
            if let Spin::Solid = spin {
                Spin::Medium
            } else {
                spin
            }
        };
        let neigh_spins = env
            .env()
            .valid_neighbours(pos_target)
            .map(|pos| env.env().cell_lattice[pos]);

        let delta_h = self.delta_hamiltonian(
            spin_source, 
            spin_target,
            neigh_spins,
            env
        ) + self.copy_biases(
            pos_source,
            pos_target,
            env
        );
        
        if !self.accept_site_copy(rng, delta_h) {
            return 0.;
        }
        let edges_update = env.grant_position(
            pos_target,
            spin_source
        );
        // Times 2 to represent the symmetric edge
        2. * (edges_update.added as f32 - edges_update.removed as f32) / env.env().neighbourhood.n_neighs() as f32
    }

    fn delta_hamiltonian(
        &self,
        spin_source: Spin,
        spin_target: Spin,
        neigh_spins: impl IntoIterator<Item = Spin>,
        env: &Self::Environment,
    ) -> f32 {
        self.delta_hamiltonian_size(spin_source, spin_target, env)
            + self.delta_hamiltonian_adhesion(spin_source, spin_target, neigh_spins, env)
    }

    fn delta_hamiltonian_adhesion(
        &self,
        spin_source: Spin,
        spin_target: Spin,
        neigh_spin: impl IntoIterator<Item = Spin>,
        env: &Self::Environment,
    ) -> f32;
}