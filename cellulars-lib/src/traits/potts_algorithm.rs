//! Contains logic associated with [`PottsAlgorithm`].

use crate::constants::FloatType;
use crate::positional::neighbourhood::Neighbourhood;
use crate::positional::pos::Pos;
use crate::spin::Spin;
use crate::traits::cellular::Cellular;
use crate::traits::habitable::Habitable;
use rand::Rng;
#[cfg(not(feature = "high-precision"))]
use std::f32::consts::E;
#[cfg(feature = "high-precision")]
use std::f64::consts::E;

/// This trait defines how a Monte Carlo [`PottsAlgorithm::step()`] of the model should modify a
///[`Habitable`] environment.
///
/// The default methods for this trait restrict lattice updates to cell borders following
/// [van Steijn, 2022](https://doi.org/10.1371/journal.pcbi.1009156) to improve computational efficiency.
pub trait PottsAlgorithm {
    /// Type of environment that is going to be modified each [`PottsAlgorithm::step()`].
    type Environment: Habitable;

    /// Returns the Boltzmann temperature of the system.
    fn boltz_t(&self) -> FloatType;

    /// Returns the scaling constant associated with the penalty given to size deviations.
    fn size_lambda(&self) -> FloatType;

    /// Returns the energy differential associated with copy biases of the model.
    ///
    /// Returns 0 by default.
    ///
    /// Overriding this method allows to easily extend the model's behaviour
    /// without having to override [`PottsAlgorithm::attempt_site_copy()`].
    fn copy_biases(&self, _pos_source: Pos<usize>, _pos_target: Pos<usize>, _env: &Self::Environment) -> FloatType {
        0.
    }

    /// Returns the energy differential associated with the size constraint if `spin_source`
    /// were to be copied into `spin_target`.
    fn delta_hamiltonian_size(
        &self,
        spin_source: Spin,
        spin_target: Spin,
        env: &Self::Environment,
    ) -> FloatType {
        let mut delta_h = 0.;
        if let Spin::Some(cell_index) = spin_source {
            let rel_cell = &env.env().cells[cell_index];
            delta_h += self.size_energy_diff(true, rel_cell.cell.area(), rel_cell.cell.target_area());
        }
        if let Spin::Some(cell_index) = spin_target {
            let rel_cell = &env.env().cells[cell_index];
            delta_h += self.size_energy_diff(false, rel_cell.cell.area(), rel_cell.cell.target_area());
        }
        delta_h
    }

    /// Returns whether a copy attempt that results in an energy differential `delta_h` should be randomly accepted
    /// by drawing from a Boltzmann distribution.
    fn accept_site_copy(&self, rng: &mut impl Rng, delta_h: FloatType) -> bool {
        delta_h < 0. || rng.random::<FloatType>() < E.powf(-delta_h / self.boltz_t())
    }

    /// Returns the energy differential resulting from an atomic increase or decrease of `area`.
    fn size_energy_diff(&self, area_increased: bool, area: u32, target_area: u32) -> FloatType {
        let delta_area = if area_increased { 1. } else { -1. };
        2. * self.size_lambda() * delta_area * (area as FloatType - target_area as FloatType) + self.size_lambda()
    }

    /// Executes a Monte Carlo step of the simulation by updating `env`.
    fn step(
        &self,
        env: &mut Self::Environment,
        rng: &mut impl Rng
    ) {
        let mut to_visit = 2. * env.env().edge_book.len() as FloatType / env.env().neighbourhood.n_neighs() as FloatType;
        while 0. < to_visit {
            let edge_i = env.env().edge_book.random_index(rng);
            let edge = env.env().edge_book.at(edge_i);
            // This is WAY faster than keeping the symmetric edge in EdgeBook (like 2x as fast!)
            // or at least, this is the case when using IndexSet, I would assume its somewhat implementation-dependent
            let (pos_source, pos_target) = if rng.random::<FloatType>() < 0.5 {
                (edge.p1, edge.p2)
            } else {
                (edge.p2, edge.p1)
            };
            to_visit += self.attempt_site_copy(pos_source, pos_target, env, rng);
            to_visit -= 1.;
        }
    }

    /// Attempts to execute the selected site copy from `pos_source` into `pos_target`.
    ///
    /// Returns the number of extra updates that the copy attempt should incur
    /// based on how many cell edges it modified.
    fn attempt_site_copy(
        &self,
        pos_source: Pos<usize>,
        pos_target: Pos<usize>,
        env: &mut Self::Environment,
        rng: &mut impl Rng
    ) -> FloatType {
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
        2. * (edges_update.added as FloatType - edges_update.removed as FloatType) / env.env().neighbourhood.n_neighs() as FloatType
    }

    /// Returns the total energy differential of the system if `spin_source` were to be copied into `spin_target`.
    fn delta_hamiltonian(
        &self,
        spin_source: Spin,
        spin_target: Spin,
        neigh_spins: impl IntoIterator<Item = Spin>,
        env: &Self::Environment,
    ) -> FloatType {
        self.delta_hamiltonian_size(spin_source, spin_target, env)
            + self.delta_hamiltonian_adhesion(spin_source, spin_target, neigh_spins, env)
    }

    /// Returns the energy differential associated with adhesion if `spin_source` were to be copied into `spin_target`.
    fn delta_hamiltonian_adhesion(
        &self,
        spin_source: Spin,
        spin_target: Spin,
        neigh_spin: impl IntoIterator<Item = Spin>,
        env: &Self::Environment,
    ) -> FloatType;
}