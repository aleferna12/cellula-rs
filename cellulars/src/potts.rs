//! Contains [`Potts`] algorithms.

use crate::constants::FloatType;
use crate::positional::neighborhood::Neighborhood;
use crate::prelude::{AdhesionSystem, Cellular, CopyBias, NoBias, Habitable, Pos, Spin};
use rand::RngExt;
#[cfg(not(feature = "f64"))]
use std::f32::consts::E;
#[cfg(feature = "f64")]
use std::f64::consts::E;

/// This type is a Potts algorithm, which can modify a [`Habitable`] `H`.
pub trait Potts<H> {
    /// Executes a Monte Carlo step of the simulation by updating `hab`.
    fn step(
        &mut self,
        hab: &mut H,
        rng: &mut impl RngExt
    );
}

/// This potts algorithm runs
#[derive(Debug, Clone, PartialEq)]
pub struct EdgePotts<A, B = NoBias> {
    /// Adhesion system used for copy attempts.
    pub adhesion: A,
    /// Bias applied to copy attempts.
    pub bias: B,
    /// Boltzmann temperature of the system.
    pub boltz_t: FloatType,
    /// Strength of the size constraint on the energy functional.
    pub size_lambda: FloatType
}

impl<A, B> EdgePotts<A, B> {
    /// Returns the energy differential associated with the size constraint if `spin_source`
    /// were to be copied into `spin_target`.
    fn delta_hamiltonian_size<H: Habitable<Cell = impl Cellular>>(
        &self,
        spin_source: Spin,
        spin_target: Spin,
        env: &H,
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
    fn accept_site_copy(&self, rng: &mut impl RngExt, delta_h: FloatType) -> bool {
        delta_h < 0. || rng.random::<FloatType>() < E.powf(-delta_h / self.boltz_t)
    }

    /// Returns the energy differential resulting from an atomic increase or decrease of `area`.
    fn size_energy_diff(&self, area_increased: bool, area: u32, target_area: u32) -> FloatType {
        let delta_area = if area_increased { 1. } else { -1. };
        2. * self.size_lambda * delta_area * (area as FloatType - target_area as FloatType) + self.size_lambda
    }

    /// Returns the total energy differential of the system if `spin_source` were to be copied into `spin_target`.
    fn delta_hamiltonian<H: Habitable<Cell = impl Cellular>>(
        &self,
        spin_source: Spin,
        spin_target: Spin,
        neigh_spins: impl IntoIterator<Item = Spin>,
        env: &H,
    ) -> FloatType
    where A: AdhesionSystem<H> {
        self.delta_hamiltonian_size(spin_source, spin_target, env)
            + self.delta_hamiltonian_adhesion(spin_source, spin_target, neigh_spins, env)
    }

    /// Returns the energy differential associated with adhesion if `spin_source` were to be copied into `spin_target`.
    fn delta_hamiltonian_adhesion<H: Habitable<Cell = impl Cellular>>(
        &self,
        spin_source: Spin,
        spin_target: Spin,
        neigh_spins: impl IntoIterator<Item = Spin>,
        env: &H,
    ) -> FloatType
    where
        A: AdhesionSystem<H> {
        let mut energy = 0.;
        for neigh in neigh_spins {
            energy -= self.adhesion.adhesion_energy(
                spin_target,
                neigh,
                env
            );
            energy += self.adhesion.adhesion_energy(
                spin_source,
                neigh,
                env
            );
        }
        energy
    }

    /// Attempts to execute the selected site copy from `pos_source` into `pos_target`.
    ///
    /// Returns the number of extra updates that the copy attempt should incur
    /// based on how many cell edges it modified.
    fn attempt_site_copy<H: Habitable<Cell = impl Cellular>>(
        &self,
        pos_source: Pos<usize>,
        pos_target: Pos<usize>,
        hab: &mut H,
        rng: &mut impl RngExt
    ) -> FloatType
    where
        A: AdhesionSystem<H>,
        B: CopyBias<H> {
        let spin_target = hab.env().cell_lattice[pos_target];
        if spin_target == Spin::Solid {
            return 0.;
        }
        let spin_source = {
            let spin = hab.env().cell_lattice[pos_source];
            // If was going to copy from a Solid, treat it as a Medium cell instead
            if let Spin::Solid = spin {
                Spin::Medium
            } else {
                spin
            }
        };
        let neigh_spins = hab
            .env()
            .valid_neighbors(pos_target)
            .map(|pos| hab.env().cell_lattice[pos]);

        let delta_h = self.delta_hamiltonian(
            spin_source,
            spin_target,
            neigh_spins,
            hab
        ) + self.bias.bias(
            pos_source,
            pos_target,
            hab
        );

        if !self.accept_site_copy(rng, delta_h) {
            return 0.;
        }
        let edges_update = hab.transfer_position(
            pos_target,
            spin_source
        );
        // Times 2 to represent the symmetric edge
        2. * (edges_update.added as FloatType - edges_update.removed as FloatType) / hab.env().neighborhood.n_neighs() as FloatType
    }
}

impl<A, B, H, C> Potts<H> for EdgePotts<A, B>
where
    A: AdhesionSystem<H>,
    B: CopyBias<H>,
    H: Habitable<Cell = C>,
    C: Cellular {
    fn step(&mut self, hab: &mut H, rng: &mut impl RngExt) {
        let mut to_visit = 2. * hab.env().edge_book.len() as FloatType / hab.env().neighborhood.n_neighs() as FloatType;
        while 0. < to_visit {
            let edge_i = hab.env().edge_book.random_index(rng);
            let edge = hab.env().edge_book.at(edge_i);
            // This is WAY faster than keeping the symmetric edge in EdgeBook (like 2x as fast!)
            // or at least, this is the case when using IndexSet, I would assume its somewhat implementation-dependent
            let (pos_source, pos_target) = if rng.random::<FloatType>() < 0.5 {
                (edge.p1, edge.p2)
            } else {
                (edge.p2, edge.p1)
            };
            to_visit += self.attempt_site_copy(pos_source, pos_target, hab, rng);
            to_visit -= 1.;
        }
    }
}