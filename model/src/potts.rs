//! Contains logic associated with [Potts].

use crate::cell::CellType;
use crate::environment::Environment;
use bon::Builder;
use cellulars_lib::positional::boundaries::Boundary;
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::spin::Spin;
use cellulars_lib::static_adhesion::StaticAdhesion;
use cellulars_lib::traits::adhesion_system::AdhesionSystem;
use cellulars_lib::traits::cellular::Cellular;
use cellulars_lib::traits::potts_algorithm::PottsAlgorithm;

// This could be a module but it's convenient to be able to access the relevant parameters
// Also we might eventually want to implement multiple CA choices, in which case I can "easily" make CA a trait 
// that just implements `step()`
/// A Potts model that implements cell migration.
#[derive(Clone, Builder)]
pub struct Potts {
    /// Boltz temperature of the model.
    pub boltz_t: f32,
    /// Scaler constant associated with the penalty for size deviations.
    pub size_lambda: f32,
    /// Scaler constant associated with the speed of migration.
    pub chemotaxis_mu: f32,
    /// Whether we allow cell migration.
    pub enable_migration: bool,
    /// Adhesion system used in [Potts::delta_hamiltonian_adhesion()].
    pub adhesion: StaticAdhesion
}

impl PottsAlgorithm for Potts {
    type Environment = Environment;

    fn boltz_t(&self) -> f32 {
        self.boltz_t
    }

    fn size_lambda(&self) -> f32 {
        self.size_lambda
    }

    fn copy_biases(&self, pos_source: Pos<usize>, pos_target: Pos<usize>, env: &Self::Environment) -> f32 {
        if !self.enable_migration {
            return 0.
        }
        let Spin::Some(cell_index) = env.base_env.cell_lattice[pos_source] else {
            return 0.;
        };
        let cell = env.base_env.cells.get_cell(cell_index);
        if let CellType::Dividing = cell.cell_type {
            return 0.;
        }

        let (dx1, dy1) = env.base_env.bounds.boundary.displacement(
            cell.center(),
            Pos::new(pos_target.x as f32, pos_target.y as f32)
        );
        let (dx2, dy2) = env.base_env.bounds.boundary.displacement(
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
            -self.chemotaxis_mu * (dot / denom)
        }
    }

    fn delta_hamiltonian_adhesion(
        &self, 
        spin_source: Spin, 
        spin_target: Spin,
        neigh_spin: impl IntoIterator<Item = Spin>,
        _env: &Self::Environment
    ) -> f32 {
        let mut energy = 0.;
        for neigh in neigh_spin {
            energy -= self.adhesion.adhesion_energy(
                spin_target,
                neigh,
                &()
            );
            energy += self.adhesion.adhesion_energy(
                spin_source,
                neigh,
                &()
            );
        }
        energy
    }
}