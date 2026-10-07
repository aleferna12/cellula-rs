//! Contains logic associated with [`Potts`].

use crate::my_cell::CellType;
use crate::my_environment::MyEnvironment;
use bon::Builder;
use cellulars_lib::constants::FloatType;
use cellulars_lib::perimeter_constraint::PerimeterConstraint;
use cellulars_lib::positional::boundaries::Boundary;
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::spin::Spin;
use cellulars_lib::traits::adhesion_system::AdhesionSystem;
use cellulars_lib::traits::cellular::Cellular;
use cellulars_lib::traits::potts_algorithm::PottsAlgorithm;
use crate::pairwise_adhesion::PairwiseAdhesion;

// This could be a module but it's convenient to be able to access the relevant parameters
// Also we might eventually want to implement multiple CA choices, in which case I can "easily" make CA a trait 
// that just implements `step()`
/// A Potts model that implements cell migration.
#[derive(Clone, Builder)]
pub struct Potts {
    /// Boltz temperature of the model.
    pub boltz_t: FloatType,
    /// Scaler constant associated with the penalty for size deviations.
    pub size_lambda: FloatType,
    /// Scaler constant associated with the speed of migration.
    pub chemotaxis_mu: FloatType,
    /// Scaler constant associated with the speed of persistent migration.
    pub persistence_mu: FloatType,
    /// Whether we allow cell migration.
    pub enable_migration: bool,
    /// Adhesion system used in [`Potts::delta_hamiltonian_adhesion()`].
    pub adhesion: PairwiseAdhesion,
    /// Penalty applied to deviations from the cells' target perimeter.
    pub perimeter: PerimeterConstraint
}

impl Potts {
    /// Returns the energy differential associated with the chemotaxis of the cell that owns `pos_source`.
    fn _chemotaxis_bias(&self, pos_source: Pos<usize>, pos_target: Pos<usize>, env: &MyEnvironment) -> FloatType {
        if !self.enable_migration {
            return 0.
        }
        let Spin::Some(cell_index) = env.env.cell_lattice[pos_source] else {
            return 0.;
        };
        let rel_cell = &env.env.cells[cell_index];
        if let CellType::Dividing = rel_cell.cell.cell_type {
            return 0.;
        }

        let (dx1, dy1) = env.env.bounds.boundary.displacement(
            rel_cell.cell.center(),
            Pos::new(pos_target.x as FloatType, pos_target.y as FloatType)
        );
        let (dx2, dy2) = env.env.bounds.boundary.displacement(
            rel_cell.cell.center(),
            rel_cell.cell.chem_center()
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

    /// Returns the energy differential associated with the persistent migration of the cell that owns
    /// `pos_source`.
    ///
    /// Protrusions that grow along the direction the cell has been travelling towards
    /// (see [`MyCell::target_vec()`](crate::my_cell::MyCell::target_vec())) are favoured, which makes the
    /// cell keep moving forward instead of drifting with the chemical gradient alone.
    ///
    /// This is the same energy term as the chemotaxis one, except that the direction the cell pushes
    /// towards is a property of the cell rather than of its perceived chemical field, following
    /// [Colizzi, 2020](https://doi.org/10.7554/eLife.56349).
    fn persistence_bias(&self, pos_source: Pos<usize>, pos_target: Pos<usize>, env: &MyEnvironment) -> FloatType {
        if !self.enable_migration || self.persistence_mu == 0. {
            return 0.
        }
        let Spin::Some(cell_index) = env.env.cell_lattice[pos_source] else {
            return 0.;
        };
        let rel_cell = &env.env.cells[cell_index];
        if let CellType::Dividing = rel_cell.cell.cell_type {
            return 0.;
        }

        // Direction in which the protrusion would grow
        let (dx, dy) = env.env.bounds.boundary.displacement(
            rel_cell.cell.center(),
            Pos::new(pos_target.x as FloatType, pos_target.y as FloatType)
        );
        let norm = dx.hypot(dy);
        // The copy would happen right on top of the cell center
        if norm <= 0. {
            return 0.;
        }
        // The target vector is already normalised, so only the protrusion needs to be
        let target_vec = rel_cell.cell.target_vec();
        -self.persistence_mu * (dx * target_vec.x + dy * target_vec.y) / norm
    }
}

impl PottsAlgorithm for Potts {
    type Environment = MyEnvironment;

    fn boltz_t(&self) -> FloatType {
        self.boltz_t
    }

    fn size_lambda(&self) -> FloatType {
        self.size_lambda
    }

    fn copy_biases(&self, pos_source: Pos<usize>, pos_target: Pos<usize>, env: &Self::Environment) -> FloatType {
        self.persistence_bias(pos_source, pos_target, env)
            // + self.chemotaxis_bias(pos_source, pos_target, env)
            + self.perimeter.energy_diff(pos_source, pos_target, env)
    }

    fn delta_hamiltonian_adhesion(
        &self,
        spin_source: Spin,
        spin_target: Spin,
        neigh_spin: impl IntoIterator<Item = Spin>,
        env: &Self::Environment
    ) -> FloatType {
        let mut energy = 0.;
        for neigh in neigh_spin {
            energy -= self.adhesion.adhesion_energy(
                spin_target,
                neigh,
                &env
            );
            energy += self.adhesion.adhesion_energy(
                spin_source,
                neigh,
                &env
            );
        }
        energy
    }
}
