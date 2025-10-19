use crate::chem_environment::ChemEnvironment;
use crate::clonal_adhesion::ClonalAdhesion;
use bon::Builder;
use cellulars_lib::adhesion::AdhesionSystem;
use cellulars_lib::basic_cell::Cellular;
use cellulars_lib::positional::boundaries::Boundary;
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::potts::Potts;
use cellulars_lib::spin::Spin;

// This could be a module but it's convenient to be able to access the relevant parameters
// Also we might eventually want to implement multiple CA choices, in which case I can "easily" make CA a trait 
// that just implements `step()`
#[derive(Clone, Builder)]
pub struct ClonalPotts {
    pub boltz_t: f32,
    pub size_lambda: f32,
    pub chemotaxis_mu: f32,
    pub enable_migration: bool,
    pub adhesion: ClonalAdhesion
}

impl Potts for ClonalPotts {
    type Environment = ChemEnvironment;

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
        let Spin::Some(cell_index) = env.cell_lattice[pos_source] else {
            return 0.;
        };
        let cell = env.cells.get_cell(cell_index);
        if !cell.is_migrating() {
            return 0.;
        }

        let (dx1, dy1) = env.bounds.boundary.displacement(
            cell.center(),
            Pos::new(pos_target.x as f32, pos_target.y as f32)
        );
        let (dx2, dy2) = env.bounds.boundary.displacement(
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
        env: &Self::Environment
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