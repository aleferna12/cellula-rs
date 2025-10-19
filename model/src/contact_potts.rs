use crate::chem_environment::ChemEnvironment;
use bon::Builder;
use cellulars_lib::adhesion::{AdhesionSystem, StaticAdhesion};
use cellulars_lib::basic_cell::Cellular;
use cellulars_lib::positional::boundaries::Boundary;
use cellulars_lib::positional::neighbourhood::Neighbourhood;
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::potts::Potts;
use cellulars_lib::spin::Spin;
use std::hint::black_box;

// This could be a module but it's convenient to be able to access the relevant parameters
// Also we might eventually want to implement multiple CA choices, in which case I can "easily" make CA a trait 
// that just implements `step()`
#[derive(Clone, Builder)]
pub struct ContactPotts {
    pub boltz_t: f32,
    pub size_lambda: f32,
    pub chemotaxis_mu: f32,
    pub act_lambda: f32,
    pub enable_migration: bool,
    pub adhesion: StaticAdhesion
}

impl ContactPotts {
    fn migration_bias(&self, pos_source: Pos<usize>, pos_target: Pos<usize>, env: &ChemEnvironment) -> f32 {
        let Spin::Some(cell_index) = env.cell_lattice[pos_source] else {
            return 0.;
        };
        let cell = env.cells.get_cell(cell_index);
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

    fn act_bias(&self, pos_source: Pos<usize>, pos_target: Pos<usize>, env: &ChemEnvironment) -> f32 {
        let act_source = Self::mean_act(pos_source, env);
        let act_target = Self::mean_act(pos_target, env);
            -self.act_lambda / env.act_max as f32 * (act_source - act_target)
    }

    fn mean_act(pos: Pos<usize>, env: &ChemEnvironment) -> f32 {
        // We do precompute these positions for pos_target in `attempt_site_copy`
        // Reusing that computation could be slightly faster
        let mut acts = env.bounds.lattice_boundary.valid_positions(
            env.neighbourhood.neighbours(pos.to_isize())
        ).map(|neigh| {
            env.act_lattice[neigh.to_usize()]
        });
        if acts.any(|x| x == 0) {
            return 0.
        }
        // This could be wrong if there are invalid pos in the neighbourhood
        let root = env.neighbourhood.n_neighs() as f32;
        let sum: f32 = acts.map(|x| (x as f32).ln()).sum();
        black_box((sum / root).exp())
    }
}

impl Potts for ContactPotts {
    type Environment = ChemEnvironment;

    fn boltz_t(&self) -> f32 {
        self.boltz_t
    }

    fn size_lambda(&self) -> f32 {
        self.size_lambda
    }

    fn copy_biases(&self, pos_source: Pos<usize>, pos_target: Pos<usize>, env: &Self::Environment) -> f32 {
        let mut biases = self.act_bias(pos_source, pos_target, env);
        if self.enable_migration {
            biases += self.migration_bias(pos_source, pos_target, env);
        }
        biases
    }

    fn delta_hamiltonian_adhesion(
        &self, 
        spin_source: Spin, 
        spin_target: Spin,
        neigh_spin: impl IntoIterator<Item=Spin>,
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