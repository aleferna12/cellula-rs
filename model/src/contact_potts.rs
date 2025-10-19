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
use rand::Rng;
use cellulars_lib::habitable::Habitable;

// This could be a module but it's convenient to be able to access the relevant parameters
// Also we might eventually want to implement multiple CA choices, in which case I can "easily" make CA a trait 
// that just implements `step()`
#[derive(Clone, Builder)]
pub struct ContactPotts {
    pub boltz_t: f32,
    pub size_lambda: f32,
    pub perimeter_lambda: f32,
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

    fn perimeter_energy_diff(&self, delta_perimeter: i32, perimeter: u32, target_perimeter: u32) -> f32 {
        2. * self.perimeter_lambda * delta_perimeter as f32 * (perimeter as f32 - target_perimeter as f32) + self.perimeter_lambda
    }

    fn delta_hamiltonian_perimeter(&self, spin_source: Spin, spin_target: Spin, env: &ChemEnvironment) -> f32 {
        let mut delta_h = 0.;
        if let Spin::Some(cell_index) = spin_source {
            let cell = env.env().cells.get_cell(cell_index);
            delta_h += self.perimeter_energy_diff(cell.delta_perimeter, cell.perimeter, cell.target_perimeter);
        }
        if let Spin::Some(cell_index) = spin_target {
            let cell = env.env().cells.get_cell(cell_index);
            delta_h += self.perimeter_energy_diff(cell.delta_perimeter, cell.perimeter, cell.target_perimeter);
        }
        delta_h
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

        let neighs_target = env.env().valid_neighbours(pos_target).map(|neigh| {
            env.env().cell_lattice[neigh]
        }).collect::<Vec<_>>();

        if let Spin::Some(cell_index) = spin_target {
            let cell_target = env.env_mut().cells.get_cell_mut(cell_index);
            cell_target.delta_perimeter = 0;
            for spin in &neighs_target {
                if spin == &spin_target {
                    cell_target.delta_perimeter += 1;
                } else {
                    cell_target.delta_perimeter -= 1;
                }
            }
        }

        if let Spin::Some(cell_index) = spin_source {
            let neighs_source = env.env().valid_neighbours(pos_source).map(|neigh| {
                env.env().cell_lattice[neigh]
            }).collect::<Vec<_>>();
            let cell_source = env.env_mut().cells.get_cell_mut(cell_index);
            cell_source.delta_perimeter = 0;
            for spin in neighs_source {
                if spin == spin_source {
                    cell_source.delta_perimeter -= 1;
                } else {
                    cell_source.delta_perimeter += 1;
                }
            }
        }

        let delta_h = self.delta_hamiltonian(
            spin_source,
            spin_target,
            neighs_target,
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
            + self.delta_hamiltonian_perimeter(spin_source, spin_target, env)
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