use crate::my_environment::MyEnvironment;
use bon::Builder;
use cellulars_lib::adhesion::{AdhesionSystem, StaticAdhesion};
use cellulars_lib::basic_cell::Cellular;
use cellulars_lib::positional::boundaries::Boundary;
use cellulars_lib::positional::neighbourhood::Neighbourhood;
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::potts::Potts;
use cellulars_lib::spin::Spin;
use rand::Rng;
use rand_distr::num_traits::Pow;
use cellulars_lib::constants::CellIndex;
use cellulars_lib::habitable::Habitable;

// This could be a module but it's convenient to be able to access the relevant parameters
// Also we might eventually want to implement multiple CA choices, in which case I can "easily" make CA a trait 
// that just implements `step()`
#[derive(Clone, Builder)]
pub struct MyPotts {
    pub boltz_t: f32,
    pub size_lambda: f32,
    pub perimeter_lambda: f32,
    pub act_lambda: f32,
    pub chemotaxis_mu: f32,
    pub enable_migration: bool,
    /// Minimum chemotaxis bias when cell experiences a chem. concentration = 0.
    #[builder(with = |min: f32| { if min > 1. { panic!("`min` must be between 0 and 1") } else { min } } )]
    pub chemotaxis_min: f32,
    pub adhesion: StaticAdhesion,
    max_chem: u32
}

impl MyPotts {

    /// This is the chemotaxis bias as defined in Colizzi, 2020.
    ///
    /// The chemotaxis bias for the contact inhibition model is built-in the Act bias.
    fn chemotaxis_bias(&self, pos_source: Pos<usize>, pos_target: Pos<usize>, env: &MyEnvironment) -> f32 {
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

    /// This includes the combined copy biases for Act and chemotaxis as per Camley, 2016.
    fn contact_biases(&self, pos_source: Pos<usize>, pos_target: Pos<usize>, env: &MyEnvironment) -> f32 {
        let act_source = self.mean_act(pos_source, env);
        let act_target = self.mean_act(pos_target, env);
        -self.act_lambda / env.act_max as f32 * (act_source - act_target)
    }

    // TODO!: confirm with Sandro that this is right (Ava's thesis doesnt agree with her code which doesnt agree
    //  with what i think is the right implementation)
    fn contact_chemotaxis(&self, cell_index: CellIndex, env: &MyEnvironment) -> f32 {
        let cell = env.cells.get_cell(cell_index);
        cell.chem_mass as f32
            / cell.area as f32
            / self.max_chem as f32
            * (1. - self.chemotaxis_min)
            + self.chemotaxis_min
    }

    /// Geometric mean of Act content of neighbourhood of `pos`, multiplied by the relative chemotaxis term.
    fn mean_act(&self, pos: Pos<usize>, env: &MyEnvironment) -> f32 {
        let cell_spin = env.cell_lattice[pos];
        let Spin::Some(cell_index) = cell_spin else {
            return 0.;
        };

        // We do precompute these positions for pos_target in `attempt_site_copy`
        // Reusing that computation could be slightly faster
        // TODO!: Turns out this is a decent performance gain
        let (count, product) = env
            .valid_neighbours(pos)
            .filter(|&pos| env.cell_lattice[pos] == cell_spin)
            .fold(
                // Use f64 throughout the calculation to prevent overflow
                // We can alternatively sum logs instead of calculating the product
                // This is a bit more costly though
                (1, env.act_lattice[pos] as f64),
                |(count, product), pos2| {
                    let act = env.act_lattice[pos2];
                    (count + 1, product * act as f64)
                }
            );
        product.pow(1. / count as f64) as f32 * self.contact_chemotaxis(cell_index, env)
    }

    fn perimeter_energy_diff(&self, delta_perimeter: i32, perimeter: u32, target_perimeter: u32) -> f32 {
        2. * self.perimeter_lambda * delta_perimeter as f32 * (perimeter as f32 - target_perimeter as f32) + self.perimeter_lambda
    }

    fn delta_hamiltonian_perimeter(&self, spin_source: Spin, spin_target: Spin, env: &MyEnvironment) -> f32 {
        self.delta_perimeter_energy(spin_source, env)
            + self.delta_perimeter_energy(spin_target, env)
    }

    #[inline]
    fn delta_perimeter_energy(&self, spin: Spin, env: &MyEnvironment) -> f32 {
        if let Spin::Some(cell_index) = spin {
            let cell = env.env().cells.get_cell(cell_index);
            self.perimeter_energy_diff(
                cell.delta_perimeter.expect("`delta_perimeter` not set"),
                cell.perimeter,
                cell.target_perimeter
            )
        } else {
            0.
        }
    }
}

impl Potts for MyPotts {
    type Environment = MyEnvironment;

    fn boltz_t(&self) -> f32 {
        self.boltz_t
    }

    fn size_lambda(&self) -> f32 {
        self.size_lambda
    }

    fn copy_biases(&self, pos_source: Pos<usize>, pos_target: Pos<usize>, env: &Self::Environment) -> f32 {
        if self.enable_migration {
            self.contact_biases(pos_source, pos_target, env)
        } else {
            0.
        }
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

        // TODO!: Benchmark vs revalidating the positions inside update_delta_perimeter
        let neighs_target = env.env().valid_neighbours(pos_target).map(|neigh| {
            env.env().cell_lattice[neigh]
        }).collect::<Vec<_>>();

        if let Spin::Some(cell_index) = spin_source {
            env.update_delta_perimeter(true, cell_index, neighs_target.iter().copied());
        }

        if let Spin::Some(cell_index) = spin_target {
            env.update_delta_perimeter(false, cell_index, neighs_target.iter().copied());
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
