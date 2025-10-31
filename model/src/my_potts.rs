use crate::my_environment::MyEnvironment;
use bon::Builder;
use cellulars_lib::adhesion::{AdhesionSystem, StaticAdhesion};
use cellulars_lib::potts::Potts;
use cellulars_lib::spin::Spin;

// This could be a module but it's convenient to be able to access the relevant parameters
// Also we might eventually want to implement multiple CA choices, in which case I can "easily" make CA a trait 
// that just implements `step()`
#[derive(Clone, Builder)]
pub struct MyPotts {
    pub boltz_t: f32,
    pub size_lambda: f32,
    pub enable_migration: bool,
    pub adhesion: StaticAdhesion
}

impl Potts for MyPotts {
    type Environment = MyEnvironment;

    fn boltz_t(&self) -> f32 {
        self.boltz_t
    }

    fn size_lambda(&self) -> f32 {
        self.size_lambda
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