//! Contains logic associated with [AdhesionSystem].

use crate::spin::Spin;

/// Trait defining a way to calculate adhesion at the interface between two spins.
pub trait AdhesionSystem {
    /// Context required to calculate the adhesion energy.
    type Context;
    /// Returns the energy at the interface between `spin1` and `spin2`, given a context.
    fn adhesion_energy(&self, spin1: Spin, spin2: Spin, context: &Self::Context) -> f32;
}