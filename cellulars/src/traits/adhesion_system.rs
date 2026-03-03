//! Contains logic associated with [`AdhesionSystem`].

use crate::constants::FloatType;
use crate::spin::Spin;

/// Trait defining a way to calculate adhesion at the interface between two spins.
pub trait AdhesionSystem<C = ()> {
    /// Returns the energy at the interface between `spin1` and `spin2`, given a context `C`.
    fn adhesion_energy(&self, spin1: Spin, spin2: Spin, context: &C) -> FloatType;
}