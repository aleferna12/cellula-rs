//! Contains logic associated with [`Pond`].

use crate::traits::potts_algorithm::PottsAlgorithm;
use crate::traits::step::Step;
use rand::Rng;

/// A pond is responsible for updating a [`Habitable`](crate::prelude::Habitable) environment
/// using a [`PottsAlgorithm`].
///
/// Comparisons using [`PartialEq`] and [`Eq`] do not compare [`Pond::rng`].
#[derive(Clone, Debug)]
pub struct Pond<P: PottsAlgorithm, R> {
    /// Environment containing the cells.
    pub env: P::Environment,
    /// Potts algorithm with which to update the CA.
    pub potts: P,
    /// Random number generator unique to this pond.
    pub rng: R,
    /// Current time-step of the pond.
    pub time_step: u32,
}

impl<P: PottsAlgorithm, R> Pond<P, R> {
    /// Makes a new pond from its constituent parts.
    pub fn new(env: P::Environment, potts: P, rng: R, time_step: u32) -> Self {
        Self {
            env,
            potts,
            rng,
            time_step
        }
    }
}

impl<P: PottsAlgorithm, R: Rng> Step for Pond<P, R> {
    fn step(&mut self) {
        self.potts.step(&mut self.env, &mut self.rng);
        self.time_step += 1;
    }
}

impl<P, R> PartialEq for Pond<P, R>
where
    P: PottsAlgorithm + PartialEq,
    P::Environment: PartialEq {
    // Dont compare rng state
    fn eq(&self, other: &Self) -> bool {
        self.time_step == other.time_step
            && self.potts == other.potts
            && self.env == other.env
    }
}

impl<P, R> Eq for Pond<P, R> where P: PottsAlgorithm + Eq, P::Environment: Eq {}
