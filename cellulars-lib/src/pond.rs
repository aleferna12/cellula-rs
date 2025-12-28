//! Contains logic associated with [Pond].

use crate::traits::potts_algorithm::PottsAlgorithm;
use crate::traits::step::Step;
use rand::Rng;

/// A pond is responsible for updating a [Habitable](crate::habitable::Habitable) using a [Potts] algorithm.
#[derive(Clone)]
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