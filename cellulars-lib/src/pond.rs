//! Contains logic associated to [Pond].

use crate::potts::Potts;
use crate::step::Step;
use rand::Rng;

/// A pond is responsible for updating a [Habitable](crate::habitable::Habitable) using a [Potts] algorithm.
#[derive(Clone)]
pub struct Pond<P: Potts, R> {
    /// Environment containing the cells.
    pub env: P::Environment,
    /// Potts algorithm with which to update the CA.
    pub potts: P,
    /// Random number generator unique to this pond.
    pub rng: R,
    /// Current time-step of the pond.
    pub time_step: u32,
}

impl<P: Potts, R> Pond<P, R> {
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

impl<P: Potts, R: Rng> Step for Pond<P, R> {
    fn step(&mut self) {
        self.potts.step(&mut self.env, &mut self.rng);
        self.time_step += 1;
    }
}