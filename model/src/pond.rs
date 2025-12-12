//! Contains logic required to run an instance of a simulation in a [Pond].

use crate::my_environment::MyEnvironment;
use crate::my_potts::MyPotts;
use bon::Builder;
use cellulars_lib::potts::Potts;
use cellulars_lib::step::Step;
use rand_xoshiro::Xoshiro256StarStar;

/// A pond is responsible for updating a [MyEnvironment] using the [MyPotts] algorithm.
///
/// All simulation logic is contained here, while [Model](crate::model::Model) is responsible for IO.
#[derive(Clone, Builder)]
pub struct Pond {
    /// Environment containing the cells.
    pub env: MyEnvironment,
    /// Potts algorithm with which to update the CA.
    pub potts: MyPotts,
    /// Random number generator unique to this pond.
    pub rng: Xoshiro256StarStar,
    /// Period with which the cells' [Cell::update()](crate::cell::Cell::update()) method should be called.
    pub update_period: u32,
    /// Whether cell division is enabled.
    pub division_enabled: bool,
    /// Current time-step of the pond.
    #[builder(default = 0)]
    time_step: u32,
}

impl Pond {
    /// Removes all cells from the pond and returns it to a clean state.
    pub fn wipe_out(&mut self) {
        self.env.wipe_out();
    }

    /// Returns the current time-step of the pond.
    ///
    /// Updated by [Pond::step()].
    pub fn time_step(&self) -> u32 {
        self.time_step
    }
}

impl Step for Pond {
    fn step(&mut self) {
        self.potts.step(&mut self.env, &mut self.rng);
        if self.time_step % self.update_period == 0 {
            self.env.env_mut().cells.iter_mut().for_each(|cell| cell.update());
            if self.division_enabled {
                self.env.reproduce();
            }
        }
        self.time_step += 1;
    }
}