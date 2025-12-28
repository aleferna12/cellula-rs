//! Contains logic required to run an instance of a simulation in a [MyPond].

use crate::my_potts::MyPotts;
use cellulars_lib::pond::Pond;
use cellulars_lib::traits::step::Step;
use rand_xoshiro::Xoshiro256StarStar;

/// A pond is responsible for updating a [MyEnvironment] using the [MyPotts] algorithm.
///
/// All simulation logic is contained here, while [Model](crate::model::Model) is responsible for IO.
#[derive(Clone)]
pub struct MyPond {
    /// Inner [Pond].
    pub pond: Pond<MyPotts, Xoshiro256StarStar>,
    /// Period with which the cells' [Cell::update()](crate::cell::Cell::update()) method should be called.
    pub update_period: u32,
    /// Whether cell division is enabled.
    pub division_enabled: bool
}

impl MyPond {
    /// Makes a new [MyPond] from an existing [Pond].
    pub fn new(
        pond: Pond<MyPotts, Xoshiro256StarStar>,
        update_period: u32,
        division_enabled: bool
    ) -> Self {
        Self {
            pond,
            update_period,
            division_enabled
        }
    }

    /// Removes all cells from the pond and returns it to a clean state.
    pub fn wipe_out(&mut self) {
        self.pond.env.wipe_out();
    }

    /// Returns the current time-step of the pond.
    ///
    /// Updated by [MyPond::step()].
    pub fn time_step(&self) -> u32 {
        self.pond.time_step
    }
}

impl Step for MyPond {
    fn step(&mut self) {
        if self.pond.time_step.is_multiple_of(self.update_period) {
            self.pond.env.env.cells.iter_mut().for_each(|cell| cell.update());
            if self.division_enabled {
                self.pond.env.reproduce();
            }
        }
        self.pond.step();
    }
}