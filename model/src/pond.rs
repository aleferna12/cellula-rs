//! Contains logic required to run an instance of a simulation in a [`Pond`].

use crate::potts::Potts;
use cellulars_lib::base::base_pond::BasePond;
use cellulars_lib::traits::step::Step;
use rand_xoshiro::Xoshiro256StarStar;
use crate::environment::Environment;

/// A pond is responsible for updating an [`Environment`](crate::environment::Environment) using the [`Potts`] algorithm.
///
/// All simulation logic is contained here, while [`Model`](crate::model::Model) is responsible for IO.
#[derive(Clone)]
pub struct Pond {
    /// Inner [`BasePond`].
    pub base_pond: BasePond<Potts, Xoshiro256StarStar>,
    /// Period with which the cells' [`Cell::update()`](crate::cell::Cell::update()) method should be called.
    pub update_period: u32,
    /// Whether cell division is enabled.
    pub division_enabled: bool
}

impl Pond {
    /// Makes a new [`Pond`] from an existing [`BasePond`].
    pub fn new(
        pond: BasePond<Potts, Xoshiro256StarStar>,
        update_period: u32,
        division_enabled: bool
    ) -> Self {
        Self {
            base_pond: pond,
            update_period,
            division_enabled
        }
    }

    /// Returns a reference to the pond's inner [`Environment`].
    pub fn env(&self) -> &Environment {
        &self.base_pond.env
    }

    /// Returns a mutable reference to the pond's inner [`Environment`].
    pub fn env_mut(&mut self) -> &mut Environment {
        &mut self.base_pond.env
    }

    /// Removes all cells from the pond and returns it to a clean state.
    pub fn wipe_out(&mut self) {
        self.env_mut().wipe_out();
    }

    /// Returns the current time-step of the pond.
    ///
    /// Updated by [`Pond::step()`].
    pub fn time_step(&self) -> u32 {
        self.base_pond.time_step
    }
}

impl Step for Pond {
    fn step(&mut self) {
        if self.base_pond.time_step.is_multiple_of(self.update_period) {
            self.env_mut().base_env.cells
                .iter_mut()
                .for_each(|rel_cell| rel_cell.cell.update());
            if self.division_enabled {
                self.env_mut().reproduce();
            }
        }
        self.base_pond.step();
    }
}