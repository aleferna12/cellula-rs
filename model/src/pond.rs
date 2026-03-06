//! Contains logic required to run an instance of a simulation in a [`Pond`].

use bon::bon;
use crate::my_environment::MyEnvironment;
use crate::potts::Potts;
use cellulars::traits::step::Step;
use rand_xoshiro::Xoshiro256StarStar;
use cellulars::prelude::PottsAlgorithm;

/// [`Pond`] is responsible for updating an [`MyEnvironment`] using the [`Potts`] algorithm.
///
/// All simulation logic is contained here, while [`Model`](crate::model::Model) is responsible for IO.
#[derive(Clone)]
pub struct Pond {
    pub env: MyEnvironment,

    pub potts: Potts,

    pub rng: Xoshiro256StarStar,

    /// Period with which the cells' [`MyCell::update()`](crate::my_cell::MyCell::update()) method should be called.
    pub update_period: u32,

    /// Whether cell division is enabled.
    pub division_enabled: bool,

    time_step: u32
}

#[bon]
impl Pond {
    /// Makes a new [`Pond`] from an existing [`Pond`].
    #[builder]
    pub fn new(
        env: MyEnvironment,
        potts: Potts,
        rng: Xoshiro256StarStar,
        update_period: u32,
        division_enabled: bool,
        time_step: u32
    ) -> Self {
        Self {
            env,
            potts,
            rng,
            update_period,
            division_enabled,
            time_step
        }
    }

    /// Removes all cells from the pond and returns it to a clean state.
    pub fn wipe_out(&mut self) {
        self.env.wipe_out();
    }

    /// Returns the current time-step of the pond.
    ///
    /// Updated by [`Pond::step()`].
    pub fn time_step(&self) -> u32 {
        self.time_step
    }
}

impl Step for Pond {
    fn step(&mut self) {
        if self.time_step.is_multiple_of(self.update_period) {
            self.env.env.cells
                .iter_mut()
                .for_each(|rel_cell| rel_cell.cell.update());
            if self.division_enabled {
                self.env.reproduce();
            }
        }
        self.potts.step(&mut self.env, &mut self.rng);
        self.time_step += 1;
    }
}