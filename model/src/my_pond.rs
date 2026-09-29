//! Contains logic required to run an instance of a simulation in a [`MyPond`].

use crate::potts::Potts;
use cellulars_lib::base::pond::Pond;
use cellulars_lib::traits::step::Step;
use rand_xoshiro::Xoshiro256StarStar;
use crate::my_environment::MyEnvironment;

/// A pond is responsible for updating an [`Environment`](crate::my_environment::MyEnvironment) using the [`Potts`] algorithm.
///
/// All simulation logic is contained here, while [`Model`](crate::model::Model) is responsible for IO.
#[derive(Clone)]
pub struct MyPond {
    /// Inner [`Pond`].
    pub pond: Pond<Potts, Xoshiro256StarStar>,
    /// Period with which the cells' [`Cell::update()`](crate::my_cell::MyCell::update()) method should be called.
    pub update_period: u32,
    /// Whether cell division is enabled.
    pub division_enabled: bool,
    pub target_move_period: u32,
}

impl MyPond {
    /// Makes a new [`MyPond`] from an existing [`Pond`].
    pub fn new(
        pond: Pond<Potts, Xoshiro256StarStar>,
        update_period: u32,
        division_enabled: bool,
        target_move_period: u32,
    ) -> Self {
        Self {
            pond,
            update_period,
            division_enabled,
            target_move_period,
        }
    }

    /// Returns a reference to the pond's inner [`MyEnvironment`].
    pub fn env(&self) -> &MyEnvironment {
        &self.pond.env
    }

    /// Returns a mutable reference to the pond's inner [`MyEnvironment`].
    pub fn env_mut(&mut self) -> &mut MyEnvironment {
        &mut self.pond.env
    }

    /// Removes all cells from the pond and returns it to a clean state.
    pub fn wipe_out(&mut self) {
        self.env_mut().wipe_out();
    }

    /// Returns the current time-step of the pond.
    ///
    /// Updated by [`MyPond::step()`].
    pub fn time_step(&self) -> u32 {
        self.pond.time_step
    }
}

impl Step for MyPond {
    fn step(&mut self) {
        let time_step = self.pond.time_step;
        // Destructured so that the environment can be updated with the pond's own rng
        let Pond { env, rng, .. } = &mut self.pond;
        if time_step.is_multiple_of(self.update_period) {
            env.env.cells
                .iter_mut()
                .for_each(|rel_cell| rel_cell.cell.update());
            if self.division_enabled {
                env.reproduce(rng);
            }
        }
        // Colizzi 2020 turns cells towards their realized direction of motion every time-step,
        // right before the CA is updated
        env.update_persistence();
        // TODO!: parameterize
        // if self.pond.time_step.is_multiple_of(self.target_move_period) {
        //     let center = self.env().target_center.cast_as();
        //     self.env_mut().target_center = self.env().env.bounds.lattice_boundary.valid_pos(Pos::new(
        //           center.x + 1,
        //           center.y,
        //     )).unwrap().cast_as();
        //     self.env_mut().update_chem_gradient();
        // }
        self.pond.step();
    }
}