use crate::my_environment::MyEnvironment;
use crate::my_potts::MyPotts;
use crate::evolution::selector::Fit;
use bon::Builder;
use rand::Rng;
use cellulars_lib::basic_cell::Cellular;
use cellulars_lib::potts::Potts;
use cellulars_lib::step::Step;
use rand_xoshiro::Xoshiro256StarStar;
use cellulars_lib::positional::pos::Pos;

#[derive(Clone, Builder)]
pub struct Pond {
    pub env: MyEnvironment,
    pub potts: MyPotts,
    pub rng: Xoshiro256StarStar,
    pub update_period: u32,
    pub cell_target_area: u32,
    pub enable_division: bool,
    pub enable_cell_updates: bool,
    pub season_duration: u32,
    #[builder(default = [
        (0, 0).into(), (env.width() - 1, 0).into(),
        (0, env.height() - 1).into(),
        (env.width() - 1, env.height() - 1).into()
    ])]
    corners: [Pos<usize>; 4],
    #[builder(default = 0)]
    next_corner: usize,
    #[builder(default = 0)]
    pub time_step: u32
}

impl Pond {
    pub fn wipe_out(&mut self) {
        self.env.wipe_out();
    }

    pub fn make_next_chem_gradient(&mut self) -> Pos<usize> {
        let curr_corner = self.next_corner;
        self.env.make_chem_gradient(self.corners[curr_corner]);
        while self.next_corner == curr_corner {
            self.next_corner = self.rng.random_range(0..self.corners.len());
        }
       self.corners[curr_corner]
    }
}

impl Step for Pond {
    fn step(&mut self) {
        self.potts.step(&mut self.env, &mut self.rng);
        if self.time_step % self.update_period == 0 {
            if self.enable_cell_updates {
                self.env.cells.iter_mut().for_each(|cell| cell.update());
            }
            if self.enable_division {
                self.env.reproduce(&mut self.rng);
            }
        }
        for val in self.env.act_lattice.iter_values_mut() {
            if *val > 0 {
                *val -= 1;
            }
        }
        self.time_step += 1;
    }
}

impl Fit for Pond {
    fn fitness(&self) -> f32 {
        let tot_fit: f32 = self
            .env
            .cells
            .iter()
            .filter(|cell| cell.is_valid())
            .map(|c| { c.fitness() })
            .sum();
        tot_fit / self.env.cells.n_valid() as f32
    }
}