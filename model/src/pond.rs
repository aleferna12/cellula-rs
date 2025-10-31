use crate::evolution::selector::Fit;
use crate::my_environment::MyEnvironment;
use crate::my_potts::MyPotts;
use bon::Builder;
use cellulars_lib::basic_cell::Cellular;
use cellulars_lib::potts::Potts;
use cellulars_lib::step::Step;
use rand_xoshiro::Xoshiro256StarStar;
use crate::cell::Cell;

// TODO: this struct can be made general if CellularAutomata is also general
#[derive(Clone, Builder)]
pub struct Pond {
    pub env: MyEnvironment,
    pub potts: MyPotts,
    pub rng: Xoshiro256StarStar,
    pub update_period: u32,
    pub cell_target_area: u32,
    pub division_enabled: bool,
    #[builder(default = 0)]
    pub time_step: u32,
}

impl Pond {
    pub fn wipe_out(&mut self) {
        self.env.wipe_out();
    }
}

impl Step for Pond {
    fn step(&mut self) {
        self.potts.step(&mut self.env, &mut self.rng);
        if self.time_step % self.update_period == 0 {
            self.env.cells.iter_mut().for_each(|cell| {
                if let Cell::Amoeba(amoeba) = &mut cell.cell {
                    // TODO: do something with this, cell must express a target vector of proteins
                    //  right before the end of a season or something like that
                    amoeba.genome.input_signals[0] = self.time_step as f32;
                    amoeba.update();
                }
            });
            if self.division_enabled {
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