use crate::clonal_potts::ClonalPotts;
use crate::chem_environment::ChemEnvironment;
use bon::Builder;
use cellulars_lib::basic_cell::Cellular;
use cellulars_lib::evolution::selector::Fit;
use rand_xoshiro::Xoshiro256StarStar;
use cellulars_lib::potts::Potts;

// TODO: this struct can be made general if CellularAutomata is also general
#[derive(Clone, Builder)]
pub struct Pond {
    pub env: ChemEnvironment,
    pub potts: ClonalPotts,
    pub rng: Xoshiro256StarStar,
    pub update_period: u32,
    pub cell_target_area: u32,
    pub division_enabled: bool,
    pub cell_search_scaler: f32,
    #[builder(default = 0)]
    pub time_step: u32,
}

impl Pond {
    pub fn step(&mut self) {
        self.potts.step(&mut self.env, &mut self.rng);
        if self.time_step % self.update_period == 0 {
            self.env.cells.iter_mut().for_each(|cell| cell.update());
            if self.division_enabled {
                self.env.reproduce(self.cell_search_scaler, &mut self.rng);
            }
        }
        self.time_step += 1;
    }

    pub fn wipe_out(&mut self) {
        self.env.wipe_out();
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