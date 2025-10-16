use crate::chem_environment::ChemEnvironment;
use crate::clonal_potts::ClonalPotts;
use bon::Builder;
use cellulars_lib::potts::Potts;
use cellulars_lib::step::Step;
use rand_xoshiro::Xoshiro256StarStar;

// TODO: this struct can be made general if CellularAutomata is also general
#[derive(Clone, Builder)]
pub struct Pond {
    pub env: ChemEnvironment,
    pub potts: ClonalPotts,
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
            self.env.cells.iter_mut().for_each(|cell| cell.update());
            if self.division_enabled {
                self.env.reproduce();
            }
            self.env.feed_cells();
        }
        self.time_step += 1;
    }
}