use crate::cellular_automata::CellularAutomata;
use crate::chem_environment::ChemEnvironment;
use crate::clonal_adhesion::ClonalAdhesion;
use crate::genetics::genome::Genome;
use bon::Builder;
use cellulars_lib::basic_cell::{Alive, Cellular};
use cellulars_lib::environment::Habitable;
use cellulars_lib::evolution::selector::Fit;
use cellulars_lib::lattice_entity::LatticeEntity;
use cellulars_lib::positional::boundaries::Boundary;
use rand::Rng;
use rand_xoshiro::Xoshiro256StarStar;

// TODO: this struct can be made general if CellularAutomata is also general
#[derive(Clone, Builder)]
pub struct Pond {
    pub env: ChemEnvironment,
    pub ca: CellularAutomata<ClonalAdhesion>,
    pub rng: Xoshiro256StarStar,
    pub update_period: u32,
    pub cell_target_area: u32,
    pub division_enabled: bool,
    pub cell_search_scaler: f32,
    #[builder(default = 0)]
    pub(crate) time_step: u32,
}

impl Pond {
    pub fn step(&mut self) {
        self.ca.step(&mut self.env, &mut self.rng);
        if self.time_step % self.update_period == 0 {
            self.env.cells.iter_mut().for_each(|cell| cell.update());
            if self.division_enabled {
                self.reproduce();
            }
        }
        self.time_step += 1;
    }
    
    pub fn reproduce(&mut self) {
        let new_spins = self.env.reproduce(self.cell_search_scaler);
        for spin in new_spins {
            self.ca.adhesion.update_clones(spin, &self.env);
            // We could also instead choose to mutate at a fix rate throughout the cell's life cycle
            if let LatticeEntity::SomeCell(cell) = self.env.cells.get_entity_mut(spin) {
                cell.genome.attempt_mutate(&mut self.rng);
            } else {
                panic!("newborn is not a cell")
            }
        }
    }

    pub fn wipe_out(&mut self) {
        self.env.wipe_out();
        self.ca.adhesion.clones_table.clear();
    }
}

impl Fit for Pond {
    fn fitness(&self) -> f32 {
        let tot_fit: f32 = self
            .env
            .cells()
            .iter()
            .filter(|cell| cell.is_valid())
            .map(|c| { c.fitness() })
            .sum();
        tot_fit / self.env.cells().n_valid() as f32
    }
}