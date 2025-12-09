use crate::cell::FitCell;
use crate::evolution::genome::Genome;
use crate::evolution::selector::{Selector, WeightedSelection};
use crate::my_environment::MyEnvironment;
use crate::my_potts::MyPotts;
use bon::Builder;
use cellulars_lib::basic_cell::Cellular;
use cellulars_lib::potts::Potts;
use cellulars_lib::step::Step;
use indexmap::IndexMap;
use rand::Rng;
use rand_xoshiro::Xoshiro256StarStar;

#[derive(Clone, Builder)]
pub struct Pond {
    pub env: MyEnvironment,
    pub potts: MyPotts,
    pub rng: Xoshiro256StarStar,
    pub cell_target_area: u32,
    pub enable_division: bool,
    pub season_duration: u32,
    pub half_fitness: f32,
    pub reproduction_steps: u32,
    #[builder(default = 0)]
    pub time_step: u32
}

impl Pond {
    pub fn wipe_out(&mut self) {
        self.env.wipe_out();
    }

    // With some unsafe code we can return Vec<&RelCell> from this function, but it would
    // require that self.divide_cell never invalidates any references to self.cells
    // we need thorough testing of self.divide_cells to make this change, and the performance
    // gain is minimal (although the ergonomic gains are significant)
    pub fn reproduce(&mut self) {
        let pop_size = self.env.cells.n_valid();
        let fit_cells = self.env.cells.iter().filter_map(|cell| {
            if !cell.is_valid() {
                return None;
            }
            Some(FitCell {
                cell,
                half_fit: self.half_fitness
            })
        }).collect::<Vec<_>>();

        let mut selector = WeightedSelection {
            select_n: pop_size,
            rng: &mut self.rng
        };
        let mut divide = IndexMap::new();
        for fit in selector.select(&fit_cells) {
            let entry = divide.entry(fit.cell.index).or_insert(0u32);
            *entry += 1;
        }

        let mut divisions_left = true;
        while divisions_left {
            divisions_left = false;
            for (cell_index, divide_n) in &mut divide {
                if !self.env.can_add_cell() {
                    return;
                }
                // Could instead remove the entries from the index map after divide_n == 0
                // but would have to store entries to remove in a vec
                if *divide_n == 0 {
                    continue;
                }

                let mom = self
                    .env
                    .cells
                    .get_cell(*cell_index);
                let new_cell = self.env.divide_cell(mom.index);
                if new_cell.is_valid() {
                    let new_index = new_cell.index;
                    // We could also instead choose to mutate at a fix rate throughout the cell's life cycle
                    self.env.cells.get_cell_mut(new_index).genome.attempt_mutate(&mut self.rng);
                }

                *divide_n -= 1;
                if *divide_n > 0 {
                    divisions_left = true;
                }
            }

            for _ in 0..self.reproduction_steps {
                self.potts.step(&mut self.env, &mut self.rng);
            }
        }
    }

    pub fn erase_half_population(&mut self) {
        let mut valid_cells = self.env.cells.iter().filter_map(|cell| {
            if !cell.is_valid() {
                return None;
            }
            Some(cell.index)
        }).collect::<Vec<_>>();
        let population_n = valid_cells.len();
        let mut kill_count = 0;
        while kill_count < population_n / 2 {
            let vec_index = self.rng.random_range(0..valid_cells.len());
            let cell_index = valid_cells[vec_index];
            self.env.erase_cell(cell_index);
            valid_cells.swap_remove(vec_index);
            kill_count += 1;
        }
    }
}

impl Step for Pond {
    fn step(&mut self) {
        self.potts.step(&mut self.env, &mut self.rng);
        if self.time_step % self.season_duration == 0 {
            if self.enable_division {
                // TODO!: is this what sandro did in his paper?
                //  its different from reproducing only cells selected multiple times and killing those missing
                self.reproduce();
                self.erase_half_population();
            }
            self.env.make_next_chem_gradient(&mut self.rng);
        }
        for val in self.env.act_lattice.iter_values_mut() {
            if *val > 0 {
                *val -= 1;
            }
        }
        self.time_step += 1;
    }
}