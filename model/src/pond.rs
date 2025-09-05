use crate::cell::Cell;
use crate::cellular_automata::CellularAutomata;
use crate::chem_space::ChemEnvironment;
use crate::genetics::genome::Genome;
use cellulars_lib::adhesion::ClonalAdhesion;
use cellulars_lib::cellular::{Cellular, RelCell};
use cellulars_lib::constants::Spin;
use cellulars_lib::evolution::selector::Fit;
use cellulars_lib::lattice_entity::LatticeEntity;
use cellulars_lib::lattice_entity::LatticeEntity::Medium;
use cellulars_lib::positional::boundary::Boundary;
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::positional::rect::Rect;
use rand::Rng;
use rand_xoshiro::Xoshiro256StarStar;

// TODO: this struct can be made general if CellularAutomata is also general
pub struct Pond {
    pub env: ChemEnvironment,
    pub ca: CellularAutomata<ClonalAdhesion>,
    pub rng: Xoshiro256StarStar,
    pub update_period: u32,
    pub cell_search_radius: f32,
    pub cell_target_area: u32,
    pub division_enabled: bool,
    pub max_cells: u32,
    population_exploded: bool,
    time_step: u32,
}

impl Pond {
    pub fn new(
        env: ChemEnvironment,
        ca: CellularAutomata<ClonalAdhesion>,
        rng: Xoshiro256StarStar,
        update_period: u32,
        cell_search_radius: f32,
        cell_target_area: u32,
        division_enabled: bool,
        max_cells: u32
    ) -> Self {
        Self {
            env,
            ca,
            rng,
            update_period,
            cell_search_radius,
            cell_target_area,
            division_enabled,
            max_cells,
            population_exploded: false,
            time_step: 0
        }
    }
    
    pub fn step(&mut self) {
        self.ca.step(&mut self.env, &mut self.rng);
        if self.time_step % self.update_period == 0 {
            self.env.cells.iter_mut().for_each(|cell| cell.update());
            let new_spins = self.reproduce();
            for spin in new_spins {
                self.ca.adhesion.update_clones(spin, self.cell_search_radius, &self.env);
                // We could also instead choose to mutate at a fix rate throughout the cell's life cycle
                if let LatticeEntity::SomeCell(cell) = self.env.cells.get_entity_mut(spin) {
                    cell.genome.attempt_mutate(&mut self.rng);
                } else { 
                    panic!("newborn is not a cell")
                }
            }
        }
        self.time_step += 1;
    }

    

    // With some unsafe code we can return Vec<&RelCell> from this function, but it would
    // require that self.divide_cell never invalidates any references to self.cells
    // we need thorough testing of self.divide_cells to make this change, and the performance
    // gain is minimal (although the ergonomic gains are significant)
    pub fn reproduce(&mut self) -> Vec<Spin> {
        let mut divide = vec![];
        if !self.division_enabled {
            return divide;
        }
        
        for cell in self.env.cells.iter() {
            if cell.is_dying() {
                continue;
            }
            // Currently cells don't need to express the dividing type to divide, they just need to be big enough
            if cell.area() >= cell.divide_area() {
                divide.push(cell.spin);
            }
        }
        divide.into_iter().filter_map(|spin| {
            if !self.can_add_cell() {
                return None;
            }
            match self.divide_cell(spin) {
                Err(e) => {
                    log::warn!("Failed to divide cell with spin {spin} with error `{e:?}`");
                    None
                },
                Ok(cell) => Some(cell.spin)
            }
        }).collect()
    }
}

impl Fit for Pond {
    fn fitness(&self) -> f32 {
        let tot_fit: f32 = self
            .env
            .cells
            .iter()
            .filter(|cell| cell.is_alive())
            .map(|c| { c.fitness() })
            .sum();
        tot_fit / self.env.cells.n_alive() as f32
    }
}

#[derive(Debug)]
pub enum DivisionError {
    NewCellTooSmall,
    NewCellTooBig
}