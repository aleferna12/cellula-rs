use rand::Rng;
use crate::adhesion::ClonalAdhesion;
use crate::cell::{Cell, Cellular, Fit, RelCell};
use crate::cellular_automata::CellularAutomata;
use crate::constants::{BoundaryType, NeighbourhoodType};
use crate::environment::{Environment, LatticeEntity};
use crate::genome::{Genome, Grn};
use rand_xoshiro::Xoshiro256StarStar;

// TODO: this struct can be made general if CellularAutomata is also general
pub struct Pond {
    pub env: Environment<Cell<Grn<1, 1>>, NeighbourhoodType, BoundaryType>,
    pub ca: CellularAutomata<ClonalAdhesion>,
    rng: Xoshiro256StarStar,
    time_step: u32
}

impl Pond {
    pub fn new(
        env: Environment<Cell<Grn<1, 1>>, NeighbourhoodType, BoundaryType>,
        ca: CellularAutomata<ClonalAdhesion>,
        rng: Xoshiro256StarStar,
    ) -> Self {
        Self {
            env,
            ca,
            rng,
            time_step: 0
        }
    }
    
    pub fn step(&mut self) {
        self.ca.step(&mut self.env, &mut self.rng);
        if self.env.time_to_update(self.time_step) {
            self.env.cells.update_cells();
            let new_spins = self.env.reproduce();
            for spin in new_spins {
                // TODO!: This function is preventing CA to be generalised in Pond
                self.ca.adhesion.update_clones(spin, &self.env);
                // We could also instead choose to mutate at a fix rate throughout the cell's life cycle
                if let LatticeEntity::SomeCell(cell) = self.env.cells.get_entity_mut(spin) {
                    cell.genome.attempt_mutate(&mut self.rng);
                } else { 
                    panic!("Newborn is not a cell")
                }
            }
        }
        self.time_step += 1;
    }

    // TODO!: i think that this function being here is an indicator that all functions that dynamically change
    //  the environment should also be in Pond
    //  What if in other impls of CA I need to do some operation on CA whenever a cell is spawned?
    //  This way Pond becomes "an environment where change happens through a CA"
    pub fn kill_cell(&mut self, cell: &mut RelCell<impl Cellular>) {
        for pos in self.env.space.search_cell_box(cell, self.env.cell_search_radius) {
            // TODO!: Parameterize chance of medium
            if self.rng.random::<f32>() < 0.1 {
                self.env.space.cell_lattice[pos] = LatticeEntity::Medium.discriminant();
            }
        }
        for i in 0..self.ca.adhesion.clone_pairs.length() {
            self.ca.adhesion.clone_pairs[(cell.spin as usize, i)] = false
        }
        cell.die();
    }
    
    pub fn wipe_out(&mut self) {
        self.env.cells.wipe_out();
        self.env.space.cell_lattice.clear();
        self.ca.adhesion.clone_pairs.clear();
    }
}

impl Fit for Pond {
    fn fitness(&self) -> f32 {
        let tot_fit: f32 = self
            .env
            .cells
            .into_iter()
            .map(|c| { c.fitness() })
            .sum();
        tot_fit / self.env.cells.n_cells() as f32
    }
}