use crate::adhesion::ClonalAdhesion;
use crate::cell::{Cell, Fit};
use crate::cellular_automata::CellularAutomata;
use crate::constants::{BoundaryType, NeighbourhoodType};
use crate::environment::{Environment, LatticeEntity};
use crate::genome::{Genome, Grn};
use rand_xoshiro::Xoshiro256StarStar;

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

#[cfg(test)]
mod tests {
    use rand::{Rng, SeedableRng};
    use rand_xoshiro::Xoshiro256StarStar;

    #[test]
    fn test_seed() {
        let mut rng = Xoshiro256StarStar::seed_from_u64(1241254152);
        let s = (0..50)
            .map(|_| rng.random_range(0..9).to_string())
            .collect::<Vec<_>>()
            .join("");
        let res = "15515320360704325727185856564110164830043067488704";
        assert_eq!(res, s);
    }
}