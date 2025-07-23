use crate::adhesion::ClonalAdhesion;
use crate::cellular_automata::CellularAutomata;
use crate::environment::{Environment, LatticeEntity};
use crate::io::io_manager::IoManager;
use crate::io::parameters::Parameters;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use std::error::Error;
use crate::genome::Genome;

pub struct Model {
    pub env: Environment,
    pub ca: CellularAutomata<ClonalAdhesion>,
    pub rng: Xoshiro256StarStar,
    pub io_manager: IoManager,
    parameters: Parameters
}

impl Model {
    // Prevent from mutating, since values might have been used to set state already
    pub fn parameters(&self) -> &Parameters {
        &self.parameters
    }
    
    pub fn run(&mut self, steps: u32) {
        log::info!("Starting simulation");
        for time_step in 0..=steps {
            let saved = self.io_manager.image_io(
                time_step,
                &self.env, 
                &self.ca.adhesion.clone_pairs
            );
            if let Err(e) = saved {
                log::warn!("Failed to save image at time step {time_step} with error `{e}`")
            }
            self.step(time_step);
        }
    }
    
    pub fn step(&mut self, time_step: u32) {
        self.ca.step(&mut self.env, &mut self.rng);
        if self.env.time_to_update(time_step) {
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
    }
}

impl TryFrom<Parameters> for Model {
    type Error = Box<dyn Error>;

    fn try_from(parameters: Parameters) -> Result<Self, Self::Error> {
        let mut rng = if parameters.general.seed == 0 {
            Xoshiro256StarStar::from_os_rng()
        } else {
            Xoshiro256StarStar::seed_from_u64(parameters.general.seed)
        };
        let model = Self {
            env: Environment::new(
                parameters.environment.clone(),
                &mut rng
            ),
            ca: CellularAutomata::new(
                parameters.cellular_automata.clone(),
                ClonalAdhesion::new(
                    parameters.cellular_automata.adhesion.clone(),
                    parameters.environment.max_cells + LatticeEntity::first_cell_spin()
                )
            ),
            io_manager: IoManager::try_from(parameters.io.clone())?,
            rng,
            parameters
        };

        log::info!("Setting model up");
        log::info!("Creating output directory and copy of parameter file");
        model.io_manager.create_directories()?;
        model.io_manager.create_parameters_file(&model.parameters)?;

        Ok(model)
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