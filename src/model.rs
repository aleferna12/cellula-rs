use std::error::Error;
use std::path::Path;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use crate::adhesion::{ClonalAdhesion};
use crate::ca::CA;
use crate::environment::Environment;
use crate::io::{CONFIG_COPY_PATH, IoManager};
use crate::parameters::Parameters;

pub struct Model {
    pub env: Environment,
    pub ca: CA<ClonalAdhesion>,
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
            self.io_manager.image_io(time_step, &self.env, &self.ca.adhesion.clone_pairs);
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
            ca: parameters.cellular_automata.clone().into(),
            io_manager: IoManager::try_from(parameters.io.clone())?,
            rng,
            parameters
        };

        log::info!("Setting model up");
        log::info!("Creating output directory");
        model.io_manager.create_directories()?;

        let params_copy = Path::new(&model.parameters.io.outdir).join(CONFIG_COPY_PATH);
        log::info!("Saving copy of parameters to `{}`", &params_copy.display());
        std::fs::write(
            params_copy,
            format!(
                "{}\n{}",
                "# This is a copy of the parameters used in the simulation",
                toml::to_string(&model.parameters)?
            )
        )?;

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