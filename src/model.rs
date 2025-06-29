use std::error::Error;
use std::path::Path;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use crate::ca::CA;
use crate::environment::Environment;
use crate::io::{create_directories, simulation_frame, IMAGES_PATH, CONFIG_COPY_PATH};
use crate::parameters::Parameters;
use crate::pos::Rect;

pub struct Model {
    pub env: Environment,
    pub ca: CA,
    pub rng: Xoshiro256StarStar,
    parameters: Parameters
}

impl Model {
    pub fn new(parameters: Parameters) -> Self {
         Self {
             env: Environment::new(
                parameters.env.width,
                parameters.env.height,
                parameters.env.neigh_r,
                parameters.env.enclose
             ),
             ca: CA::new(
                 parameters.ca.boltz_t,
                 parameters.ca.size_lambda,
                 parameters.ca.cell_energy,
                 parameters.ca.medium_energy,
                 parameters.ca.solid_energy
             ), 
             rng: if parameters.general.seed == 0 { 
                 Xoshiro256StarStar::from_os_rng()
             } else {
                 Xoshiro256StarStar::seed_from_u64(parameters.general.seed) 
             },
             parameters 
         }
    }
    
    // Prevent from mutating, since values might have been used to set state already
    pub fn parameters(&self) -> &Parameters {
        &self.parameters
    }
    
    pub fn setup(&mut self) -> Result<(), Box<dyn Error>> {
        log::info!("Setting model up");
        log::info!("Creating output directory");
        create_directories(&self.parameters.io.outdir, self.parameters.io.replace_outdir)?;

        let params_copy = Path::new(&self.parameters.io.outdir).join(CONFIG_COPY_PATH);
        log::info!("Saving copy of parameters to `{}`", &params_copy.display());
        std::fs::write(
            params_copy, 
            format!(
                "{}\n{}", 
                "# This is a copy of the parameters used in the simulation",
                toml::to_string(&self.parameters)?
            )
        )?;
        
        let mut cell_count = 0;
        let cell_side = (self.parameters.env.cell_start_area as f32).sqrt() as usize;
        for _ in 0..self.parameters.env.n_cells {
            let pos = self.env.cell_lattice.random_pos(&mut self.rng);
            let cell = self.env.spawn_rect_cell(
                Rect::new(
                    pos,
                    (pos.x + cell_side, pos.y + cell_side).into()
                ),
                self.parameters.ca.cell_target_area
            );
            if cell.is_some() {
                cell_count += 1;
            }
        }
        log::info!("Created {} out of the {} cells requested", cell_count, self.parameters.env.n_cells);
        
        Ok(())
    }
    
    pub fn run(&mut self, steps: u32) {
        log::info!("Starting simulation");
        for i in 0..=steps {
            if i % self.parameters.io.image_period == 0 {
                let saved_image = simulation_frame(&self.env)
                    .save(Path::new(&self.parameters.io.outdir)
                    .join(IMAGES_PATH)
                    .join(format!("{i}.{}", &self.parameters.io.image_format.to_lowercase())));
                if let Err(e) = saved_image {
                    log::warn!("Failed to save simulation frame at time step {} with error `{}`", i, e);
                }
            }
            self.step();
        }
    }
    
    pub fn step(&mut self) {
        self.ca.step(&mut self.env, &mut self.rng);
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