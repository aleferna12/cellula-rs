use std::error::Error;
use std::path::Path;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use crate::adhesion::{ClonalAdhesion};
use crate::ca::CA;
use crate::environment::Environment;
use crate::io::{create_directories, simulation_image, IMAGES_PATH, CONFIG_COPY_PATH, MovieMaker};
use crate::parameters::Parameters;

pub struct Model {
    pub env: Environment,
    pub ca: CA<ClonalAdhesion>,
    pub rng: Xoshiro256StarStar,
    movie_maker: Option<MovieMaker>,
    parameters: Parameters
}

impl Model {
    // Prevent from mutating, since values might have been used to set state already
    pub fn parameters(&self) -> &Parameters {
        &self.parameters
    }
    
    pub fn run(&mut self, steps: u32) {
        log::info!("Starting simulation");

        let mut issued_image_warning = false;
        let mut issued_movie_warning = false;
        for i in 0..=steps {
            if i % self.parameters.io.image_period == 0 {
                let saved = simulation_image(&self.env).save(
                    Path::new(&self.parameters.io.outdir)
                        .join(IMAGES_PATH)
                        .join(format!("{i}.{}", &self.parameters.io.image_format.to_lowercase()))
                );
                if let Err(e) = saved {
                    if !issued_image_warning {
                        log::warn!("Failed to save simulation frame at time step {} with error `{}`", i, e);
                        issued_image_warning = true;
                    }
                }
            }

            if let Some(mm) = &mut self.movie_maker {
                if i % self.parameters.movie.frame_period == 0 && mm.window_works() {
                    let resized = image::imageops::resize(
                        &simulation_image(&self.env), 
                        self.parameters.movie.width,
                        self.parameters.movie.height, 
                        image::imageops::Nearest
                    );
                    if let Err(e) = mm.update(&resized) {
                        if !issued_movie_warning {
                            log::warn!("Failed to display simulation frame at time step {} with error `{}`", i, e);
                            issued_movie_warning = true;
                        }
                    }
                }
            }

            self.ca.step(&mut self.env, &mut self.rng);
            if self.env.time_to_update(i) {
                self.env.cells.update_cells();
                self.env.reproduce();
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
        let mut model = Self {
            env: Environment::new(
                parameters.environment.clone(),
                &mut rng
            ),
            ca: parameters.cellular_automata.clone().into(),
            rng,
            // Initialised in setup
            movie_maker: None,
            parameters
        };

        log::info!("Setting model up");
        log::info!("Creating output directory");
        create_directories(&model.parameters.io.outdir, model.parameters.io.replace_outdir)?;

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

        // Hopefully this prevents most compatibility problems
        // A more extreme solution is to make minifb an optional dependency
        if model.parameters.movie.show {
            log::info!("Creating window for real-time movie display");
            match MovieMaker::new(model.parameters.movie.width, model.parameters.movie.height) {
                Ok(maker) => model.movie_maker = Some(maker),
                Err(e) => log::warn!("Failed to initialise movie maker with error `{}`", e),
            }
        };
        
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