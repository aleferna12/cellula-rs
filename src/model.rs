use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use crate::ca::CA;
use crate::environment::Environment;
use crate::parameters::Parameters;
use crate::pos::Rect;

pub struct Model {
    pub env: Environment,
    pub ca: CA,
    pub rng: Xoshiro256StarStar,
    pub parameters: Parameters
}
impl Model {
    pub fn new(parameters: Parameters) -> Self {
         Self {
             env: Environment::new(
                parameters.width,
                parameters.height,
                parameters.neigh_r
             ),
             ca: CA::new(
                 parameters.boltz_t,
                 parameters.size_lambda
             ), 
             rng: if parameters.seed == 0 { 
                 Xoshiro256StarStar::from_os_rng()
             } else {
                 Xoshiro256StarStar::seed_from_u64(parameters.seed) 
             },
             parameters 
         }
    }
    
    pub fn setup(&mut self) {
        log::info!("Setting model up");
        self.env.spawn_rect_cell(
            Rect::new(
                (10, 10).into(),
                (20, 20).into()
            ),
            self.parameters.target_area
        );
    }
    
    pub fn run(&mut self, steps: u32) {
        log::info!("Starting simulation");
        for _ in 0..steps {
            self.step();
        }
    }
    
    pub fn step(&mut self) {
        self.ca.step(&mut self.env, &mut self.rng);
    }
}

#[cfg(test)]
mod tests {
    use clap::Parser;
    use rand::Rng;
    use crate::model::Model;
    use crate::parameters::Parameters;
    use crate::utils::TEST_SEED;

    #[test]
    fn test_xoshiro() {
        let mut model = Model::new(Parameters::parse_from(["", "--seed", &TEST_SEED.to_string()]));
        let s = (0..50)
            .map(|_| model.rng.random_range(0..9).to_string())
            .collect::<Vec<_>>()
            .join("");
        let res = "15515320360704325727185856564110164830043067488704";
        assert_eq!(res, s);
    }
    
    #[test]
    fn test_run() {
        todo!()
    }
}