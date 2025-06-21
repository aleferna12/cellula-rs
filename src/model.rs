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
                 parameters.size_lambda,
                 parameters.cell_energy,
                 parameters.medium_energy,
                 parameters.solid_energy
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
        let mut cell_count = 0;
        let cell_side = (self.parameters.cell_start_area as f32).sqrt() as usize;
        for _ in 0..self.parameters.n_cells {
            let pos = self.env.cell_lattice.random_pos(&mut self.rng);
            let cell = self.env.spawn_rect_cell(
                Rect::new(
                    pos,
                    (pos.x + cell_side, pos.y + cell_side).into()
                ),
                self.parameters.cell_target_area
            );
            if cell.is_some() {
                cell_count += 1;
            }
        }
        log::info!("Created {} out of the {} cells requested", cell_count, self.parameters.n_cells);
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

    #[test]
    fn test_xoshiro() {
        let mut model = Model::new(Parameters::parse_from(["", "--seed", &1241254152.to_string()]));
        let s = (0..50)
            .map(|_| model.rng.random_range(0..9).to_string())
            .collect::<Vec<_>>()
            .join("");
        let res = "15515320360704325727185856564110164830043067488704";
        assert_eq!(res, s);
    }
}