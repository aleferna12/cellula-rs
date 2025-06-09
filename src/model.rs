use std::env::args_os;
use rand::SeedableRng;
use rand_xoshiro::Xoshiro256StarStar;
use crate::dish::Dish;
use crate::parameters::Parameters;

pub struct Model {
    pub dish: Dish,
    pub rng: Xoshiro256StarStar,
    pub parameters: Parameters
}
impl Model {
    pub fn new(parameters: Parameters) -> Self {
        Self {
            dish: Dish::new(
                parameters.width,
                parameters.height
            ),
            rng: if parameters.seed == 0 { 
                Xoshiro256StarStar::from_os_rng()
            } else {
                Xoshiro256StarStar::seed_from_u64(parameters.seed) 
            },
            parameters
        }
    }
    pub fn welcome(&self) {
        let command = args_os()
            .map(|s| s.into_string().unwrap())
            .collect::<Vec<_>>()
            .join(" ");
        println!("Model initialised with command: {}", command);
        println!("Model parameters:");
        println!("{:?}", self.parameters);
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
        let mut model = Model::new(Parameters::parse_from(["", "--seed", "1241254152"]));
        let s = (0..50)
            .map(|_| model.rng.random_range(0..9).to_string())
            .collect::<Vec<_>>()
            .join("");
        let res = "15515320360704325727185856564110164830043067488704";
        assert_eq!(res, s);
    }
}