use std::env::args_os;
use crate::lattice::Lattice;
use crate::parameters::Parameters;

pub struct Model {
    pub lattice: Lattice<u32>,
    pub parameters: Parameters
}
impl Model {
    pub fn new(parameters: Parameters) -> Self {
        Self {
            lattice: Lattice::new(
                parameters.width,
                parameters.height
            ),
            parameters
        }
    }
    pub fn welcome(&self) {
        let command = args_os()
            .map(|s| s.into_string().unwrap())
            .collect::<Vec<_>>()
            .join(" ");
        println!("Model initialised with command: {}", command);
    }
}