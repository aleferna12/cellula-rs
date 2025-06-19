use std::env::args;
use std::hash::{DefaultHasher, Hash, Hasher};
use image::RgbImage;
use log;
use crate::environment::Environment;
use crate::parameters::Parameters;

// TODO: license, hello etc
/// Welcomes the user and spits out information about the model parameters.
///
/// This should not require the model to be correctly initiated.
pub fn welcome(parameters: &Parameters) {
    let command = args()
        .collect::<Vec<_>>()
        .join(" ");
    log::info!("Command used: {}", command);
    log::info!("Model parameters:");
    log::info!("{:?}", parameters);
}

pub fn simulation_frame(env: &Environment) -> RgbImage {
    fn hash_sigma(sigma: i16) -> [u8; 3] {
        if sigma == 0 {
            return [255, 255, 255]
        }

        let mut hasher = DefaultHasher::new();
        sigma.hash(&mut hasher);
        let hashed = hasher.finish();
        [
            (hashed & 0xFF).try_into().unwrap(),
            (hashed >> 8 & 0xFF) as u8,
            (hashed >> 16 & 0xFF) as u8,
        ]
    }
    
    let sigmas: Vec<_> = env.cell_lattice
        .iter_values()
        .flat_map(|val| { hash_sigma(val) })
        .collect();

    RgbImage::from_vec(
        env.width() as u32,
        env.height() as u32,
        sigmas
    ).unwrap()
}