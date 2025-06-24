use std::env::args;
use std::fs::{create_dir_all, remove_dir_all};
use std::hash::{DefaultHasher, Hash, Hasher};
use std::path::Path;
use image::RgbImage;
use log;
use std::io;
use crate::environment::{Environment, Sigma};
use crate::parameters::Parameters;

pub(crate) static IMAGES_PATH: &str = "images";

// TODO: license, hello etc
/// Welcomes the user and spits out information about the model parameters.
///
/// This should not require the model to be correctly initialised 
/// (initialising the parameters beforehand is ergonomic with how `clap` is set up).
pub fn welcome(parameters: &Parameters) {
    let command = args()
        .collect::<Vec<_>>()
        .join(" ");
    log::info!("Command used: {}", command);
    log::info!("Model parameters:");
    log::info!("{:?}", parameters);
}

pub fn create_directories<T>(outpath: T, replace_outdir: bool) -> io::Result<()>
where
    T: AsRef<Path> {
    let outpath = outpath.as_ref(); // Convert to &Path
    let outdir_exists = outpath.try_exists()?;
    if outdir_exists {
        if replace_outdir {
            log::info!("Cleaning contents of '{}'", outpath.display());
            remove_dir_all(outpath)?;
        } else {
            return Err(io::Error::new(
                io::ErrorKind::AlreadyExists, 
                "`outdir` already exists and `replace_outdir` is `false`"
            ));
        }
    }
    create_dir_all(outpath.join(IMAGES_PATH))
}

pub fn simulation_frame(env: &Environment) -> RgbImage {
    let sigmas: Vec<_> = env.cell_lattice
        .iter_values()
        .flat_map(|val| { sigma_to_rgb(val) })
        .collect();

    RgbImage::from_vec(
        env.width() as u32,
        env.height() as u32,
        sigmas
    ).unwrap()
}

/// Converts a sigma into a unique color.
///
/// This method guarantees 5232 unique colors, starting from this sigma the colors will repeat.
fn sigma_to_rgb(sigma: Sigma) -> [u8; 3] {
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

#[cfg(test)]
mod tests {
    use std::collections::HashSet;
    use super::*;

    #[test]
    fn test_hash_sigma() {
        let mut tested = HashSet::<[u8; 3]>::default();
        // We can guarantee at least 5232 unique colors with this method
        for i in 0..5232 as Sigma {
            let rgb = sigma_to_rgb(i);
            assert!(!tested.contains(&rgb));
            tested.insert(rgb);
        }
    }
}