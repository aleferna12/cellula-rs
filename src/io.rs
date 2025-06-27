// TODO!:
//  Its time that this module gets a set of traits
//  The outdir parameter of parameters need to be stored somewhere as an absolute path
//  I also want to be able to both save images or display them in real time

use std::hash::{DefaultHasher, Hash, Hasher};
use std::path::Path;
use image::RgbImage;
use std::io;
use crate::environment::{Environment, Sigma};

pub(crate) static IMAGES_PATH: &str = "images";
pub(crate) static CONFIG_COPY_PATH: &str = "config.toml";

pub fn create_directories(outpath: impl AsRef<Path>, replace_outdir: bool) -> io::Result<()> {
    let outpath = outpath.as_ref();
    let outdir_exists = outpath.try_exists()?;
    if outdir_exists {
        if replace_outdir {
            log::info!("Cleaning contents of '{}'", outpath.display());
            std::fs::remove_dir_all(outpath)?;
        } else {
            return Err(io::Error::new(
                io::ErrorKind::AlreadyExists,
                "`outdir` already exists and `replace_outdir` is `false`"
            ));
        }
    }
    std::fs::create_dir_all(outpath)?;
    std::fs::create_dir(outpath.join(IMAGES_PATH))
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
    fn test_sigma_to_rgb() {
        let mut tested = HashSet::<[u8; 3]>::default();
        // We can guarantee at least 5232 unique colors with this method
        for i in 0..5232 as Sigma {
            let rgb = sigma_to_rgb(i);
            assert!(!tested.contains(&rgb));
            tested.insert(rgb);
        }
    }
}
