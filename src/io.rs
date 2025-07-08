use std::collections::HashSet;
use std::hash::{DefaultHasher, Hash, Hasher};
use std::path::Path;
use image::{EncodableLayout, RgbImage};
use std::io;
use imageproc::drawing::draw_line_segment_mut;
use crate::boundary::Boundary;
use crate::cell::RelCell;
use crate::constants::Spin;
use crate::environment::{Environment, LatticeEntity};
use crate::parameters::{IoParameters, MovieParameters};
use crate::pos::Pos2D;

pub(crate) static IMAGES_PATH: &str = "images";
pub(crate) static CONFIG_COPY_PATH: &str = "config.toml";

pub struct IoManager {
    pub outdir: String,
    pub replace_outdir: bool,
    pub show_attached_cells: bool,
    pub show_cell_centers: bool,
    pub image_period: u32,
    pub image_format: String,
    pub movie_maker: Option<MovieMaker>
}

impl IoManager {
    pub fn image_io(
        &mut self,
        time_step: u32,
        env: &Environment,
        clone_pairs: &HashSet<(Spin, Spin)>
    ) {
        let mut frame = None;
        let movie_update = if let Some(mm) = &self.movie_maker {
            time_step % mm.frame_period == 0 && mm.window_works()
        } else {
            false
        };
        if movie_update {
            frame = Some(self.simulation_image(env, clone_pairs));
            let mm = self.movie_maker.as_mut().unwrap();
            let resized = image::imageops::resize(
                frame.as_ref().unwrap(),
                mm.width,
                mm.height,
                image::imageops::Nearest,
            );
            if let Err(e) = mm.update(&resized) {
                log::warn!(
                "Failed to display simulation frame at time step {} with error `{}`",
                time_step,
                e
            );
            }
        }

        if time_step % self.image_period == 0 {
            if frame.is_none() {
                frame = Some(self.simulation_image(env, clone_pairs));
            }
            let saved = frame.unwrap().save(
                Path::new(&self.outdir)
                    .join(IMAGES_PATH)
                    .join(format!("{time_step}.{}", self.image_format.to_lowercase()))
            );
            if let Err(e) = saved {
                log::warn!("Failed to save simulation frame at time step {} with error `{}`", time_step, e);
            }
        }
    }
    
    pub fn create_directories(&self) -> io::Result<()> {
        let outpath = Path::new(&self.outdir);
        let outdir_exists = outpath.try_exists()?;
        if outdir_exists {
            if self.replace_outdir {
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

    pub fn simulation_image(&self, env: &Environment, clone_pairs: &HashSet<(Spin, Spin)>) -> RgbImage {
        let spins: Vec<_> = env
            .cell_lattice
            .iter_values()
            .flat_map(Self::spin_to_rgb)
            .collect();

        let mut image = RgbImage::from_vec(
            env.width() as u32,
            env.height() as u32,
            spins
        ).unwrap();

        if self.show_cell_centers {
            for cell in &env.cells {
                let center = env.cell_lattice.bound.valid_pos(Pos2D::new(
                    cell.center.pos().x as isize,
                    cell.center.pos().y as isize,
                ));
                if let Some(pos) = center {
                    image.put_pixel(pos.x as u32, pos.y as u32, [0, 255, 0].into());
                }
            }
        }
        
        if self.show_attached_cells {
            for (spin1, spin2) in clone_pairs.iter().copied() {
                let message = "non-cell stored as clone";
                let center1 = env.cells.get_entity(spin1).expect_cell(message).center.pos;
                let center2 = env.cells.get_entity(spin2).expect_cell(message).center.pos;
                draw_line_segment_mut(
                    &mut image,
                    (center1.x, center1.y),
                    (center2.x, center2.y),
                    [255, 0, 0].into()
                )
            }
        }
        image
    }

    /// Converts a spin into a unique color.
    ///
    /// This method guarantees 5232 unique colors, starting from this spin the colors will repeat.
    fn spin_to_rgb(spin: Spin) -> [u8; 3] {
        if spin == LatticeEntity::Medium::<&RelCell>.spin() {
            return [255, 255, 255];
        } else if spin == LatticeEntity::Solid::<&RelCell>.spin() {
            return [0, 0, 0]
        }

        let mut hasher = DefaultHasher::new();
        spin.hash(&mut hasher);
        let hashed = hasher.finish();
        [
            (hashed & 0xFF).try_into().unwrap(),
            (hashed >> 8 & 0xFF) as u8,
            (hashed >> 16 & 0xFF) as u8,
        ]
    }
}

impl TryFrom<IoParameters> for IoManager {
    type Error = <MovieMaker as TryFrom<MovieParameters>>::Error;
    fn try_from(params: IoParameters) -> Result<Self, Self::Error> {
        Ok(Self {
            outdir: params.outdir,
            replace_outdir: params.replace_outdir,
            show_attached_cells: params.show_attached_cells,
            show_cell_centers: params.show_cell_centers,
            image_period: params.image_period,
            image_format: params.image_format,
            movie_maker: if params.movie.show { 
                match MovieMaker::try_from(params.movie.clone()) { 
                    Ok(mm) => {
                        log::info!("Creating window for real-time movie display");
                        Some(mm)
                    },
                    Err(e) => {
                        log::warn!("Failed to initialise movie maker with error `{}`", e);
                        None
                    }
                }
            } else {
                None 
            },
        })
    }
}

pub struct MovieMaker {
    pub width: u32,
    pub height: u32,
    pub frame_period: u32,
    pub window: minifb::Window,
}

impl MovieMaker {
    pub fn window_works(&self) -> bool {
        self.window.is_open() && !self.window.is_key_down(minifb::Key::Escape)
    }

    pub fn update(&mut self, image: &RgbImage) -> minifb::Result<()> {
        let buffer: Vec<_> = image
            .as_bytes()
            .chunks_exact(3)
            .map(|rgb| {
                u32::from_le_bytes([rgb[2], rgb[1], rgb[0], 255])
            })
            .collect();
        self.window.update_with_buffer(&buffer, self.width as usize, self.height as usize)
    }
}

impl TryFrom<MovieParameters> for MovieMaker {
    type Error = minifb::Error;

    fn try_from(params: MovieParameters) ->  Result<Self, Self::Error> {
        let window = minifb::Window::new(
            "Evo-CPM",
            params.width as usize,
            params.height as usize,
            minifb::WindowOptions::default()
        )?;
        Ok(Self {
            width: params.width,
            height: params.height,
            frame_period: params.frame_period,
            window
        })
    }
}

#[cfg(test)]
mod tests {
    use std::collections::HashSet;
    use super::*;

    #[test]
    fn test_spin_to_rgb() {
        let mut tested = HashSet::<[u8; 3]>::default();
        // We can guarantee 5232 unique colors with this method, after that colors repeat
        for i in 0..5232 as Spin {
            let rgb = IoManager::spin_to_rgb(i);
            assert!(!tested.contains(&rgb));
            tested.insert(rgb);
        }
    }
}
