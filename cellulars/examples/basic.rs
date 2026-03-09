use cellulars::io::write::image::movie_window::MovieWindow;
use cellulars::io::write::image::plot::{Plot, SpinPlot};
use cellulars::prelude::*;
use image::RgbaImage;
use rand::rngs::ThreadRng;

const W: usize = 300;
const H: usize = 300;
const SIDE: usize = 100;

fn main() -> Result<(), minifb::Error> {
    // Initialize periodic boundary conditions
    let boundary = UnsafePeriodicBoundary::new(Rect::new(
        Pos::new(0., 0.),
        Pos::new(W as f64, H as f64)
    ));
    // Initialize an empty environment
    let mut env = Environment::new(
        cell_container![],
        Lattice::new(W, H),
        MooreNeighborhood::new(1),
        Boundaries::new(boundary)
    );
    // Spawn a cell in a rectangular region of the environment
    let cell_rect = Rect::new(
        Pos::new(W / 2 - SIDE / 2, H / 2 - SIDE / 2),
        Pos::new(W / 2 + SIDE / 2, H / 2 + SIDE / 2)
    );
    env.spawn_cell(Cell::new_empty(cell_rect.area() as u32), cell_rect.iter_positions());

    // Initialize the Potts algorithm used to update the environment
    let mut potts = EdgePotts {
        boltz_t: 10.,
        size_lambda: 10.,
        adhesion: StaticAdhesion {
            cell_energy: 10.,
            medium_energy: 10.,
            solid_energy: 10.
        },
        bias: NoBias
    };
    let mut rng = ThreadRng::default();
    let spin_plot = SpinPlot {
        solid_color: Default::default(),
        medium_color: None,
    };
    let mut image = RgbaImage::new(W as u32, H as u32);
    let mut window = MovieWindow::new(W, H)?;

    // Run 1m steps of the simulation
    for step in 0..1_000_000u32 {
        if step.is_multiple_of(10) {
            spin_plot.plot(&env, &mut image);
            window.update(&image)?;
            image.fill(0);
        }
        potts.step(&mut env, &mut rng);
    }
    Ok(())
}