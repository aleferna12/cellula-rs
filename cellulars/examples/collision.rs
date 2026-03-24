use cellulars::io::write::image::movie_window::MovieWindow;
use cellulars::io::write::image::plot::{Plot, SpinPlot};
use cellulars::prelude::*;
use image::RgbaImage;
use palette::Srgba;
use rand::rngs::ThreadRng;
use cellulars::copy_bias::DirectionBias;

const W: usize = 400;
const H: usize = 400;
const SIDE: usize = 50;

fn main() -> Result<(), minifb::Error> {
    // Initialize fixed boundary conditions so that the cell can collide with the wall
    let boundary = FastPeriodicBoundary::new(Rect::new(
        Pos::new(0., 0.),
        Pos::new(W as f64, H as f64)
    ));
    // Initialize an empty environment
    let mut env = Environment::new(
        cell_container![],
        Lattice::new(W, H),
        MooreNeighborhood::new(1),
        Boundaries::new(boundary.clone())
    );
    // Spawn a cell in a rectangular region of the environment
    let cell_rect = Rect::new(
        Pos::new(W / 8 - SIDE / 2, H / 2 - SIDE / 2),
        Pos::new(W / 8 + SIDE / 2, H / 2 + SIDE / 2)
    );
    env.spawn_cell(Cell::new_empty(cell_rect.area() as u32), cell_rect.iter_positions());

    // Spawn the object against which the cell will collide
    let solid_rect = Rect::new(
        Pos::new(W - 10, 0),
        Pos::new(W, H)
    );
    env.spawn_solid(solid_rect.iter_positions());

    // Initialize the Potts algorithm used to update the environment
    let mut potts = EdgePotts {
        boltz_t: 10.,
        size_lambda: 10.,
        adhesion: StaticAdhesion {
            cell_energy: 10.,
            medium_energy: 10.,
            solid_energy: 10.
        },
        bias: Biases {
            dir: DirectionBias {
                // Crank this value up to see some wacky physics
                lambda: 10.,
                dir_params: DirectionalOptions {
                    protrusions: true,
                    retractions: true,
                    contact_inhibition: false
                },
                boundary,
            }
        }
    };
    let mut rng = ThreadRng::default();
    let spin_plot = SpinPlot {
        solid_color: Srgba::new(1., 1., 1., 1.),
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

struct Biases {
    dir: DirectionBias<FastPeriodicBoundary<FloatType>>
}

impl CopyBias<Environment<Cell>> for Biases {
    fn bias(&self, pos_source: Pos<usize>, pos_target: Pos<usize>, context: &Environment<Cell>) -> FloatType {
        // The cell is commited to going to the right
        self.dir.bias(
            pos_source,
            pos_target,
            &DirectionContext{
                cell_lattice: &context.cell_lattice,
                angle: 0.0f64.to_radians()
            })
    }
}