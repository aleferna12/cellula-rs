use cellulars::io::write::image::movie_window::MovieWindow;
use cellulars::io::write::image::plot::{Plot, SpinPlot};
use cellulars::prelude::*;
use image::{GrayImage, ImageReader, Rgba, RgbaImage};
use rand::rngs::ThreadRng;

const W: usize = 500;
const H: usize = 500;
const SIDE: usize = 215;

fn main() -> Result<(), minifb::Error> {
    // Initialize periodic boundary conditions
    let boundary = FastPeriodicBoundary::new(Rect::new(
        Pos::new(0., 0.),
        Pos::new(W as f64, H as f64)
    ));
    let img = ImageReader::open("cellulars/examples/data/butterfly.png")
        .expect("failed to open img file")
        .decode()
        .expect("failed to decode img")
        .into_luma8();

    // Initialize an empty environment
    let mut env = ShapeEnvironment::new(
        Environment::new(
            cell_container![],
            Lattice::new(W, H),
            MooreNeighborhood::new(1),
            Boundaries::new(boundary),
        ),
        &img
    );
    // Spawn a cell in a rectangular region of the environment
    let cell_rect = Rect::new(
        Pos::new(W - SIDE, H / 2 - SIDE / 2),
        Pos::new(W, H / 2 + SIDE / 2)
    );
    env.env.spawn_cell(Cell::new_empty(cell_rect.area() as u32), cell_rect.iter_positions());

    // Initialize the Potts algorithm used to update the environment
    let mut potts = EdgePotts {
        boltz_t: 10.,
        size_lambda: 10.,
        adhesion: StaticAdhesion {
            cell_energy: 10.,
            medium_energy: 10.,
            solid_energy: 10.
        },
        bias: Biases { shape_lambda: 10. }
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
            ShapePlot {}.plot(&env.shape_lat, &mut image);
            spin_plot.plot(&env.env, &mut image);
            window.update(&image)?;
            image.fill(0);
        }
        potts.step(&mut env, &mut rng);
    }
    Ok(())
}

struct ShapeEnvironment {
    env: Environment<Cell>,
    shape_lat: Lattice<f64>,
    fill_score: f64
}

impl ShapeEnvironment {
    fn new(env: Environment<Cell>, shape_img: &GrayImage) -> Self {
        let mut lat = Lattice::new(shape_img.width() as usize, shape_img.height() as usize);
        for x in 0..lat.width() {
            for y in 0..lat.height() {
                lat[(x, y).into()] = shape_img[(x as u32, y as u32)].0[0] as f64 / 255.;
            }
        }
        ShapeEnvironment { env, shape_lat: lat, fill_score: 0.0 }
    }
}

impl TransferPosition for ShapeEnvironment {
    fn transfer_position(&mut self, pos: Pos<usize>, to: Spin) -> EdgesUpdate {
        match to {
            Spin::Some(_) => self.fill_score += self.shape_lat[pos],
            Spin::Medium => self.fill_score -= self.shape_lat[pos],
            Spin::Solid => unreachable!()
        }
        self.env.transfer_position(pos, to)
    }
}

impl AsEnv for ShapeEnvironment {
    type Cell = Cell;

    fn env(&self) -> &Environment<Self::Cell, impl Neighborhood, impl ToLatticeBoundary> {
        &self.env
    }

    fn env_mut(&mut self) -> &mut Environment<Self::Cell, impl Neighborhood, impl ToLatticeBoundary> {
        &mut self.env
    }
}

struct Biases {
    shape_lambda: f64
}

impl CopyBias<ShapeEnvironment> for Biases {
    fn bias(&self, pos_source: Pos<usize>, pos_target: Pos<usize>, context: &ShapeEnvironment) -> f64 {
        let spin_source = context.env.cell_lattice[pos_source];
        let value_target = context.shape_lat[pos_target];
        // We only check source because assuming one cell
        self.shape_lambda * match spin_source {
            Spin::Some(_) => value_target - context.fill_score / context.env.cells[0].cell.area() as f64,
            Spin::Medium => 0.,
            Spin::Solid => 0.
        }
    }
}

struct ShapePlot;

impl Plot<Lattice<f64>> for ShapePlot {
    fn plot(&self, plottable: &Lattice<f64>, image: &mut RgbaImage) {
        for pos in plottable.iter_positions() {
            let val = (255. * plottable[pos]).round() as u8;
            image.put_pixel(pos.x as u32, pos.y as u32, Rgba::from([val, val, val, 255]));
        }
    }
}