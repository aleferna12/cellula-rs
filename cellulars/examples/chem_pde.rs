use cellulars::io::write::image::lerper::Lerper;
use cellulars::io::write::image::movie_window::MovieWindow;
use cellulars::io::write::image::plot::{srgba_to_rgba, Plot, SpinPlot};
use cellulars::prelude::*;
use image::RgbaImage;
use palette::{IntoColor, Oklab, Srgba};
use rand::RngExt;
use std::mem;

const W: usize = 300;
const H: usize = 300;
const SIDE: u32 = 20;

fn main() -> Result<(), minifb::Error> {
    let mut pond = Pond {
        env: ChemEnvironment {
            env: Environment::new(
                cell_container![],
                Lattice::new(W, H),
                MooreNeighborhood::new(1),
                Boundaries::new(UnsafePeriodicBoundary::new(
                    Rect::new(
                        Pos::new(0., 0.),
                        Pos::new(W as f64, H as f64),
                    )
                ))
            ),
            current_chem: Lattice::new(W, H),
            new_chem: Lattice::new(W, H),
            // Max is 1 / n_neighbors = 0.125
            diffuse_rate: 0.125,
            secrete_rate: 0.025,
            decay_rate: 0.0002,
        },
        potts: EdgePotts {
            boltz_t: 5.,
            size_lambda: 10.,
            adhesion: StaticAdhesion {
                cell_energy: 10.,
                medium_energy: 10.,
                solid_energy: 10.,
            },
            bias: Biases {
                chem_lambda: 50.,
            }
        },
        rng: Default::default(),
    };

    for _ in 0..25 {
        let x = pond.rng.random_range(0..W);
        let y = pond.rng.random_range(0..H);
        let cell_rect = Rect::new(
            Pos::new(x, y),
            Pos::new(x + SIDE as usize, y + SIDE as usize),
        );
        let valid: Vec<_> = cell_rect.iter_positions().filter_map(|pos| {
            let pos_isize = pos.cast_as();
            pond.env.env.bounds.lattice_boundary.valid_pos(pos_isize).map(|valid_pos| valid_pos.cast_as())
        }).collect();
        pond.env.env.spawn_cell(Cell::new_empty(SIDE * SIDE), valid);
    }

    let mut window = MovieWindow::new(W, H)?;
    let mut image = RgbaImage::new(W as u32, H as u32);
    let spin_plot = SpinPlot {
        solid_color: Default::default(),
        medium_color: None
    };
    let chem_plot = ChemPlot {
        lerper: Lerper {
            min_color: Default::default(),
            max_color: Srgba::new(1., 1., 1., 1.).into_color(),
        }
    };

    for step in 0..1_000_000u32 {
        if step.is_multiple_of(10) {
            chem_plot.plot(&mut pond.env.current_chem, &mut image);
            spin_plot.plot(&mut pond.env.env, &mut image);
            window.update(&mut image)?;
            image.fill(0);
        }
        pond.step();
    }
    Ok(())
}

struct ChemEnvironment {
    env: Environment<Cell>,
    current_chem: Lattice<f64>,
    new_chem: Lattice<f64>,
    diffuse_rate: f64,
    secrete_rate: f64,
    decay_rate: f64
}

impl ChemEnvironment {
    // Updates the chemical signal on the environment with simple PDE system
    // In a real model you would parallelize this ofc
    fn update_chem(&mut self) {
        for pos in self.current_chem.iter_positions() {
            let current_chem = self.current_chem[pos];
            let mut neigh_chem = 0.;
            for neigh in self.env.valid_neighbors(pos) {
                neigh_chem += self.current_chem[neigh];
            }
            let diffusion_diff = neigh_chem - self.env.neighborhood.n_neighs() as f64 * current_chem;
            let mut new_chem = (current_chem + self.diffuse_rate * diffusion_diff) * (1. - self.decay_rate);
            if let Spin::Some(_) = self.env.cell_lattice[pos] {
                new_chem += self.secrete_rate;
            }
            self.current_chem[pos] = new_chem;
        }
        mem::swap(&mut self.current_chem, &mut self.new_chem);
    }
}

impl AsEnv for ChemEnvironment {
    type Cell = Cell;

    fn env(&self) -> &Environment<Self::Cell, impl Neighborhood, impl ToLatticeBoundary> {
        &self.env
    }

    fn env_mut(&mut self) -> &mut Environment<Self::Cell, impl Neighborhood, impl ToLatticeBoundary> {
        &mut self.env
    }
}

impl TransferPosition for ChemEnvironment {
    fn transfer_position(&mut self, pos: Pos<usize>, to: Spin) -> EdgesUpdate {
        self.env.transfer_position(pos, to)
    }
}

pub struct Biases {
    chem_lambda: f64
}

impl CopyBias<ChemEnvironment> for Biases {
    fn bias(&self, pos_source: Pos<usize>, pos_target: Pos<usize>, env: &ChemEnvironment) -> f64 {
        -ChemotaxisBias { lambda: self.chem_lambda }.bias(pos_source, pos_target, &env.current_chem)
    }
}

struct Pond {
    env: ChemEnvironment,
    potts: EdgePotts<StaticAdhesion, Biases>,
    rng: rand::rngs::ThreadRng
}

impl Step for Pond {
    fn step(&mut self) {
        self.potts.step(&mut self.env, &mut self.rng);
        self.env.update_chem();
    }
}

struct ChemPlot {
    lerper: Lerper<Oklab<f64>>
}

impl Plot<Lattice<f64>> for ChemPlot {
    fn plot(&self, plottable: &Lattice<f64>, image: &mut RgbaImage) {
        for pos in plottable.iter_positions() {
            let chem = plottable[pos];
            let Ok(color) = self.lerper.lerp((chem / 30.).min(1.)) else {
                continue;
            };
            image.put_pixel(pos.x as u32, pos.y as u32, srgba_to_rgba(color.into_color()))
        }
    }
}
