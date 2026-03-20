use cellulars::io::write::image::plot::{BorderPlot, Plot, SpinPlot};
use cellulars::prelude::*;
use image::RgbaImage;
use rand::{RngExt, SeedableRng};
use rand_xoshiro::Xoshiro256StarStar;

const W: usize = 1000;
const H: usize = 1000;
const SIDE: u32 = 10;
const GROWTH_PERIOD: u32 = 32;
const SAVE_PERIOD: u32 = 1_000;
const MAX_POP: u32 = 2000;

pub fn main() {
    let mut pond = Pond {
        env: Environment::new(
            cell_container![],
            Lattice::new(W, H),
            MooreNeighborhood::new(1),
            Boundaries::new(FastPeriodicBoundary::new(
                Rect::new(
                    Pos::new(0., 0.),
                    Pos::new(W as f64, H as f64),
                )
            ))
        ),
        potts: EdgePotts {
            boltz_t: 10.,
            size_lambda: 4.,
            adhesion: StaticAdhesion {
                // Morpheus does this differently to us (cell adh is not calculated per cell)
                cell_energy: 5.,
                medium_energy: 10.,
                solid_energy: 10.,
            },
            bias: NoBias,
        },
        rng: Xoshiro256StarStar::seed_from_u64(rand::random()),
        step: 0,
    };

    let cell_rect = Rect::new(
        Pos::new(W / 2, H / 2),
        Pos::new(W / 2 + SIDE as usize, H / 2 + SIDE as usize),
    );
    pond.env.spawn_cell(DividingCell::new_empty(SIDE * SIDE), cell_rect.iter_positions());

    let mut image = RgbaImage::new(W as u32, H as u32);
    let spin_plot = SpinPlot {
        solid_color: Default::default(),
        medium_color: None
    };
    let border_plot = BorderPlot {
        color: Default::default(),
    };

    for step in 0..1_000_000u32 {
        if pond.env.cells.n_cells() >= MAX_POP {
            println!("Population size reached {} on time step {step}", pond.env.cells.n_cells());
            break;
        }
        if step.is_multiple_of(SAVE_PERIOD) {
            println!("{step}: {} cells", pond.env.cells.n_cells());
            spin_plot.plot(&mut pond.env, &mut image);
            border_plot.plot(&mut pond.env, &mut image);
            let res = image.save(format!("benchmark/out/cellulars/{step}.png"));
            if let Err(e) = res {
                println!("Failed to save image with error: {e}");
            };
            image.fill(0);
        }
        pond.step();
    }
    println!("Reached end of simulation");
}

struct DividingCell {
    newborn_target_area: u32,
    cell: Cell
}

impl DividingCell {
    fn new_empty(target_area: u32) -> EmptyCell<Self> {
        EmptyCell::new_unchecked(DividingCell {
            newborn_target_area: target_area,
            cell: Cell::new_empty(target_area).into_cell()
        })
    }

    fn reset_target_area(&mut self) {
        self.cell.target_area = self.newborn_target_area;
    }
}

impl Cellular for DividingCell {
    fn target_area(&self) -> u32 {
        self.cell.target_area()
    }

    fn area(&self) -> u32 {
        self.cell.area()
    }

    fn shift_position(
        &mut self,
        pos: Pos<usize>,
        add: bool,
        bound: &impl Boundary<Coord = f64>
    ) -> Result<(), ShiftError> {
        self.cell.shift_position(pos, add, bound)
    }
}

impl Alive for DividingCell {
    fn is_alive(&self) -> bool {
        self.cell.is_alive()
    }

    fn apoptosis(&mut self) {
        self.cell.apoptosis()
    }

    fn birth(&self) -> EmptyCell<Self> {
        let mut cell = self.cell.birth().into_cell();
        cell.target_area = self.newborn_target_area;
        EmptyCell::new_unchecked(Self { newborn_target_area: self.newborn_target_area, cell })
    }
}

impl Empty for DividingCell {
    fn empty_default() -> EmptyCell<Self> {
        EmptyCell::new_unchecked(Self { newborn_target_area: 0, cell: Cell::empty_default().into_cell() })
    }

    fn is_empty(&self) -> bool {
        self.cell.is_empty()
    }
}

impl HasCenter for DividingCell {
    fn center(&self) -> Pos<f64> {
        self.cell.center()
    }
}

struct Pond {
    env: Environment<DividingCell>,
    potts: EdgePotts<StaticAdhesion, NoBias>,
    rng: Xoshiro256StarStar,
    step: u32
}

impl Pond {
    // Divides the by tracing a line through its center in a random angle
    fn divide_cell(&mut self, cell_index: u32) {
        let rel_cell = &self.env.cells[cell_index];
        let angle = self.rng.random::<f64>() * std::f64::consts::PI;
        let (sin, cos) = angle.sin_cos();
        let center = rel_cell.cell.center();
        let positions: Vec<_> = self
            .env
            .search_cell_box(rel_cell, 1.5)
            .iter()
            .copied()
            .filter(|pos| {
                sin * (pos.y as f64 - center.y) - cos * (pos.x as f64 - center.x) < 0.
            })
            .collect();

        let new = rel_cell.cell.birth();
        let new_index = self.env.cells.add(new).index;
        for pos in positions {
            self.env.transfer_position(pos, Spin::Some(new_index));
        }
        self.env.cells[cell_index].cell.reset_target_area();
    }
}

impl Step for Pond {
    fn step(&mut self) {
        self.potts.step(&mut self.env, &mut self.rng);
        let mut to_divide = vec![];
        for rel_cell in self.env.cells.iter_non_empty() {
            // If they have grown enough, it's time to divide!
            if rel_cell.cell.cell.area() > SIDE * SIDE * 2 {
                to_divide.push(rel_cell.index);
            }
        }
        for index in to_divide {
            self.divide_cell(index);
        }
        if self.step.is_multiple_of(GROWTH_PERIOD) {
            for rel_cell in self.env.cells.iter_non_empty_mut() {
                rel_cell.cell.cell.target_area += 1;
            }
        }
        self.step += 1;
    }
}

