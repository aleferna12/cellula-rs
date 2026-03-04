use cellulars::cell_container;
use cellulars::io::write::image::movie_window::MovieWindow;
use cellulars::io::write::image::plot::{Plot, SpinPlot};
use cellulars::positional::boundaries::Boundaries;
use cellulars::prelude::*;
use image::RgbaImage;
use palette::Srgba;
use rand::RngExt;

const W: usize = 600;
const H: usize = 600;
const SIDE: u32 = 10;

fn main() -> Result<(), minifb::Error> {
    let mut pond = Pond {
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
        potts: Potts {
            boltz_t: 5.,
            size_lambda: 20.,
            adhesion: StaticAdhesion {
                cell_energy: 5.,
                medium_energy: 10.,
                solid_energy: 10.,
            }
        },
        rng: Default::default(),
        step: 0,
    };

    let cell_rect = Rect::new(
        Pos::new(W / 2, H / 2),
        Pos::new(W / 2 + SIDE as usize, H / 2 + SIDE as usize),
    );
    pond.env.spawn_cell(DividingCell::new_empty(SIDE * SIDE), cell_rect.iter_positions());

    let mut window = MovieWindow::new(W, H)?;
    let mut image = RgbaImage::new(W as u32, H as u32);
    let plot = SpinPlot {
        solid_color: Srgba::new(1., 1., 1., 1.),
        medium_color: None
    };

    for step in 0..1_000_000u32 {
        if step.is_multiple_of(10) {
            plot.plot(&mut pond.env, &mut image);
            window.update(&mut image)?;
            image.fill(0);
        }
        pond.step();
    }
    Ok(())
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

struct Potts {
    boltz_t: f64,
    size_lambda: f64,
    adhesion: StaticAdhesion,
}

impl PottsAlgorithm for Potts {
    type Environment = Environment<DividingCell>;

    fn boltz_t(&self) -> f64 {
        self.boltz_t
    }

    fn size_lambda(&self) -> f64 {
        self.size_lambda
    }

    fn delta_hamiltonian_adhesion(
        &self,
        spin_source: Spin,
        spin_target: Spin,
        neigh_spin: impl IntoIterator<Item=Spin>,
        _env: &Self::Environment
    ) -> f64 {
        let mut energy = 0.0;
        for neigh in neigh_spin {
            energy -= self.adhesion.adhesion_energy(neigh, spin_target, &());
            energy += self.adhesion.adhesion_energy(neigh, spin_source, &());
        }
        energy
    }
}

struct Pond {
    env: Environment<DividingCell>,
    potts: Potts,
    rng: rand::rngs::ThreadRng,
    step: u32
}

impl Pond {
    fn divide_cell(&mut self, cell_index: u32) {
        let rel_cell = &self.env.cells[cell_index];
        let angle = self.rng.random::<f64>() * std::f64::consts::PI;
        let center = rel_cell.cell.center();
        let positions: Vec<_> = self
            .env
            .search_cell_box(rel_cell, 1.5)
            .iter()
            .copied()
            .filter(|pos| {
                angle.sin() * (pos.y as f64 - center.y) - angle.cos() * (pos.x as f64 - center.x) < 0.
            })
            .collect();

        let new = rel_cell.cell.birth();
        let new_index = self.env.cells.add(new).index;
        for pos in positions {
            self.env.grant_position(pos, Spin::Some(new_index));
        }
        self.env.cells[cell_index].cell.reset_target_area();
    }
}

impl Step for Pond {
    fn step(&mut self) {
        self.potts.step(&mut self.env, &mut self.rng);
        if self.step.is_multiple_of(32) {
            let mut to_divide = vec![];
            for rel_cell in self.env.cells.iter_non_empty_mut() {
                rel_cell.cell.cell.target_area += 1;
                if rel_cell.cell.cell.area() > SIDE.pow(2) * 2 {
                    to_divide.push(rel_cell.index);
                }
            }
            for index in to_divide {
                self.divide_cell(index);
            }
        }
        self.step += 1;
    }
}
