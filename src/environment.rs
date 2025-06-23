use crate::boundary::{Boundary, UnsafePeriodicBoundary};
use crate::cell::Cell;
use crate::edge::{Edge, EdgeBook};
use crate::environment::LatticeEntity::*;
use crate::lattice::Lattice;
use crate::neighbourhood::{MooreNeighbourhood, Neighbourhood};
use crate::pos::{Pos2D, Rect};

pub type Sigma = i32;

pub struct Environment {
    pub cell_lattice: Lattice<Sigma, UnsafePeriodicBoundary<isize>>,
    cell_vec: Vec<Cell>,
    pub edge_book: EdgeBook,
    pub neighbourhood: MooreNeighbourhood
}
impl Environment {
    pub fn new(width: usize, height: usize, neigh_r: u8, enclose: bool) -> Self {
        let rect = Rect::new(
            (0, 0).into(),
            (width as isize, height as isize).into()
        );
        let mut env = Self {
            cell_lattice: Lattice::new(UnsafePeriodicBoundary::new(rect)),
            cell_vec: vec![],
            edge_book: EdgeBook::new(),
            neighbourhood: MooreNeighbourhood::new(neigh_r)
        };
        if enclose {
            env.make_border();
        }
        env
    }

    pub fn width(&self) -> usize {
        self.cell_lattice.width()
    }

    pub fn height(&self) -> usize {
        self.cell_lattice.height()
    }

    pub fn get_entity(&self, sigma: Sigma) -> LatticeEntity<&Cell> {
        if sigma == Medium.discriminant() {
            return Medium;
        }
        if sigma == Solid.discriminant() {
            return Solid;
        }
        SomeCell(&self.cell_vec[sigma as usize - LatticeEntity::first_sigma() as usize])
    }

    pub fn get_entity_mut(&mut self, sigma: Sigma) -> LatticeEntity<&mut Cell> {
        if sigma == Medium.discriminant() {
            return Medium;
        }
        if sigma == Solid.discriminant() {
            return Solid;
        }
        SomeCell(&mut self.cell_vec[sigma as usize - LatticeEntity::first_sigma() as usize])
    }
    
    pub fn n_cells(&self) -> usize {
        self.cell_vec.len()
    }

    pub fn spawn_rect_cell(&mut self, rect: Rect<usize>, target_area: u32) -> Option<&Cell> {
        let mut cell_area = 0u32;
        let sigma = self.n_cells() as Sigma + LatticeEntity::first_sigma();
        for pos in rect.iter_positions() {
            let trans_pos = self.cell_lattice.bound.valid_pos(pos.into());
            if trans_pos.is_none() {
                continue;
            }
            let valid_pos: Pos2D<usize> = trans_pos.unwrap().into();
            if self.cell_lattice[valid_pos] != Medium.discriminant() {
                continue
            }
            self.cell_lattice[valid_pos] = sigma;
            let valid_neighs = self.cell_lattice
                .bound
                .valid_positions(self.neighbourhood.neighbours(pos.into()))
                .map(|neigh| neigh.into());
            for neigh in valid_neighs {
                let edge = Edge::new(valid_pos, neigh);
                if self.cell_lattice[neigh] != sigma {
                    self.edge_book.insert(edge);
                } else { 
                    self.edge_book.remove(&edge);
                }
            }
            cell_area += 1;
        }
        if cell_area == 0 { 
            return None;
        }
        self.cell_vec.push(Cell::new(cell_area, target_area));
        Some(self.get_entity(sigma).unwrap_cell())
    }
    
    pub fn spawn_solid(&mut self, positions: impl Iterator<Item = Pos2D<usize>>) -> usize {
        let mut area = 0;
        for pos in positions {
            if self.cell_lattice[pos] != Medium.discriminant() {
                continue
            }
            self.cell_lattice[pos] = Solid.discriminant();
            area += 1;
        }
        area
    }
    
    pub fn make_border(&mut self) {
        let mut border_positions = Vec::<Pos2D<usize>>::new();
        for x in 0..self.width() {
            border_positions.push((x, 0).into());
        }
        for y in 1..self.height() {
            border_positions.push((self.width() - 1, y).into());
        }
        if self.width() > 1 {
            for y in (1..self.height() - 1).rev() {
                border_positions.push((0, y).into());
            }
        }
        if self.height() > 1 {
            for x in (0..self.width() - 1).rev() {
                border_positions.push((x, self.height() - 1).into());
            }
        }
        
        self.spawn_solid(border_positions.into_iter());
    }
    
    pub fn update_edges(&mut self, pos: Pos2D<usize>) -> (u16, u16) {
        let mut removed = 0;
        let mut added = 0;
        let sigma = self.cell_lattice[pos];
        let valid_neighs = self.cell_lattice
            .bound
            .valid_positions(self.neighbourhood.neighbours(pos.into()))
            .map(|neigh| neigh.into());
        for neigh in valid_neighs {
            let edge = Edge::new(pos, neigh);
            let sigma_neigh = self.cell_lattice[neigh];
            if sigma == sigma_neigh {
                self.edge_book.remove(&edge);
                // Also representing the inverse edge
                removed += 2;
            // Since we filtered Medium, Medium before, this should only be 0 when one sigma is 1 and the other -1
            // Ideally we should test for the cases more explicitly, but I couldnt figure out an easy way to do that
            } else if sigma + sigma_neigh >= 0 && self.edge_book.insert(edge) {
                // Also representing the inverse edge
                added += 2;
            }
        }
        (removed, added)
    }
}

/// This enum represents anything that can be on the cell lattice.
#[derive(Debug, Copy, Clone)]
pub enum LatticeEntity<C> {
    Solid,
    Medium,
    SomeCell(C),
}
impl<C> LatticeEntity<C> {
    pub fn map<D, F: FnOnce(C) -> D>(self, f: F) -> LatticeEntity<D> {
        match self {
            SomeCell(c) => SomeCell(f(c)),
            Medium => Medium,
            Solid => Solid,
        }
    }
}

impl LatticeEntity<()> {
    pub fn first_sigma() -> Sigma {
        SomeCell(()).discriminant()
    }

    // There is another way to obtain these according to the docs:
    // https://doc.rust-lang.org/core/mem/fn.discriminant.html
    // I've benchmarked and it doesnt make a difference
    // If in the future sigma becomes a cell property, we can implement `LatticeEntity<&Cell>::as_sigma()` and replace
    // most references to this function with that.
    /// This returns a unique `i16` discriminant for each possible type of `LatticeEntity`.
    ///
    /// These values are used as sigmas in the cell lattice, except for the discriminant for `SomeCell`.
    pub fn discriminant(&self) -> Sigma {
        match self {
            SomeCell(_) => 1,
            Medium => 0,
            Solid => -1
        }
    }
}

impl<C: std::fmt::Debug> LatticeEntity<C> {
    pub fn unwrap_cell(self) -> C {
        match self {
            SomeCell(cell) => cell,
            _ => panic!("called `LatticeEntity::unwrap_cell()` on a `{:?}` value", self)
        }
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;

    // Setup functions
    pub fn empty_env() -> Environment {
        Environment::new(100, 100, 1, false)
    }
    
    fn env_with_cell() -> Environment {
        let mut env = empty_env();
        env.spawn_rect_cell(
            Rect::new(
                Pos2D::new(10, 10),
                Pos2D::new(20, 20)
            ),
            100
        );
        env
    }
    
    #[test]
    fn test_spawn_rect_cell() {
        let mut env = env_with_cell();
        assert_eq!(env.edge_book.len(), 8 * 4 * 3 + 4 * 5);
        env.spawn_rect_cell(
            Rect::new(
                Pos2D::new(15, 15),
                Pos2D::new(25, 25)
            ),
            10
        );
        assert_eq!(env.get_entity(1).unwrap_cell().area, 100);
        assert_eq!(env.get_entity(2).unwrap_cell().area, 75);
    }

    #[test]
    fn test_lattice_entity_discriminant() {
        assert_eq!(1, SomeCell(()).discriminant());
        assert_eq!(0, Medium.discriminant());
        assert_eq!(-1, Solid.discriminant());
    }
}
