use crate::boundary::{FixedBoundary, Boundary};
use crate::cell::Cell;
use crate::edge::{Edge, EdgeBook};
use crate::environment::LatticeEntity::*;
use crate::lattice::Lattice;
use crate::pos::{Pos2D, Rect};

pub struct Environment {
    pub cell_lattice: Lattice<i16, FixedBoundary<usize>>,
    cell_vec: Vec<Cell>,
    pub edge_book: EdgeBook,
    // TODO: this should be a MooreNeighbourhood field that implements Neighbourhood
    pub neigh_r: u8
}
impl Environment {
    pub fn new(width: usize, height: usize, neigh_r: u8) -> Self {
        let mut me = Self {
            cell_lattice: Lattice::new(FixedBoundary::new(Rect::new(
                (0, 0).into(),
                (width, height).into()
            ))),
            cell_vec: vec![],
            edge_book: EdgeBook::new(),
            neigh_r
        };
        me.make_border();
        me
    }

    pub fn width(&self) -> usize {
        self.cell_lattice.width()
    }

    pub fn height(&self) -> usize {
        self.cell_lattice.height()
    }

    pub fn get_entity(&self, sigma: i16) -> LatticeEntity<&Cell> {
        if sigma == Medium.discriminant() {
            return Medium;
        }
        if sigma == Solid.discriminant() {
            return Solid;
        }
        SomeCell(&self.cell_vec[sigma as usize - LatticeEntity::first_sigma() as usize])
    }

    pub fn get_entity_mut(&mut self, sigma: i16) -> LatticeEntity<&mut Cell> {
        if sigma == Medium.discriminant() {
            return Medium;
        }
        if sigma == Solid.discriminant() {
            return Solid;
        }
        SomeCell(&mut self.cell_vec[sigma as usize - LatticeEntity::first_sigma() as usize])
    }

    // TODO: ensure this makes sense for neigh_r > 1
    pub fn edge_per_pos(&self) -> f64 {
        self.neigh_r as f64
    }
    
    pub fn n_cells(&self) -> usize {
        self.cell_vec.len()
    }

    pub fn spawn_rect_cell(&mut self, rect: Rect<usize>, target_area: u32) -> Option<&Cell> {
        let mut cell_area = 0u32;
        let sigma = self.n_cells() as i16 + LatticeEntity::first_sigma();
        for p in rect.iter_positions() {
            if !self.cell_lattice.bound.inbounds(p) || self.cell_lattice[p] != Medium.discriminant() {
                continue;
            }
            self.cell_lattice[p] = sigma;
            for neigh in self.cell_lattice.bound.validate_positions(p.moore_neighs(self.neigh_r)) {
                let edge = Edge::new(p, neigh, self.neigh_r).unwrap();
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
        for neigh in self.cell_lattice.bound.validate_positions(pos.moore_neighs(self.neigh_r)) {
            let edge = Edge::new(pos, neigh, self.neigh_r).unwrap();
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
    pub fn first_sigma() -> i16 {
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
    pub fn discriminant(&self) -> i16 {
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
mod tests {
    use super::*;

    // Setup functions
    fn empty_env() -> Environment {
        Environment::new(100, 100, 1)
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
