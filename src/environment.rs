use crate::cell::Cell;
use crate::edge::{Edge, EdgeBookkeeper};
use crate::lattice::{Lattice, LatticeEntity};
use crate::pos::{Pos2D, Rect};

pub struct Environment {
    pub cell_lattice: Lattice<usize>,
    cell_vec: Vec<Cell>,
    pub edge_bk: EdgeBookkeeper,
    pub neigh_r: u8
}
impl Environment {
    pub fn new(width: usize, height: usize, neigh_r: u8) -> Self {
        Self {
            cell_lattice: Lattice::new(width, height),
            cell_vec: vec![],
            edge_bk: EdgeBookkeeper::new(),
            neigh_r
        }
    }

    pub fn width(&self) -> usize {
        self.cell_lattice.width()
    }

    pub fn height(&self) -> usize {
        self.cell_lattice.height()
    }

    // Cell population functions
    pub fn get_cell(&self, sigma: usize) -> LatticeEntity<&Cell> {
        match sigma { 
            0 => LatticeEntity::Medium,
            _ => LatticeEntity::SomeCell(&self.cell_vec[sigma - 1])
        }
    }
    
    pub fn get_cell_mut(&mut self, sigma: usize) -> LatticeEntity<&mut Cell> {
        match sigma {
            0 => LatticeEntity::Medium,
            _ => LatticeEntity::SomeCell(&mut self.cell_vec[sigma - 1])
        }
    }
    
    pub fn n_cells(&self) -> usize {
        self.cell_vec.len()
    }

    pub fn spawn_rect_cell(&mut self, rect: Rect<usize>, target_area: u32) -> Option<usize> {
        let mut cell_area = 0u32;
        let sigma = self.n_cells() + 1;
        for p in rect.iter_positions() {
            if self.cell_lattice[p] != 0 {
                continue;
            }
            self.cell_lattice[p] = sigma;
            for neigh in self.cell_lattice.validate(p.moore_neighs(self.neigh_r)) {
                let edge = Edge::new(p, neigh, self.neigh_r).unwrap();
                if self.cell_lattice[neigh] != sigma {
                    self.edge_bk.insert(edge);
                } else { 
                    self.edge_bk.remove(&edge);
                }
            }
            cell_area += 1;
        }
        if cell_area == 0 { 
            return None;
        }
        self.cell_vec.push(Cell::new(cell_area, target_area));
        Some(sigma)
    }
    
    pub fn update_edges(&mut self, pos: Pos2D<usize>) -> (u16, u16) {
        let mut removed = 0;
        let mut added = 0;
        let sigma = self.cell_lattice[pos];
        for neigh in self.cell_lattice.validate(pos.moore_neighs(self.neigh_r)) {
            let edge = Edge::new(pos, neigh, self.neigh_r).unwrap();
            let sigma_neigh = self.cell_lattice[neigh];
            if sigma == sigma_neigh {
                self.edge_bk.remove(&edge);
                removed += 1;
            } else if self.edge_bk.insert(edge) { 
                added += 1;
            }
        }
        (removed, added)
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
        assert_eq!(env.edge_bk.len(), 8 * 4 * 3 + 4 * 5);
        env.spawn_rect_cell(
            Rect::new(
                Pos2D::new(15, 15),
                Pos2D::new(25, 25)
            ),
            10
        );
        assert_eq!(env.get_cell(1).unwrap().area, 100);
        assert_eq!(env.get_cell(2).unwrap().area, 75);
    }
}
