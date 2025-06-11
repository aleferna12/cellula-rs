use std::cmp::min;
use indexmap::IndexSet;
use rand::Rng;
use crate::cell::Cell;
use crate::lattice::Lattice;
use crate::pos::{Edge, Pos2D, Rect};

pub struct Dish {
    pub cell_lattice: Lattice<usize>,
    cell_vec: Vec<Cell>,
    // TODO: profile using this struct, I have no clue of whether it's fast enough
    edge_set: IndexSet<Edge>,
    pub neigh_r: u8
}
impl Dish {
    pub fn new(width: usize, height: usize, neigh_r: u8) -> Self {
        Self {
            cell_lattice: Lattice::new(width, height),
            cell_vec: vec![],
            edge_set: IndexSet::new(),
            neigh_r
        }
    }
    
    pub fn get_cell(&self, sigma: usize) -> &Cell {
        &self.cell_vec[sigma - 1]
    }
    
    pub fn n_cells(&self) -> usize {
        self.cell_vec.len()
    }

    pub fn spawn_rect_cell(&mut self, rect: Rect<usize>) -> Option<usize> {
        let mut cell_area = 0usize;
        let sigma = self.n_cells() + 1;
        for p in rect.iterate_pos() {
            if self.cell_lattice[p] != 0 {
                continue;
            }
            self.cell_lattice[p] = sigma;
            for neigh in self.cell_lattice.moore_neighs(&p, self.neigh_r) {
                let edge = Edge::new(p, neigh, self.neigh_r).unwrap();
                if self.cell_lattice[neigh] != sigma {
                    self.insert_edge(edge);
                } else { 
                    self.remove_edge(&edge);
                }
            }
            cell_area += 1;
        }
        if cell_area == 0 { 
            return None;
        }
        self.cell_vec.push(Cell::new(cell_area as u32));
        Some(sigma)
    }

    pub fn n_edges(&self) -> usize { self.edge_set.len() }

    pub fn insert_edge(&mut self, edge: Edge) -> bool {
        self.edge_set.insert(edge)
    }

    fn remove_edge(&mut self, edge: &Edge) {
        self.edge_set.swap_remove(edge);
    }

    pub fn remove_random_edge(&mut self, rng: &mut impl Rng) -> Edge {
        let index = rng.random_range(0..self.edge_set.len() - 1);
        self.edge_set.swap_remove_index(index).unwrap()
    }

    pub fn random_neighbour(&self, p: &Pos2D<usize>, neigh_r: u8, rng: &mut impl Rng) -> Pos2D<usize> {
        let oldp = (p.x as i32, p.y as i32);
        let mut newp = oldp;
        let dist = neigh_r as i32;
        while oldp == newp {
            newp.0 = oldp.0 + rng.random_range(
                -min(dist, oldp.0)..min(dist + 1, self.cell_lattice.width as i32 - oldp.0)
            );
            newp.1 = oldp.1 + rng.random_range(
                -min(dist, oldp.1)..min(dist + 1, self.cell_lattice.height as i32 - oldp.1)
            );
        }
        Pos2D::new(newp.0 as usize, newp.1 as usize)
    }
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256StarStar;
    use super::*;

    #[test]
    fn test_random_neighbour() {
        let dish = Dish::new(100, 100, 1);
        let mut rng = Xoshiro256StarStar::from_os_rng();
        for neigh_r in 1..4 {
            let mut too_far = false;
            for _ in 0..1000 {
                let p1 = dish.cell_lattice.random_pos(&mut rng);
                let p2 = dish.random_neighbour(&p1, neigh_r, &mut rng);
                assert!(Edge::new(p1, p2, neigh_r).is_ok());
                if !too_far {
                    too_far = Edge::new(p1, p2, neigh_r - 1).is_err()
                }
            }
            assert!(too_far)
        }
    }
    
    #[test]
    fn test_spawn_rect_cell() {
        let mut dish = Dish::new(100, 100, 1);
        let sigma1 = dish.spawn_rect_cell(
            Rect::new(
                Pos2D::new(10, 10), 
                Pos2D::new(20, 20)
            )
        );
        assert_eq!(dish.edge_set.len(), 8 * 4 * 3 + 4 * 5);
        let sigma2 = dish.spawn_rect_cell(
            Rect::new(
                Pos2D::new(15, 15),
                Pos2D::new(25, 25)
            )
        );
        assert_eq!(dish.get_cell(sigma1.unwrap()).area, 100);
        assert_eq!(dish.get_cell(sigma2.unwrap()).area, 75);
    }
}