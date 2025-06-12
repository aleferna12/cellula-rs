use std::cmp::min;
use indexmap::IndexSet;
use rand::Rng;
use crate::cell::Cell;
use crate::lattice::Lattice;
use crate::pos::{Edge, Pos2D, Rect};

enum NotCell {
    Medium
}

pub struct Environment {
    pub cell_lattice: Lattice<usize>,
    cell_vec: Vec<Cell>,
    // TODO: profile using this crate, I have no clue of whether it's fast enough
    edge_set: IndexSet<Edge>,
    pub neigh_r: u8
}
impl Environment {
    pub fn new(width: usize, height: usize, neigh_r: u8) -> Self {
        Self {
            cell_lattice: Lattice::new(width, height),
            cell_vec: vec![],
            edge_set: IndexSet::new(),
            neigh_r
        }
    }

    // Cell population functions
    pub fn get_cell(&self, sigma: usize) -> Option<&Cell> {
        if self.is_cell(sigma) { 
            Some(&self.cell_vec[sigma - 1] )
        } else { 
            None
        }
    }
    
    pub fn n_cells(&self) -> usize {
        self.cell_vec.len()
    }

    pub fn spawn_rect_cell(&mut self, rect: Rect<usize>, target_area: u32) -> Option<usize> {
        let mut cell_area = 0usize;
        let sigma = self.n_cells() + 1;
        for p in rect.iterate_pos() {
            if self.cell_lattice[p] != 0 {
                continue;
            }
            self.cell_lattice[p] = sigma;
            for neigh in self.cell_lattice.filter_inbounds(p.moore_neighs(self.neigh_r)) {
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
        self.cell_vec.push(Cell::new(cell_area as u32, target_area));
        Some(sigma)
    }
    
    pub fn is_cell(&self, sigma: usize) -> bool {
        sigma != NotCell::Medium as usize && sigma <= self.n_cells()
    }

    // Edge bookkeeping functions
    pub fn n_edges(&self) -> usize { self.edge_set.len() }

    pub fn insert_edge(&mut self, edge: Edge) -> bool {
        self.edge_set.insert(edge)
    }

    fn remove_edge(&mut self, edge: &Edge) -> bool {
        self.edge_set.swap_remove(edge)
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
                -min(dist, oldp.0)..min(dist + 1, self.cell_lattice.width() as i32 - oldp.0)
            );
            newp.1 = oldp.1 + rng.random_range(
                -min(dist, oldp.1)..min(dist + 1, self.cell_lattice.height() as i32 - oldp.1)
            );
        }
        Pos2D::new(newp.0 as usize, newp.1 as usize)
    }

    // TODO: think about whether these should be in a trait CA that we implement for environment
    // CA functions
    pub fn ca_step() {
        todo!()
    }

    pub fn delta_hamiltonian(&self, copy_attempt: &Edge, size_lambda: f32) -> f32 {
        let mut delta_h = 0f32;
        let sigma_from = self.cell_lattice[copy_attempt.p1];
        let sigma_into = self.cell_lattice[copy_attempt.p2];
        if let Some(cell) = self.get_cell(sigma_from) {
            delta_h += self.delta_hamiltonian_size(
                cell.area + 1,
                cell.target_area,
                size_lambda
            )
        }
        if let Some(cell) = self.get_cell(sigma_into) {
            delta_h += self.delta_hamiltonian_size(
                cell.area - 1,
                cell.target_area,
                size_lambda
            )
        }
        delta_h
    }

    pub fn delta_hamiltonian_size(&self, area: u32, target_area: u32, size_lambda: f32) -> f32 {
        let da = area.abs_diff(target_area) as f32;
        size_lambda * da * da
    }
}

#[cfg(test)]
mod tests {
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256StarStar;
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
    fn test_random_neighbour() {
        let env = empty_env();
        let mut rng = Xoshiro256StarStar::from_os_rng();
        for neigh_r in 1..4 {
            let mut too_far = false;
            for _ in 0..1000 {
                let p1 = env.cell_lattice.random_pos(&mut rng);
                let p2 = env.random_neighbour(&p1, neigh_r, &mut rng);
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
        let mut env = env_with_cell();
        assert_eq!(env.edge_set.len(), 8 * 4 * 3 + 4 * 5);
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

    #[test]
    fn test_delta_hamiltonian() {
        let env = env_with_cell();
        let cp_att = Edge::new((10, 10).into(), (9, 9).into(), 1).unwrap();
        let dh = env.delta_hamiltonian(&cp_att, 1f32);
        assert_eq!(dh, 1f32);
    }
}
