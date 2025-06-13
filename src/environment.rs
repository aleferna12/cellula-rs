use crate::cell::Cell;
use crate::edge::{Edge, EdgeBookkeeper};
use crate::lattice::Lattice;
use crate::pos::{Pos2D, Rect};
use rand::Rng;
use std::f32::consts::E;

enum NotCell {
    Medium
}

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

    // Cell population functions
    pub fn get_cell(&self, sigma: usize) -> Option<&Cell> {
        if self.is_cell(sigma) {
            Some(&self.cell_vec[sigma - 1] )
        } else {
            None
        }
    }
    
    pub fn get_cell_mut(&mut self, sigma: usize) -> Option<&mut Cell> {
        if self.is_cell(sigma) { 
            Some(&mut self.cell_vec[sigma - 1] )
        } else { 
            None
        }
    }
    
    pub fn n_cells(&self) -> usize {
        self.cell_vec.len()
    }

    pub fn spawn_rect_cell(&mut self, rect: Rect<usize>, target_area: u32) -> Option<usize> {
        let mut cell_area = 0u32;
        let sigma = self.n_cells() + 1;
        for p in rect.iterate_pos() {
            if self.cell_lattice[p] != 0 {
                continue;
            }
            self.cell_lattice[p] = sigma;
            for neigh in self.cell_lattice.filter_inbounds(p.moore_neighs(self.neigh_r)) {
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
    
    pub fn is_cell(&self, sigma: usize) -> bool {
        sigma != NotCell::Medium as usize && sigma <= self.n_cells()
    }

    // TODO: think about whether these should be in a trait CA that we implement for environment (or model)
    //   Or: make a CA struct that has all necessary parameters and store it in model
    //   Also write tests
    // CA functions
    pub fn ca_step(&mut self, rng: &mut impl Rng, size_lambda: f32, boltz_t: f32) {
        // TODO: ensure this makes sense for neigh_r > 1
        let edge_per_pos = self.neigh_r as f32 / 2f32;
        let mut to_visit = self.edge_bk.len() as f32 / edge_per_pos;
        while 0f32 < to_visit {
            let edge_i = self.edge_bk.random_index(rng);
            let edge = self.edge_bk.at(edge_i);
            // TODO: is this really faster than just keeping both edges in the IndexSet? Benchmark
            let (p1, p2) = if rng.random::<f32>() < 0.5 {
                (edge.p1, edge.p2)
            } else {
                (edge.p2, edge.p1)
            };
            let sigma_from = self.cell_lattice[p1];
            let sigma_to = self.cell_lattice[p2];
            let delta_h = self.delta_hamiltonian(sigma_from, sigma_to, size_lambda);
            if Environment::accept_copy(rng, delta_h, boltz_t) {
                self.cell_lattice[p2] = sigma_from;
                if let Some(cell) = self.get_cell_mut(sigma_from) {
                    cell.area += 1;
                }
                if let Some(cell) = self.get_cell_mut(sigma_to) {
                    cell.area -= 1;
                }
                let (removed, added) = self.update_edges(p2);
                // TODO: ensure this makes sense for neigh_r > 1
                to_visit += (added as f32 - removed as f32) / edge_per_pos;
            }
            to_visit -= 1f32;
        }
    }
    
    pub fn update_edges(&mut self, pos: Pos2D<usize>) -> (u16, u16) {
        let mut removed = 0;
        let mut added = 0;
        let sigma = self.cell_lattice[pos];
        for neigh in self.cell_lattice.filter_inbounds(pos.moore_neighs(self.neigh_r)) {
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
    
    pub fn accept_copy(rng: &mut impl Rng, delta_h: f32, boltz_t: f32) -> bool {
        delta_h < 0f32 || rng.random::<f32>() < E.powf(-delta_h / boltz_t)
    }

    pub fn delta_hamiltonian(&self, sigma_from: usize, sigma_to: usize, size_lambda: f32) -> f32 {
        let mut delta_h = 0f32;
        if let Some(cell) = self.get_cell(sigma_from) {
            delta_h += self.delta_hamiltonian_size(
                1,
                cell.area,
                cell.target_area,
                size_lambda
            )
        }
        if let Some(cell) = self.get_cell(sigma_to) {
            delta_h += self.delta_hamiltonian_size(
                -1,
                cell.area,
                cell.target_area,
                size_lambda
            )
        }
        delta_h
    }

    pub fn delta_hamiltonian_size(&self, delta_area: i32, area: u32, target_area: u32, size_lambda: f32) -> f32 {
        2f32 * size_lambda * delta_area as f32 * (area as f32 - target_area as f32) + size_lambda
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

    #[test]
    fn test_delta_hamiltonian() {
        let env = env_with_cell();
        let cp_att = Edge::new((10, 10).into(), (9, 9).into(), 1).unwrap();
        let sigma_from = env.cell_lattice[cp_att.p1];
        let sigma_to = env.cell_lattice[cp_att.p2];
        let dh = env.delta_hamiltonian(sigma_from, sigma_to, 1f32);
        assert_eq!(dh, 1f32);
    }
}
