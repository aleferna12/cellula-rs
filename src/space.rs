use crate::cell::RelCell;
use crate::constants::{BoundaryType, LatticeBoundaryType, Spin};
use crate::lattice::Lattice;
use crate::positional::boundary::Boundary;
use crate::positional::neighbourhood::Neighbourhood;
use crate::positional::pos::Pos;
use crate::positional::rect::Rect;
use std::cmp::max;
use std::collections::{HashSet, VecDeque};
use std::f32::consts::PI;

pub struct Space {
    pub bound: BoundaryType,
    pub lat_bound: LatticeBoundaryType,
    pub cell_lattice: Lattice<Spin>,
    pub light_lattice: Lattice<usize>,
}

impl Space {
    pub fn new(width: usize, height: usize) -> Self {
        let rect = Rect::new(
            (0, 0).into(),
            (width, height).into()
        );
        Self {
            bound: BoundaryType::new(Rect::new(
                (0., 0.).into(),
                (width as f32, height as f32).into()
            )),
            lat_bound: LatticeBoundaryType::new(Rect::new(
                (0, 0).into(),
                (width as isize, height as isize).into()
            )),
            cell_lattice: Lattice::new(rect.clone()),
            light_lattice: Lattice::new(rect),
        }
    }

    /// This is the fastest cell search function possible, but it is NOT SAFE.
    ///
    /// Prefer `box_cell_positions()`, which warns about missing values.
    /// This function should only be used when not all positions are required to be found.
    pub fn iter_box_cell_positions(&self, cell: &RelCell, radius_scaler: f32) -> impl Iterator<Item = Pos<usize>> {
        let search_radius = (radius_scaler * (max(cell.target_area, cell.area) as f32 / PI).sqrt()) as isize;
        let center = Pos::new(
            cell.center.pos.x as isize,
            cell.center.pos.y as isize
        );
        let rect = Rect::new(
            (center.x - search_radius, center.y - search_radius).into(),
            (center.x + search_radius, center.y + search_radius).into(),
        );
        self.lat_bound
            .valid_positions(rect.iter_positions())
            .filter_map(|pos| {
                let lat_pos = Pos::<usize>::from(pos);
                if self.cell_lattice[lat_pos] == cell.spin {
                    return Some(lat_pos);
                }
                None
            })
    }

    // This function returns a Vec so that we can check that the site number matches
    /// Searches for all cell positions by creating a box around the cell and iterating all the positions inside of it.
    ///
    /// May fail if `radius_scaler` is too small.
    pub fn box_cell_positions(&self, cell: &RelCell, radius_scaler: f32) -> Vec<Pos<usize>> {
        let found: Vec<_> = self.iter_box_cell_positions(cell, radius_scaler).collect();
        if found.len() != cell.area as usize {
            log::warn!(
                "Only found {} positions out of the {} expected for cell with spin {} \
                (try to increase `search-radius`)", 
                found.len(),
                cell.area,
                cell.spin
            )
        }
        found
    }

    /// Searches for all cell positions with a BFS algorithm to traverse the lattice sites.
    ///
    /// Is considerably slower than `box_cell_positions()` and may fail if cell is not contiguous 
    /// or if the cell center is not a cell position.
    pub fn contiguous_cell_positions<N: Neighbourhood>(&self, cell: &RelCell, neighbourhood: &N) -> Vec<Pos<usize>> {
        let mut visited = Lattice::<bool>::new(self.cell_lattice.rect.clone());
        let mut found = Vec::with_capacity(cell.area as usize);
        let mut queue = VecDeque::from([Pos::new(
            cell.center.pos.x as isize,
            cell.center.pos.y as isize
        )]);

        while let Some(pos) = queue.pop_front() {
            let lat_pos = Pos::from(pos);
            if cell.spin != self.cell_lattice[lat_pos] {
                continue;
            }
            self.lat_bound
                .valid_positions(neighbourhood.neighbours(pos))
                .for_each(|neigh| {
                    let lat_neigh = Pos::from(neigh);
                    if !visited[lat_neigh] {
                        visited[lat_neigh] = true;
                        queue.push_back(neigh);
                    }
                });
            visited[lat_pos] = true;
            found.push(lat_pos);
        }

        if found.len() != cell.area as usize {
            log::warn!(
                "Only found {} positions out of the {} expected for cell with spin {} \
                (cell might be discontiguous)", 
                found.len(),
                cell.area,
                cell.spin
            )
        }
        found
    }

    pub fn cell_neighbours<N: Neighbourhood>(
        &self,
        cell: &RelCell,
        radius_scaler: f32,
        neighbourhood: &N
    ) -> HashSet<Spin> {
        let mut neighs = HashSet::default();
        let mut border_pos = None;
        for pos in self.iter_box_cell_positions(cell, radius_scaler) {
            if let Some(neigh) = self.lat_bound.valid_pos(Pos::new(pos.x as isize - 1, pos.y as isize)) {
                if self.cell_lattice[pos] != self.cell_lattice[Pos::from(neigh)] {
                    border_pos = Some(pos);
                    break
                }
            }
        }
        if border_pos.is_none() {
            return neighs;
        }

        let mut visited = Lattice::<bool>::new(self.cell_lattice.rect.clone());
        let mut queue = VecDeque::from([border_pos.unwrap().into()]);
        while let Some(pos) = queue.pop_front() {
            let spin = self.cell_lattice[Pos::from(pos)];
            let mut same_spin_neighs = Vec::new();
            let mut has_diff_neighbor = false;
            for neigh in self.lat_bound.valid_positions(neighbourhood.neighbours(pos)) {
                let neigh_pos = Pos::from(neigh);

                let neigh_spin = self.cell_lattice[neigh_pos];
                if neigh_spin != spin {
                    has_diff_neighbor = true;
                    neighs.insert(neigh_spin);
                } else if !visited[neigh_pos] {
                    visited[neigh_pos] = true;
                    same_spin_neighs.push(neigh);
                }
            }

            if has_diff_neighbor {
                queue.extend(same_spin_neighs);
            }
        }
        neighs
    }
}