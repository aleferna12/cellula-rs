use std::cmp::max;
use std::collections::{HashSet, VecDeque};
use std::f32::consts::PI;
use rand::Rng;
use std::ops::{Deref, DerefMut, Index, IndexMut};
use crate::cell::RelCell;
use crate::constants::Spin;
use crate::positional::boundary::LatticeBoundary;
use crate::positional::neighbourhood::Neighbourhood;
use crate::positional::pos::Pos2D;
use crate::positional::rect::Rect;

pub struct Lattice<T, B> {
    array: Box<[T]>,
    pub bound: B
}
// Since the lattice size is naturally usize, boundary coord should be isize to avoid overflow errors
// Although technically it only has to be slightly larger than its defined size
impl<T: Default + Copy, B: LatticeBoundary> Lattice<T, B> {
    pub fn new(bound: B) -> Self {
        Self {
            array: vec![T::default(); bound.rect().width() as usize * bound.rect().height() as usize]
                .into_boxed_slice(),
            bound,
        }
    }

    pub fn width(&self) -> usize {
        self.bound.rect().width() as usize
    }

    pub fn height(&self) -> usize {
        self.bound.rect().height() as usize
    }

    pub fn random_pos(&self, rng: &mut impl Rng) -> Pos2D<usize> {
        Pos2D::new(
            rng.random_range(0..self.width()),
            rng.random_range(0..self.height())
        )
    }

    pub fn iter_positions(&self) -> impl Iterator<Item = Pos2D<usize>> {
        self.bound.rect().iter_positions().map(|p| Pos2D::new(
            p.x as usize,
            p.y as usize
        ))
    }

    pub fn iter_values(&self) -> impl Iterator<Item = T> {
        self.iter_positions()
            .map(|pos| {
                self[pos]
            })
    }
}

impl<T: Copy + Default, B: LatticeBoundary> Index<Pos2D<usize>> for Lattice<T, B> {
    type Output = T;

    fn index(&self, pos: Pos2D<usize>) -> &Self::Output {
        &self.array[pos.row_major(self.height())]
    }
}

impl<T: Copy + Default, B: LatticeBoundary> IndexMut<Pos2D<usize>> for Lattice<T, B> {
    fn index_mut(&mut self, pos: Pos2D<usize>) -> &mut Self::Output {
        &mut self.array[pos.row_major(self.height())]
    }
}

pub struct CellLattice<B> {
    pub lat: Lattice<Spin, B>
}

impl<B: LatticeBoundary + Clone> CellLattice<B> {
    pub fn new(bound: B) -> Self {
        Self { lat: Lattice::new(bound) }
    }
    
    /// This is the fastest cell search function possible, but it is NOT SAFE.
    /// 
    /// Prefer `box_cell_positions()`, which warns about missing values.
    /// This function should only be used when not all positions are required to be found.
    pub fn iter_box_cell_positions(&self, cell: &RelCell, radius_scaler: f32) -> impl Iterator<Item = Pos2D<usize>> {
        let search_radius = (radius_scaler * (max(cell.target_area, cell.area) as f32 / PI).sqrt()) as isize;
        let center = Pos2D::new(
            cell.center.pos.x as isize,
            cell.center.pos.y as isize
        );
        let rect = Rect::new(
            (center.x - search_radius, center.y - search_radius).into(),
            (center.x + search_radius, center.y + search_radius).into(),
        );
        self.bound
            .valid_positions(rect.iter_positions())
            .filter_map(|pos| {
                let lat_pos = Pos2D::<usize>::from(pos);
                if self[lat_pos] == cell.spin {
                    return Some(lat_pos);
                }
                None
            })
    }
    
    // This function returns a Vec so that we can check that the site number matches
    /// Searches for all cell positions by creating a box around the cell and iterating all the positions inside of it.
    ///
    /// May fail if `radius_scaler` is too small.
    pub fn box_cell_positions(&self, cell: &RelCell, radius_scaler: f32) -> Vec<Pos2D<usize>> {
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
    pub fn contiguous_cell_positions<N: Neighbourhood>(&self, cell: &RelCell, neighbourhood: &N) -> Vec<Pos2D<usize>> {
        let mut visited = Lattice::<bool, _>::new(self.bound.clone());
        let mut found = Vec::with_capacity(cell.area as usize);
        let mut queue = VecDeque::from([Pos2D::new(
            cell.center.pos.x as isize,
            cell.center.pos.y as isize
        )]);

        while let Some(pos) = queue.pop_front() {
            let lat_pos = Pos2D::from(pos);
            if cell.spin != self[lat_pos] {
                continue;
            }
            self.bound
                .valid_positions(neighbourhood.neighbours(pos))
                .for_each(|neigh| {
                    let lat_neigh = Pos2D::from(neigh);
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
            if let Some(neigh) = self.bound.valid_pos(Pos2D::new(pos.x as isize - 1, pos.y as isize)) {
                if self[pos] != self[Pos2D::from(neigh)] {
                    border_pos = Some(pos);
                    break
                }
            }
        }
        if border_pos.is_none() {
            return neighs;
        }

        let mut visited = Lattice::<bool, _>::new(self.bound.clone());
        let mut queue = VecDeque::from([border_pos.unwrap().into()]);
        while let Some(pos) = queue.pop_front() {
            let spin = self[Pos2D::from(pos)];
            let mut same_spin_neighs = Vec::new();
            let mut has_diff_neighbor = false;
            for neigh in self.bound.valid_positions(neighbourhood.neighbours(pos)) {
                let neigh_pos = Pos2D::from(neigh);

                let neigh_spin = self[neigh_pos];
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

impl<B> Deref for CellLattice<B> {
    type Target = Lattice<Spin, B>;

    fn deref(&self) -> &Self::Target {
        &self.lat
    }
}

impl<B> DerefMut for CellLattice<B> {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.lat
    }
}
