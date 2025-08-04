use crate::cell::{CellLike, RelCell};
use crate::constants::Spin;
use crate::lattice::Lattice;
use crate::positional::boundary::{AsLatticeBoundary, Boundary};
use crate::positional::neighbourhood::Neighbourhood;
use crate::positional::pos::Pos;
use crate::positional::rect::Rect;
use std::cmp::max;
use std::collections::{HashSet, VecDeque};
use std::error::Error;
use std::f32::consts::PI;

pub struct Space<B: AsLatticeBoundary> {
    pub bound: B,
    pub lat_bound: B::LatticeBoundary,
    pub cell_lattice: Lattice<Spin>,
    pub chem_lattice: Lattice<u32>,
}

impl<B: AsLatticeBoundary> Space<B> {
    pub fn new(bound: B) -> Result<Self, Box<dyn Error>> 
    where 
        B: AsLatticeBoundary<Coord = f32>,
        B::Error: 'static + Error {
        let rect: Rect<usize> = bound.rect().clone().try_into()?;
        Ok(Self{
            lat_bound: bound.as_lattice_boundary()?,
            cell_lattice: Lattice::<Spin>::new(rect.clone()),
            chem_lattice: Lattice::<u32>::new(rect),
            bound,
        })
    }

    // TODO!: These should most definitely be implemented as part of a Lattice<impl PartialEq> API
    //  or moved into Environment such that we dont need to pass in the Neighbourhood
    //  (but this would be painful due to borrow semantics)
    /// This is the fastest cell search function possible, but it is NOT SAFE.
    ///
    /// Prefer `box_cell_positions()`, which warns about missing values.
    /// This function should only be used when not all positions are required to be found.
    pub fn iter_box_cell_positions(
        &self, 
        cell: &RelCell<impl CellLike>,
        radius_scaler: f32
    ) -> impl Iterator<Item = Pos<usize>> {
        let search_radius = (radius_scaler * (max(cell.target_area(), cell.area()) as f32 / PI).sqrt()) as isize;
        let center = Pos::new(
            cell.center().x as isize,
            cell.center().y as isize
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
    pub fn box_cell_positions(&self, cell: &RelCell<impl CellLike>, radius_scaler: f32) -> Vec<Pos<usize>> {
        let found: Vec<_> = self.iter_box_cell_positions(cell, radius_scaler).collect();
        if found.len() != cell.area() as usize {
            log::warn!(
                "Only found {} positions out of the {} expected for cell with spin {} \
                (try to increase `search-radius`)", 
                found.len(),
                cell.area(),
                cell.spin
            )
        }
        found
    }

    /// Searches for all cell positions with a BFS algorithm to traverse the lattice sites.
    ///
    /// Is considerably slower than `box_cell_positions()` and may fail if cell is not contiguous 
    /// or if the cell center is not a cell position.
    pub fn contiguous_cell_positions<N: Neighbourhood>(
        &self,
        cell: &RelCell<impl CellLike>, 
        neighbourhood: &N
    ) -> Vec<Pos<usize>> {
        let mut visited = Lattice::<bool>::new(self.cell_lattice.rect.clone());
        let mut found = Vec::with_capacity(cell.area() as usize);
        let mut queue = VecDeque::from([Pos::new(
            cell.center().x as isize,
            cell.center().y as isize
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

        if found.len() != cell.area() as usize {
            log::warn!(
                "Only found {} positions out of the {} expected for cell with spin {} \
                (cell might be discontiguous)", 
                found.len(),
                cell.area(),
                cell.spin
            )
        }
        found
    }

    pub fn cell_neighbours<N: Neighbourhood>(
        &self,
        cell: &RelCell<impl CellLike>,
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cell::{Cell, RelCell};
    use crate::constants::BoundaryType;
    use crate::genome::MockGenome;
    use crate::positional::neighbourhood::MooreNeighbourhood;
    use crate::positional::pos::Pos;

    fn space_cell_pair(positions: &[Pos<usize>]) -> (Space<BoundaryType>, RelCell<impl CellLike>) {
        let mut cell = RelCell::mock(Cell::new_empty(
            10,
            20,
            MockGenome::new(0)
        ));
        let mut space = Space::new(BoundaryType::new(Rect::new(
            (0., 0.).into(),
            (10., 10.).into()
        ))).unwrap();
        for pos in positions {
            space.cell_lattice[*pos] = cell.spin;
            cell.shift_position(*pos, true, &space.bound)
        }
        (space, cell)
    }

    #[test]
    fn test_box_cell_positions() {
        let positions = [
            Pos::new(5, 5),
            Pos::new(5, 6),
            Pos::new(6, 5),
            Pos::new(6, 6),
        ];
        let (space, cell) = space_cell_pair(&positions);
        let boxed_positions = space.box_cell_positions(&cell, 2.0);
        assert_eq!(boxed_positions.len(), positions.len());
        for pos in &positions {
            assert!(boxed_positions.contains(pos));
        }
    }

    #[test]
    fn test_contiguous_cell_positions() {
        let positions = [
            Pos::new(5, 5),
            Pos::new(5, 6),
            Pos::new(6, 5),
            Pos::new(6, 6),
        ];
        let (space, cell) = space_cell_pair(&positions);
        let neighbourhood = MooreNeighbourhood::new(1);
        let contiguous_positions = space.contiguous_cell_positions(&cell, &neighbourhood);

        // Should find all 4 contiguous positions
        assert_eq!(contiguous_positions.len(), positions.len());
        for pos in &positions {
            assert!(contiguous_positions.contains(pos));
        }
    }

    #[test]
    fn test_contiguous_cell_positions_discontiguous() {
        let positions = [
            Pos::new(5, 5),
            Pos::new(5, 6),
            Pos::new(7, 7), // discontiguous point
        ];
        let (space, cell) = space_cell_pair(&positions);
        let neighbourhood = MooreNeighbourhood::new(1);
        let result = space.contiguous_cell_positions(&cell, &neighbourhood);

        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_cell_neighbours() {
        let positions = [
            Pos::new(2, 2),
            Pos::new(2, 1),
        ];
        let (mut space, cell) = space_cell_pair(&positions);
        let neighbourhood = MooreNeighbourhood::new(1);

        let neighbour_spins = [cell.spin + 1, cell.spin + 2];
        space.cell_lattice[Pos::new(1, 2)] = neighbour_spins[0];
        space.cell_lattice[Pos::new(2, 0)] = neighbour_spins[1];

        let neighs = space.cell_neighbours(&cell, 1.0, &neighbourhood);

        assert!(neighs.contains(&neighbour_spins[0]));
        assert!(neighs.contains(&neighbour_spins[1]));
        assert!(!neighs.contains(&cell.spin));
    }
}