use crate::cell::{Cellular, RelCell};
use crate::constants::Spin;
use crate::lattice::Lattice;
use crate::positional::boundary::AsLatticeBoundary;
use crate::positional::neighbourhood::Neighbourhood;
use crate::positional::pos::Pos;
use crate::positional::rect::Rect;
use std::cmp::max;
use std::collections::{HashSet};
use std::error::Error;
use std::f32::consts::PI;

// TODO: generalise this
//  a lot of calls in env and cellular automata currently rely on this having a chem gradient layer
//  we can make a trait HasChem to make this general
//  there should also be a trait SpaceLike, that implements a method to shift the cell positions appropriately
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

    /// This is the fastest cell search function possible, but it is NOT SAFE.
    ///
    /// <div class="warning">
    ///
    /// This function should only be used when not all positions are required to be found.
    ///
    /// Prefer `search_cell_box()`, which warns about missing values.
    ///
    /// </div>
    pub fn search_cell_box_iter(
        &self,
        cell: &RelCell<impl Cellular>,
        diameter_scaler: f32
    ) -> impl Iterator<Item = Pos<usize>> {
        let search_diam = (
            diameter_scaler
                * 2.
                * (max(cell.target_area(), cell.area()) as f32 / PI)
                .sqrt()
        ) as usize;
        self.cell_lattice.search_box(
            &cell.spin,
            cell.center().to_usize(),
            search_diam,
            &self.lat_bound,
        )
    }

    // This function returns a Vec so that we can check that the site number matches
    /// Searches for all cell positions by creating a box around the cell and iterating all the positions inside it.
    ///
    /// May fail if `radius_scaler` is too small, in which case logs a warning.
    pub fn search_cell_box(&self, cell: &RelCell<impl Cellular>, diameter_scaler: f32) -> Vec<Pos<usize>> {
        let found: Vec<_> = self.search_cell_box_iter(cell, diameter_scaler).collect();
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
    /// Is considerably slower than `search_cell_box()` and may fail if the cell is not contiguous
    /// or if the cell centre is not a cell position.
    pub fn search_cell_contiguous<N: Neighbourhood>(
        &self,
        cell: &RelCell<impl Cellular>,
        neighbourhood: &N
    ) -> Vec<Pos<usize>> {
        let found = self.cell_lattice.search_contiguous(
            &cell.spin,
            cell.center().to_usize(),
            &self.lat_bound,
            neighbourhood
        );

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

    pub fn search_cell_outline<N: Neighbourhood>(
        &self,
        cell: &RelCell<impl Cellular>,
        diameter_scaler: f32,
        neighbourhood: &N
    ) -> Vec<Pos<usize>> {
        let search_diam = (
            diameter_scaler
                * 2.
                * (max(cell.target_area(), cell.area()) as f32 / PI)
                .sqrt()
        ) as usize;
        self.cell_lattice.search_outline(
            &cell.spin,
            cell.center().to_usize(),
            search_diam,
            &self.lat_bound,
            neighbourhood
        )
    }

    pub fn cell_neighbours<N: Neighbourhood>(
        &self,
        cell: &RelCell<impl Cellular>,
        diameter_scaler: f32,
        neighbourhood: &N
    ) -> HashSet<Spin> {
        let outline = self.search_cell_outline(
            cell,
            diameter_scaler,
            neighbourhood,
        ).into_iter().map(|pos| { self.cell_lattice[pos] });
        HashSet::from_iter(outline)
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

    fn space_cell_pair(positions: &[Pos<usize>]) -> (Space<BoundaryType>, RelCell<impl Cellular>) {
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
        let boxed_positions = space.search_cell_box(&cell, 2.0);
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
        let contiguous_positions = space.search_cell_contiguous(&cell, &neighbourhood);

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
        let result = space.search_cell_contiguous(&cell, &neighbourhood);

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