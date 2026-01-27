//! Contains logic associated with [`Environment`].

use crate::base::cell::Cell;
use crate::cell_container::{CellContainer, RelCell};
use crate::constants::FloatType;
use crate::lattice::Lattice;
use crate::positional::boundaries::{Boundaries, Boundary, ToLatticeBoundary};
use crate::positional::edge::Edge;
use crate::positional::edge_book::EdgeBook;
use crate::positional::neighborhood::Neighborhood;
use crate::positional::pos::{CastCoords, Pos};
use crate::prelude::{MooreNeighborhood, UnsafePeriodicBoundary};
use crate::spin::Spin;
use crate::traits::cellular::{Alive, Cellular, HasCenter};
use crate::traits::habitable::Habitable;
use core::fmt;
use std::cmp::max;
use std::collections::HashSet;
#[cfg(not(feature = "high-precision"))]
use std::f32::consts::PI;
#[cfg(feature = "high-precision")]
use std::f64::consts::PI;

// Has manual implementations for PartialEq, Debug and Clone (needed due to ToLatticeBoundary)
// If adding fields, remember to also change those!!!
/// An environment where cells are spatially and relationally localised.
pub struct Environment<C, N = MooreNeighborhood, B: ToLatticeBoundary = UnsafePeriodicBoundary<FloatType>> {
    /// Boundaries of the environment.
    ///
    /// These are used to validate that a position can be used to access information in the environment.
    pub bounds: Boundaries<B>,
    /// Cell container with spins matching those in [`Environment::cell_lattice`].
    pub cells: CellContainer<C>,
    /// Cell lattice with spins matching those in [`Environment::cells`].
    pub cell_lattice: Lattice<Spin>,
    /// Edge book containing all positions at cell-cell interfaces.
    pub edge_book: EdgeBook,
    /// Neighborhood of the environment.
    pub neighborhood: N
}

impl<C, N, B: ToLatticeBoundary<Coord = FloatType>> Environment<C, N, B> {
    /// Makes a new empty environment with no cells.
    pub fn new_empty(
        neighborhood: N,
        bounds: Boundaries<B>,
    ) -> Self {
        Self {
            cell_lattice: Lattice::from(bounds.boundary.rect().cast_coords()),
            edge_book: EdgeBook::new(),
            cells: CellContainer::new(),
            bounds,
            neighborhood
        }
    }
}

impl<C, N, B: ToLatticeBoundary> Environment<C, N, B> {
    /// Makes a new environment from its components.
    pub fn new(
        cells: CellContainer<C>,
        cell_lattice: Lattice<Spin>,
        neighborhood: N,
        bounds: Boundaries<B>,
    ) -> Self {
        Self {
            cell_lattice,
            bounds,
            cells,
            neighborhood,
            edge_book: EdgeBook::new()
        }
    }

    /// Returns the width of the environment.
    pub fn width(&self) -> usize {
        self.cell_lattice.width()
    }

    /// Returns the height of the environment.
    pub fn height(&self) -> usize {
        self.cell_lattice.height()
    }

    /// Returns an iterator over all sites in the neighborhood of `pos` that are within lattice boundaries.
    pub fn valid_neighbors(&self, pos: Pos<usize>) -> impl Iterator<Item = Pos<usize>>
    where N: Neighborhood {
        valid_neighbors(&self.bounds.lattice_boundary, &self.neighborhood, pos)
    }
}

impl<C: Cellular + HasCenter, N: Neighborhood, B: ToLatticeBoundary> Environment<C, N, B> {
    // This function returns a Box (not an Iter) so that we can check that the site number matches
    /// Searches for all cell positions by creating a box around the cell and iterating all the positions inside it.
    ///
    /// `search_scaler` multiplies the expected radius of `rel_cell` (which is calculated from its area).
    ///
    /// May fail to find all cell positions if `search_scaler` is too small, in which case logs a warning.
    pub fn search_cell_box(&self, rel_cell: &RelCell<C>, search_scaler: FloatType) -> Box<[Pos<usize>]> {
        let search_diam = (
            search_scaler
                * 2.
                * (max(rel_cell.cell.target_area(), rel_cell.cell.area()) as FloatType / PI)
                .sqrt()
        ) as usize;

        let found: Box<_> = self.cell_lattice.search_box(
            &Spin::Some(rel_cell.index),
            rel_cell.cell.center().round().cast_as(),
            search_diam,
            &self.bounds.lattice_boundary,
        ).collect();

        if found.len() != rel_cell.cell.area() as usize {
            log::warn!(
                "Only found {} positions out of the {} expected for cell with index {} \
                (try to increase `search-radius`)",
                found.len(),
                rel_cell.cell.area(),
                rel_cell.index
            )
        }
        found
    }

    /// Searches for all cell positions with a BFS algorithm to traverse the lattice sites.
    ///
    /// Is considerably slower than [`Environment::search_cell_box()`] and may not return all positions
    /// if the cell is not contiguous or if the cell center is not a cell position.
    pub fn search_cell_contiguous(
        &self,
        rel_cell: &RelCell<C>,
    ) -> Box<[Pos<usize>]> {
        let found = self.cell_lattice.search_contiguous(
            &Spin::Some(rel_cell.index),
            rel_cell.cell.center().round().cast_as(),
            &self.bounds.lattice_boundary,
            &self.neighborhood
        );

        if found.len() != rel_cell.cell.area() as usize {
            log::warn!(
                "Only found {} positions out of the {} expected for cell with index {} \
                (cell might be discontiguous)",
                found.len(),
                rel_cell.cell.area(),
                rel_cell.index
            )
        }
        found
    }

    /// Searches for all cell positions that interface a different spin.
    pub fn search_cell_outline(
        &self,
        rel_cell: &RelCell<C>,
        search_scaler: FloatType
    ) -> Box<[Pos<usize>]> {
        let search_diam = (
            search_scaler
                * 2.
                * (max(rel_cell.cell.target_area(), rel_cell.cell.area()) as FloatType / PI)
                .sqrt()
        ) as usize;

        self.cell_lattice.search_outline(
            &Spin::Some(rel_cell.index),
            rel_cell.cell.center().round().cast_as(),
            search_diam,
            &self.bounds.lattice_boundary,
            &self.neighborhood
        )
    }

    /// Finds all cells adjacent to `cell` using a search algorithm.
    pub fn cell_neighbors(
        &self,
        cell: &RelCell<C>,
        search_scaler: FloatType
    ) -> HashSet<Spin> {
        let outline = self
            .search_cell_outline(cell, search_scaler)
            .into_iter()
            .map(|pos| self.cell_lattice[pos]);
        HashSet::from_iter(outline)
    }

    /// Removes all cells from the environment and returns it to a clean state.
    pub fn wipe_out(&mut self) {
        self.cells.wipe_out();
        self.cell_lattice.iter_values_mut().for_each(|value| {
            if let Spin::Some(_) = value {
                *value = Spin::Medium;
            }
        });
        self.edge_book.clear();
    }

    /// Updates the edges around the position `pos` and counts how many were added/removed.
    pub fn update_edges(&mut self, pos: Pos<usize>) -> EdgesUpdate {
        let mut removed = 0;
        let mut added = 0;
        let spin = self.cell_lattice[pos];
        let valid_neighs = valid_neighbors(
            &self.bounds.lattice_boundary,
            &self.neighborhood,
            pos
        );
        for neigh in valid_neighs {
            let edge = Edge::new(pos, neigh);
            let spin_neigh = self.cell_lattice[neigh];
            // The order of these if statements matter A LOT, dont mess with it
            if spin == spin_neigh {
                if self.edge_book.remove(&edge) {
                    removed += 1;
                }
                continue;
            }
            if (matches!(spin, Spin::Some(_))
                || matches!(spin_neigh, Spin::Some(_)))
                && self.edge_book.insert(edge) {
                added += 1;
            }
        }
        EdgesUpdate { added, removed }
    }
}

impl<N: Neighborhood, B: ToLatticeBoundary<Coord = FloatType>> Habitable for Environment<Cell, N, B> {
    type Cell = Cell;

    fn env(&self) -> &Environment<Self::Cell, impl Neighborhood, impl ToLatticeBoundary> {
        self
    }

    fn env_mut(&mut self) -> &mut Environment<Self::Cell, impl Neighborhood, impl ToLatticeBoundary> {
        self
    }

    fn grant_position(
        &mut self,
        pos: Pos<usize>,
        to: Spin
    ) -> EdgesUpdate {
        if let Spin::Some(index) = to {
            self.cells[index].cell.shift_position(pos, true, &self.bounds.boundary);
        }
        if let Spin::Some(index) = self.cell_lattice[pos] {
            let from_rel_cell = &mut self.cells[index];
            from_rel_cell.cell.shift_position(pos, false, &self.bounds.boundary);
            if from_rel_cell.cell.area() == 0 {
                from_rel_cell.cell.apoptosis();
            }
        }
        // Executes the copy
        self.cell_lattice[pos] = to;
        self.update_edges(pos)
    }
}

impl<C, N, B: ToLatticeBoundary> Clone for Environment<C, N, B>
where
    C: Clone,
    N: Clone,
    B: Clone,
    B::LatticeBoundary: Clone {
    fn clone(&self) -> Self {
        Environment {
            bounds: self.bounds.clone(),
            cells: self.cells.clone(),
            cell_lattice: self.cell_lattice.clone(),
            edge_book: self.edge_book.clone(),
            neighborhood: self.neighborhood.clone(),
        }
    }
}

impl<C, N, B:ToLatticeBoundary> fmt::Debug for Environment<C, N, B>
where
    C: fmt::Debug,
    N: fmt::Debug,
    B: fmt::Debug,
    B::LatticeBoundary: fmt::Debug {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("Environment")
            .field("bounds", &self.bounds)
            .field("cells", &self.cells)
            .field("cell", &self.cell_lattice)
            .field("edge_book", &self.edge_book)
            .field("neighborhood", &self.neighborhood)
            .finish()
    }
}

impl<C, N, B:ToLatticeBoundary> PartialEq for Environment<C, N, B>
where
    C: PartialEq,
    N: PartialEq,
    B: PartialEq,
    B::LatticeBoundary: PartialEq {
    fn eq(&self, other: &Self) -> bool {
        // Assume correctness: no need to check edge_books
        self.cells == other.cells
            && self.cell_lattice == other.cell_lattice
            && self.neighborhood == other.neighborhood
            && self.bounds == other.bounds
    }
}

/// Counts the number of changed cell-cell edges after modifying the environment with [Habitable::grant_position()].
#[derive(Clone, Debug, PartialEq)]
pub struct EdgesUpdate {
    /// Number of cell-cell edges added by granting the position.
    pub added: u16,
    /// Number of cell-cell edges removed by granting the position.
    pub removed: u16
}

/// Helper function used in [Environment::valid_neighbors].
fn valid_neighbors(
    lattice_boundary: &impl Boundary<Coord = isize>,
    neighborhood: &impl Neighborhood,
    pos: Pos<usize>
) -> impl Iterator<Item = Pos<usize>> {
    lattice_boundary.valid_positions(
        neighborhood.neighbors(
            pos.cast_as()
        )
    ).map(|pos| pos.cast_as())
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use crate::base::cell::Cell;
    use crate::cell_container::RelCell;
    use crate::positional::boundaries::UnsafePeriodicBoundary;
    use crate::positional::neighborhood::MooreNeighborhood;
    use crate::positional::pos::Pos;
    use crate::positional::rect::Rect;

    fn make_test_env() -> Environment<Cell, MooreNeighborhood, UnsafePeriodicBoundary<FloatType>> {
        let rect = Rect::new(
            (0., 0.,).into(),
            (10., 10.).into()
        );
        Environment::new_empty(
            MooreNeighborhood::new(1),
            Boundaries::new(UnsafePeriodicBoundary::new(rect.clone()))
        )
    }

    fn add_cell(
        positions: &[Pos<usize>],
        env: &mut Environment<Cell, MooreNeighborhood, UnsafePeriodicBoundary<FloatType>>
    ) -> RelCell<Cell> {
        let mut rel_cell = RelCell{ index: 0, cell: Cell::new_empty(2).into_cell() };
        for &pos in positions {
            rel_cell.cell.shift_position(pos, true, &env.bounds.boundary);
            env.cell_lattice[pos] = Spin::Some(rel_cell.index);
        }
        rel_cell
    }

    #[test]
    fn test_spawn_solid() {
        let mut env = make_test_env();
        let positions = vec![
            Pos::new(1, 1),
            Pos::new(2, 2),
            Pos::new(3, 3),
            Pos::new(1, 1), // duplicate to test deduplication
        ];
        env.spawn_solid(positions.into_iter());
        let solid_count = env
            .cell_lattice
            .iter_values()
            .filter(|&&val| matches!(val, Spin::Solid))
            .count();
        assert_eq!(solid_count, 3); // One was a duplicate
        for pos in &[
            Pos::new(1, 1),
            Pos::new(2, 2),
            Pos::new(3, 3),
        ] {
            assert_eq!(env.cell_lattice[*pos], Spin::Solid);
        }
    }

    #[test]
    fn test_update_edges_adds_and_removes() {
        let mut env = make_test_env();
        let spin = Spin::Some(0);
        env.cell_lattice[Pos::new(5, 5)] = spin;
        let mut edges_update = env.update_edges(Pos::new(5, 5));
        assert_eq!(edges_update.removed, 0);
        assert_eq!(edges_update.added, 8);

        env.cell_lattice[Pos::new(6, 5)] = spin;
        edges_update = env.update_edges(Pos::new(5, 5));
        assert_eq!(edges_update.removed, 1);
        assert_eq!(edges_update.added, 0);

        env.cell_lattice[Pos::new(6, 5)] = Spin::Some(1);
        edges_update = env.update_edges(Pos::new(5, 5));
        assert_eq!(edges_update.removed, 0);
        assert_eq!(edges_update.added, 1);
    }

    #[test]
    fn test_box_cell_positions() {
        let positions = [
            Pos::new(5, 5),
            Pos::new(5, 6),
            Pos::new(6, 5),
            Pos::new(6, 6),
        ];
        let mut env = make_test_env();
        let cell = add_cell(&positions, &mut env);
        let boxed_positions = env.search_cell_box(&cell, 2.);
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
        let mut env = make_test_env();
        let cell = add_cell(&positions, &mut env);
        let contiguous_positions = env.search_cell_contiguous(&cell);

        // Should find all 4 contiguous positions
        assert_eq!(contiguous_positions.len(), positions.len());
        for pos in &positions {
            assert!(contiguous_positions.contains(pos));
        }
    }

    #[test]
    fn test_contiguous_cell_positions_discontiguous() {
        let positions = [
            Pos::new(4, 5),
            Pos::new(5, 4),
            Pos::new(5, 5),
            Pos::new(6, 5),
            Pos::new(5, 6),
            Pos::new(7, 7), // discontiguous point
        ];
        let mut env = make_test_env();
        let cell = add_cell(&positions, &mut env);
        let result = env.search_cell_contiguous(&cell);
        assert_eq!(result.len(), 5);
    }

    #[test]
    fn test_cell_neighbors() {
        let positions = [
            Pos::new(2, 2),
            Pos::new(2, 1),
        ];
        let mut env = make_test_env();
        let rel_cell = add_cell(&positions, &mut env);

        let neighbor_spins = [Spin::Some(rel_cell.index + 1), Spin::Some(rel_cell.index + 2)];
        env.cell_lattice[Pos::new(1, 2)] = neighbor_spins[0];
        env.cell_lattice[Pos::new(2, 0)] = neighbor_spins[1];

        let neighs = env.cell_neighbors(&rel_cell, 2.);

        assert!(neighs.contains(&neighbor_spins[0]));
        assert!(neighs.contains(&neighbor_spins[1]));
        assert!(!neighs.contains(&Spin::Some(rel_cell.index)));
    }
}

