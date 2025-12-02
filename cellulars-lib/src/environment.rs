//! Contains logic associated with [Environment].

use crate::basic_cell::{Alive, BasicCell, Cellular, RelCell};
use crate::cell_container::CellContainer;
use crate::constants::CellIndex;
use crate::habitable::Habitable;
use crate::lattice::Lattice;
use crate::positional::boundaries::{Boundaries, Boundary, ToLatticeBoundary};
use crate::positional::edge::Edge;
use crate::positional::edge_book::EdgeBook;
use crate::positional::neighbourhood::Neighbourhood;
use crate::positional::pos::Pos;
use crate::positional::rect::RectConversionError;
use crate::spin::Spin;
use rustworkx_core::petgraph::prelude::UnGraph;
use std::cmp::max;
use std::collections::{HashMap, HashSet};
use std::f32::consts::PI;

/// An environment where cells are spatially and relationally localised.
pub struct Environment<C, N, B: ToLatticeBoundary> {
    /// Boundaries of the environment.
    ///
    /// These are used to validate that a position can be used to access information in the environment.
    pub bounds: Boundaries<B>,
    /// Cell container with spins matching those in [Environment::cell_lattice].
    pub cells: CellContainer<C>,
    /// Cell lattice with spins matching those in [Environment::cells].
    pub cell_lattice: Lattice<Spin>,
    /// Edge book containing all positions at cell-cell interfaces.
    pub edge_book: EdgeBook,
    /// Neighbourhood of the environment.
    pub neighbourhood: N
}

impl<C, N, B: ToLatticeBoundary<Coord = f32>> Environment<C, N, B> {
    /// Makes a new empty environment with no cells.
    pub fn new_empty(
        neighbourhood: N,
        bounds: Boundaries<B>,
    ) -> Result<Self, RectConversionError> {
        Ok(Self {
            cell_lattice: Lattice::new(bounds.boundary.rect().clone().try_into()?),
            edge_book: EdgeBook::new(),
            cells: CellContainer::new(),
            bounds,
            neighbourhood
        })
    }
}

impl<C, N, B: ToLatticeBoundary> Environment<C, N, B> {
    /// Makes a new environment from its components.
    pub fn new(
        cells: CellContainer<C>,
        cell_lattice: Lattice<Spin>,
        neighbourhood: N,
        bounds: Boundaries<B>,
    ) -> Self {
        Self {
            cell_lattice,
            bounds,
            cells,
            neighbourhood,
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

    /// Returns an iterator over all sites in the neighbourhood of `pos` that are within lattice boundaries.
    pub fn valid_neighbours(&self, pos: Pos<usize>) -> impl Iterator<Item = Pos<usize>>
    where N: Neighbourhood {
        self.bounds.lattice_boundary.valid_positions(
            self.neighbourhood.neighbours(pos.to_isize())
        ).map(|pos| pos.to_usize())
    }
}

impl<C: Cellular, N: Neighbourhood, B: ToLatticeBoundary> Environment<C, N, B> {
    // This function returns a Vec so that we can check that the site number matches
    /// Searches for all cell positions by creating a box around the cell and iterating all the positions inside it.
    ///
    /// May fail if `radius_scaler` is too small, in which case logs a warning.
    pub fn search_cell_box(&self, cell: &RelCell<impl Cellular>, search_scaler: f32) -> Vec<Pos<usize>> {
        let search_diam = (
            search_scaler
                * 2.
                * (max(cell.target_area(), cell.area()) as f32 / PI)
                .sqrt()
        ) as usize;

        let found: Vec<_> = self.cell_lattice.search_box(
            &Spin::Some(cell.index),
            cell.center().to_usize(),
            search_diam,
            &self.bounds.lattice_boundary,
        ).collect();

        if found.len() != cell.area() as usize {
            log::warn!(
                "Only found {} positions out of the {} expected for cell with index {} \
                (try to increase `search-radius`)",
                found.len(),
                cell.area(),
                cell.index
            )
        }
        found
    }

    /// Searches for all cell positions with a BFS algorithm to traverse the lattice sites.
    ///
    /// Is considerably slower than [Environment::search_cell_box()] and may fail if the cell is not contiguous
    /// or if the cell centre is not a cell position.
    pub fn search_cell_contiguous(
        &self,
        cell: &RelCell<impl Cellular>,
    ) -> Vec<Pos<usize>> {
        let found = self.cell_lattice.search_contiguous(
            &Spin::Some(cell.index),
            cell.center().to_usize(),
            &self.bounds.lattice_boundary,
            &self.neighbourhood
        );

        if found.len() != cell.area() as usize {
            log::warn!(
                "Only found {} positions out of the {} expected for cell with index {} \
                (cell might be discontiguous)",
                found.len(),
                cell.area(),
                cell.index
            )
        }
        found
    }

    /// Searches for all cell positions that interface a different spin.
    pub fn search_cell_outline(
        &self,
        cell: &RelCell<impl Cellular>,
        search_scaler: f32
    ) -> Vec<Pos<usize>> {
        let search_diam = (
            search_scaler
                * 2.
                * (max(cell.target_area(), cell.area()) as f32 / PI)
                .sqrt()
        ) as usize;

        self.cell_lattice.search_outline(
            &Spin::Some(cell.index),
            cell.center().to_usize(),
            search_diam,
            &self.bounds.lattice_boundary,
            &self.neighbourhood
        )
    }

    /// Finds all cells adjacent to `cell` using a search algorithm.
    pub fn cell_neighbours(
        &self,
        cell: &RelCell<impl Cellular>,
        search_scaler: f32
    ) -> HashSet<Spin> {
        let outline = self
            .search_cell_outline(cell, search_scaler)
            .into_iter()
            .map(|pos| self.cell_lattice[pos]);
        HashSet::from_iter(outline)
    }

    /// Builds an undirected graph of cell-cell neighbours (cells that are adjacent).
    pub fn build_neighbours_graph(&self, search_scaler: f32) -> UnGraph<CellIndex, ()> {
        let mut graph = UnGraph::new_undirected();
        let mut node_map = HashMap::new();

        for cell in self.cells.iter() {
            if !cell.is_valid() {
                continue;
            }

            // Add or retrieve the node for this cell's index
            let cell_idx = *node_map.entry(cell.index)
                .or_insert_with(|| graph.add_node(cell.index));

            let neighs = self.cell_neighbours(cell, search_scaler);

            for neigh_spin in neighs {
                let neigh_index = match neigh_spin {
                    Spin::Some(cell_index) => cell_index,
                    _ => continue,
                };
                // Add or retrieve the node for the neighbor index
                let neigh_idx = *node_map.entry(neigh_index)
                    .or_insert_with(|| graph.add_node(neigh_index));

                graph.update_edge(cell_idx, neigh_idx, ());
            }
        }
        graph
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
        let valid_neighs = self
            .bounds
            .lattice_boundary
            .valid_positions(self.neighbourhood.neighbours(pos.to_isize()))
            .map(|neigh| neigh.to_usize());
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

impl<N: Neighbourhood, B: ToLatticeBoundary<Coord = f32>> Habitable for Environment<BasicCell, N, B> {
    type Cell = BasicCell;

    fn env(&self) -> &Environment<Self::Cell, impl Neighbourhood, impl ToLatticeBoundary> {
        self
    }

    fn env_mut(&mut self) -> &mut Environment<Self::Cell, impl Neighbourhood, impl ToLatticeBoundary> {
        self
    }

    fn grant_position(
        &mut self,
        pos: Pos<usize>,
        to: Spin
    ) -> EdgesUpdate {
        if let Spin::Some(index) = to {
            self.cells.get_cell_mut(index).shift_position(pos, true, &self.bounds.boundary);
        }
        if let Spin::Some(index) = self.cell_lattice[pos] {
            let from_cell = self.cells.get_cell_mut(index);
            from_cell.shift_position(pos, false, &self.bounds.boundary);
            if from_cell.area() == 0 {
                from_cell.apoptosis();
            }
        }
        // Executes the copy
        self.cell_lattice[pos] = to;
        self.update_edges(pos)
    }
}

impl<C, N, B: ToLatticeBoundary> Clone for Environment<C, N, B>
where
    B: Clone,
    B::LatticeBoundary: Clone,
    CellContainer<C>: Clone,
    N: Clone,
{
    fn clone(&self) -> Self {
        Environment {
            bounds: self.bounds.clone(),
            cells: self.cells.clone(),
            cell_lattice: self.cell_lattice.clone(),
            edge_book: self.edge_book.clone(),
            neighbourhood: self.neighbourhood.clone(),
        }
    }
}

/// Counts the number of changed cell-cell edges after modifying the environment with [Habitable::grant_position()].
pub struct EdgesUpdate {
    /// Number of cell-cell edges added by granting the position.
    pub added: u16,
    /// Number of cell-cell edges removed by granting the position.
    pub removed: u16
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use crate::basic_cell::BasicCell;
    use crate::positional::boundaries::UnsafePeriodicBoundary;
    use crate::positional::neighbourhood::MooreNeighbourhood;
    use crate::positional::pos::Pos;
    use crate::positional::rect::Rect;

    fn make_test_env() -> Environment<BasicCell, MooreNeighbourhood, UnsafePeriodicBoundary<f32>> {
        let rect = Rect::new(
            (0., 0.,).into(),
            (10., 10.).into()
        );
        Environment::new_empty(
            MooreNeighbourhood::new(1),
            Boundaries::new(UnsafePeriodicBoundary::new(rect.clone())).expect("failed to make boundaries"),
        ).expect("failed to make test environment")
    }

    fn add_cell(
        positions: &[Pos<usize>],
        env: &mut Environment<BasicCell, MooreNeighbourhood, UnsafePeriodicBoundary<f32>>
    ) -> RelCell<BasicCell> {
        let mut cell = RelCell::mock(BasicCell::new_empty(2));
        for &pos in positions {
            cell.shift_position(pos, true, &env.bounds.boundary);
            env.cell_lattice[pos] = Spin::Some(cell.index);
        }
        cell
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
            Pos::new(5, 5),
            Pos::new(5, 6),
            Pos::new(7, 7), // discontiguous point
        ];
        let mut env = make_test_env();
        let cell = add_cell(&positions, &mut env);
        let result = env.search_cell_contiguous(&cell);
        assert_eq!(result.len(), 2);
    }

    #[test]
    fn test_cell_neighbours() {
        let positions = [
            Pos::new(2, 2),
            Pos::new(2, 1),
        ];
        let mut env = make_test_env();
        let cell = add_cell(&positions, &mut env);

        let neighbour_spins = [Spin::Some(cell.index + 1), Spin::Some(cell.index + 2)];
        env.cell_lattice[Pos::new(1, 2)] = neighbour_spins[0];
        env.cell_lattice[Pos::new(2, 0)] = neighbour_spins[1];

        let neighs = env.cell_neighbours(&cell, 2.);

        assert!(neighs.contains(&neighbour_spins[0]));
        assert!(neighs.contains(&neighbour_spins[1]));
        assert!(!neighs.contains(&Spin::Some(cell.index)));
    }
}

