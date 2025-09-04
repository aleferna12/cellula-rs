use crate::cell_container::CellContainer;
use crate::cellular::{Cellular, RelCell};
use crate::constants::Spin;
use crate::lattice_entity::LatticeEntity;
use crate::lattice_entity::LatticeEntity::*;
use crate::positional::boundary::Boundary;
use crate::positional::edge::Edge;
use crate::positional::edge_book::EdgeBook;
use crate::positional::neighbourhood::Neighbourhood;
use crate::positional::pos::Pos;
use crate::spatial::Spatial;
use rustworkx_core::petgraph::prelude::UnGraph;
use std::cmp::max;
use std::collections::{HashMap, HashSet};
use std::f32::consts::PI;
use std::fmt::Debug;

pub struct Environment<C, N, S> {
    pub space: S,
    pub cells: CellContainer<C>,
    pub edge_book: EdgeBook,
    pub neighbourhood: N
}

impl<C, N, S: Spatial> Environment<C, N, S> {
    pub fn new(
        cells: CellContainer<C>,
        space: S,
        neighbourhood: N
    ) -> Self {
        Self {
            space,
            cells,
            neighbourhood,
            edge_book: EdgeBook::new(),
        }
    }
    
    pub fn width(&self) -> usize {
        self.space.cell_lattice().width()
    }

    pub fn height(&self) -> usize {
        self.space.cell_lattice().height()
    }
    
    pub fn spawn_solid(&mut self, positions: impl Iterator<Item = Pos<usize>>) -> usize {
        let mut area = 0;
        for pos in positions {
            if self.space.cell_lattice()[pos] != Medium.discriminant() {
                continue
            }
            self.space.cell_lattice_mut()[pos] = Solid.discriminant();
            area += 1;
        }
        area
    }
    
    pub fn make_border(
        &mut self,
        bottom: bool,
        top: bool,
        left: bool,
        right: bool,
    ) {
        let mut border_positions = Vec::<Pos<usize>>::new();
        if bottom {
            for x in 0..self.width() {
                border_positions.push((x, 0).into());
            }
        }
        if top {
            if self.height() > 1 {
                for x in (0..self.width() - 1).rev() {
                    border_positions.push((x, self.height() - 1).into());
                }
            }
        }
        if left {
            if self.width() > 1 {
                for y in (1..self.height() - 1).rev() {
                    border_positions.push((0, y).into());
                }
            }
        }
        if right {
            for y in 1..self.height() {
                border_positions.push((self.width() - 1, y).into());
            }
        }
        
        self.spawn_solid(border_positions.into_iter());
    }

    pub fn update_edges(&mut self, pos: Pos<usize>) -> EdgesUpdate
    where N: Neighbourhood {
        let mut removed = 0;
        let mut added = 0;
        let spin = self.space.cell_lattice()[pos];
        let valid_neighs = self
            .space
            .lattice_boundary()
            .valid_positions(self.neighbourhood.neighbours(pos.to_isize()))
            .map(|neigh| neigh.to_usize());
        for neigh in valid_neighs {
            let edge = Edge::new(pos, neigh);
            let spin_neigh = self.space.cell_lattice()[neigh];
            // The order of these if statements matter A LOT, dont mess with it
            if spin == spin_neigh {
                if self.edge_book.remove(&edge) {
                    removed += 1;
                }
                continue;
            }
            if spin < LatticeEntity::first_cell_spin() && spin_neigh < LatticeEntity::first_cell_spin() {
                continue;
            }
            if self.edge_book.insert(edge) {
                added += 1;
            }
        }
        EdgesUpdate { added, removed }
    }
}

impl<C: Cellular, N: Neighbourhood, S: Spatial> Environment<C, N, S> {
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
        self.space.cell_lattice().search_box(
            &cell.spin,
            cell.center().to_usize(),
            search_diam,
            self.space.lattice_boundary(),
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
    pub fn search_cell_contiguous(
        &self,
        cell: &RelCell<impl Cellular>,
    ) -> Vec<Pos<usize>> {
        let found = self.space.cell_lattice().search_contiguous(
            &cell.spin,
            cell.center().to_usize(),
            self.space.lattice_boundary(),
            &self.neighbourhood
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

    pub fn search_cell_outline(
        &self,
        cell: &RelCell<impl Cellular>,
        diameter_scaler: f32,
    ) -> Vec<Pos<usize>> {
        let search_diam = (
            diameter_scaler
                * 2.
                * (max(cell.target_area(), cell.area()) as f32 / PI)
                .sqrt()
        ) as usize;
        self.space.cell_lattice().search_outline(
            &cell.spin,
            cell.center().to_usize(),
            search_diam,
            self.space.lattice_boundary(),
            &self.neighbourhood
        )
    }

    pub fn cell_neighbours(
        &self,
        cell: &RelCell<impl Cellular>,
        diameter_scaler: f32,
    ) -> HashSet<Spin> {
        let outline = self.search_cell_outline(
            cell,
            diameter_scaler,
        ).into_iter().map(|pos| { self.space.cell_lattice()[pos] });
        HashSet::from_iter(outline)
    }

    pub fn build_neighbours_graph(&self, diameter_scaler: f32) -> UnGraph<Spin, ()> {
        let mut graph = UnGraph::new_undirected();
        let mut node_map = HashMap::new();

        for cell in self.cells.iter() {
            if !cell.is_alive() {
                continue;
            }

            // Add or retrieve the node for this cell's spin
            let cell_idx = *node_map.entry(cell.spin)
                .or_insert_with(|| graph.add_node(cell.spin));

            let neighs = self.cell_neighbours(
                cell,
                diameter_scaler,
            );

            for neigh_spin in neighs {
                if neigh_spin < LatticeEntity::first_cell_spin() {
                    continue;
                }
                // Add or retrieve the node for the neighbor spin
                let neigh_idx = *node_map.entry(neigh_spin)
                    .or_insert_with(|| graph.add_node(neigh_spin));

                graph.update_edge(cell_idx, neigh_idx, ());
            }
        }
        graph
    }
}

pub struct EdgesUpdate {
    pub added: u16,
    pub removed: u16
}

#[derive(Debug)]
pub enum DivisionError {
    NewCellTooSmall,
    NewCellTooBig
}

// #[cfg(test)]
// pub mod tests {
//     use super::*;
//     use crate::cell::Cell;
//     use crate::positional::pos::Pos;
//     use crate::positional::rect::Rect;
// 
//     fn make_env_for_division() -> Environment<Cell<MockGenome>, NeighbourhoodType, BoundaryType> {
//         let env = Environment::new(
//             2.0,
//             10,
//             CellContainer::new(
//                 4,
//                 true,
//                 false
//             ),
//             Space::new(BoundaryType::new(Rect::new(
//                 (0., 0.,).into(),
//                 (100., 100.).into()
//             ))).expect("failed to make test `Space`"),
//             NeighbourhoodType::new(1)
//         );
//         env
//     }
// 
//     #[test]
//     fn test_spawn_solid() {
//         let mut env = Environment::new_empty_test(10, 10);
//         let positions = vec![
//             Pos::new(1, 1),
//             Pos::new(2, 2),
//             Pos::new(3, 3),
//             Pos::new(1, 1), // duplicate to test deduplication
//         ];
//         let area = env.spawn_solid(positions.into_iter());
//         assert_eq!(area, 3); // One was a duplicate
//         for pos in &[
//             Pos::new(1, 1),
//             Pos::new(2, 2),
//             Pos::new(3, 3),
//         ] {
//             assert_eq!(env.space.cell_lattice[*pos], Solid.discriminant());
//         }
//     }
// 
//     #[test]
//     fn test_update_edges_adds_and_removes() {
//         let mut env = Environment::new_empty_test(10, 10);
//         let spin = LatticeEntity::first_cell_spin();
//         env.space.cell_lattice[Pos::new(5, 5)] = spin;
//         let mut edges_update = env.update_edges(Pos::new(5, 5));
//         assert_eq!(edges_update.removed, 0);
//         assert_eq!(edges_update.added, 8);
//         
//         env.space.cell_lattice[Pos::new(6, 5)] = spin;
//         edges_update = env.update_edges(Pos::new(5, 5));
//         assert_eq!(edges_update.removed, 1);
//         assert_eq!(edges_update.added, 0);
// 
//         env.space.cell_lattice[Pos::new(6, 5)] = spin + 1;
//         edges_update = env.update_edges(Pos::new(5, 5));
//         assert_eq!(edges_update.removed, 0);
//         assert_eq!(edges_update.added, 1);
//     }
// 
//     #[test]
//     fn test_divide_cell() {
//         let mut env = make_env_for_division();
// 
//         let rect = Rect::new(Pos::new(20, 20), Pos::new(23, 23));
//         env.spawn_rect_cell(
//             rect, 
//             Cell::new_empty(4, 8, MockGenome::new(0))
//         );
// 
//         let spin = LatticeEntity::first_cell_spin();
//         let result = env.divide_cell(spin);
//         assert!(result.is_ok());
//         let new_cell = result.unwrap();
//         assert_ne!(new_cell.spin, spin);
//     }
// 
//     #[test]
//     fn test_reproduce() {
//         let mut env = make_env_for_division();
//         env.spawn_rect_cell(
//             Rect::new(Pos::new(30, 30), Pos::new(33, 33)),
//             Cell::new_empty(4, 8, MockGenome::new(0))
//         );
// 
//         let divided_spins = env.reproduce();
//         assert_eq!(divided_spins.len(), 1);
//     }
// 
//     #[test]
//     fn test_reproduce_limit_population() {
//         let mut env = make_env_for_division();
//         env.max_cells = 1;
//         env.spawn_rect_cell(
//             Rect::new(Pos::new(30, 30), Pos::new(33, 33)),
//             Cell::new_empty(4, 8, MockGenome::new(0))
//         );
// 
//         let divided_spins = env.reproduce();
//         assert_eq!(divided_spins.len(), 0);
//     }
// }
// 
