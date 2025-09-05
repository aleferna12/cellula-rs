use crate::cell_container::CellContainer;
use crate::cellular::{Cellular, RelCell};
use crate::constants::Spin;
use crate::lattice_entity::LatticeEntity;
use crate::lattice_entity::LatticeEntity::*;
use crate::positional::boundary::{ToLatticeBoundary, Boundary};
use crate::positional::edge::Edge;
use crate::positional::edge_book::EdgeBook;
use crate::positional::neighbourhood::Neighbourhood;
use crate::positional::pos::Pos;
use crate::space::{Spatial};
use rustworkx_core::petgraph::prelude::UnGraph;
use std::cmp::max;
use std::collections::{HashMap, HashSet};
use std::f32::consts::PI;
use rand::Rng;
use rand_xoshiro::Xoshiro256StarStar;
use crate::lattice::Lattice;
use crate::positional::rect::Rect;

pub trait Habitable {
    type Cell;
    fn cell_lattice(&self) -> &Lattice<Spin>;
    fn cell_lattice_mut(&mut self) -> &mut Lattice<Spin>;
    fn update_edges(&mut self, pos: Pos<usize>) -> EdgesUpdate;
}

pub struct Environment<C, N, S> {
    pub space: S,
    pub cells: CellContainer<C>,
    pub cell_lattice: Lattice<Spin>,
    pub edge_book: EdgeBook,
    pub neighbourhood: N,
    pub rng: Xoshiro256StarStar,
    pub cell_search_scaler: f32
}

impl<C, N, S: Spatial> Environment<C, N, S> {
    pub fn new(
        cell_search_scaler: f32,
        cells: CellContainer<C>,
        neighbourhood: N,
        space: S,
        cell_lattice: Lattice<Spin>,
        rng: Xoshiro256StarStar,
    ) -> Self {
        Self {
            space,
            cells,
            neighbourhood,
            cell_lattice,
            cell_search_scaler,
            rng,
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
    pub fn grant_position(
        &mut self,
        pos: Pos<usize>,
        to: &mut RelCell<C>,
    ) -> EdgesUpdate {
        to.shift_position(pos, true, self.space.boundary());
        if let SomeCell(from) = self.cells.get_entity_mut(self.cell_lattice()[pos]) {
            from.shift_position(pos, false, self.space.boundary());
        }
        // Executes the copy
        self.cell_lattice_mut()[pos] = to.spin;
        self.update_edges(pos)
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
        cell: &RelCell<impl Cellular>
    ) -> impl Iterator<Item = Pos<usize>> {
        let search_diam = (
            self.cell_search_scaler
                * 2.
                * (max(cell.target_area(), cell.area()) as f32 / PI)
                .sqrt()
        ) as usize;
        self.cell_lattice.search_box(
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
    pub fn search_cell_box(&self, cell: &RelCell<impl Cellular>) -> Vec<Pos<usize>> {
        let found: Vec<_> = self.search_cell_box_iter(cell).collect();
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
        cell: &RelCell<impl Cellular>
    ) -> Vec<Pos<usize>> {
        let search_diam = (
            self.cell_search_scaler
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
        cell: &RelCell<impl Cellular>
    ) -> HashSet<Spin> {
        let outline = self
            .search_cell_outline(cell)
            .into_iter()
            .map(|pos| self.space.cell_lattice()[pos]);
        HashSet::from_iter(outline)
    }

    pub fn build_neighbours_graph(&self) -> UnGraph<Spin, ()> {
        let mut graph = UnGraph::new_undirected();
        let mut node_map = HashMap::new();

        for cell in self.cells.iter() {
            if !cell.is_alive() {
                continue;
            }

            // Add or retrieve the node for this cell's spin
            let cell_idx = *node_map.entry(cell.spin)
                .or_insert_with(|| graph.add_node(cell.spin));

            let neighs = self.cell_neighbours(cell);

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
    }pub fn spawn_cell_random(
        &mut self,
        cell_area: u32,
        empty_cell: C,
        rng: &mut impl Rng,
    ) -> Option<&RelCell<C>> {
        let cell_side = (cell_area as f32).sqrt() as usize;
        let pos = self.cell_lattice().random_pos(rng);
        self.spawn_rect_cell(
            Rect::new(
                pos,
                (pos.x + cell_side, pos.y + cell_side).into()
            ),
            empty_cell.clone()
        )
    }

    pub fn spawn_rect_cell(&mut self, rect: Rect<usize>, mut empty_cell: C) -> Option<&RelCell<C>> {
        if !self.can_add_cell() {
            return None;
        }

        let spin = self.cells.next_spin();
        for pos in rect.iter_positions() {
            if let Some(valid_pos) = self.space.lattice_boundary().valid_pos(pos.to_isize()) {
                let lat_pos = valid_pos.to_usize();
                if self.cell_lattice_mut()[lat_pos] != Medium.discriminant() {
                    continue
                }
                self.cell_lattice_mut()[lat_pos] = spin;
                self.update_edges(lat_pos);
                empty_cell.shift_position(
                    lat_pos,
                    true,
                    self.space.boundary()
                );
            }
        }
        if empty_cell.area() == 0 {
            return None;
        }
        Some(self.cells.push(empty_cell, None))
    }

    pub fn kill_cell(&mut self, cell: &mut RelCell<C>) {
        for pos in self.search_cell_box(cell) {
            // TODO!: Parameterize chance of medium
            if self.rng.random::<f32>() < 0.1 {
                self.cell_lattice_mut()[pos] = Medium.discriminant();
            }
        }
        cell.apoptosis();
    }

    pub fn wipe_out(&mut self) {
        self.cells.wipe_out();
        self.cell_lattice_mut().iter_values_mut().for_each(|value| {
            if *value >= LatticeEntity::first_cell_spin() {
                *value = Medium.discriminant();
            }
        });
    }

    pub fn divide_cell(&mut self, mother: &RelCell<C>, empty_cell: C) -> Result<&RelCell<C>, DivisionError> {
        let mut mom_clone = mother.clone();

        let new_positions: Vec<_> = self
            .search_cell_box(&mom_clone)
            .into_iter()
            .filter(|pos| {
                // TODO!: use principal component to determine division axis
                //  current algorithm hands out all x positions to the right of the cell centre to the new cell
                self.space.boundary().displacement(Pos::new(pos.x as f32, pos.y as f32), mom_clone.center()).0 > 0.
            })
            .collect();
        if new_positions.is_empty() {
            return  Err(DivisionError::NewCellTooSmall);
        }
        if new_positions.len() >= mother.area() as usize {
            return Err(DivisionError::NewCellTooBig);
        }

        let mut new_cell = self.cells.push(empty_cell, Some(mother.spin));
        for pos in new_positions {
            if mom_clone.area() == 1 {
                return Err(DivisionError::NewCellTooBig);
            }
            self.grant_position(
                pos,
                &mut new_cell,
            );
        }
        Ok(new_cell)
    }
}

impl<C, N, S: Spatial> Spatial for Environment<C, N, S> {
    type Boundary = S::Boundary;

    fn boundary(&self) -> &Self::Boundary {
        self.space.boundary()
    }

    fn lattice_boundary(&self) -> &<Self::Boundary as ToLatticeBoundary>::LatticeBoundary {
        self.space.lattice_boundary()
    }
}

impl <C, N, S: Spatial> Habitable for Environment<C, N, S> {
    type Cell = C;

    fn cell_lattice(&self) -> &Lattice<Spin> {
        &self.cell_lattice
    }

    fn cell_lattice_mut(&mut self) -> &mut Lattice<Spin> {
        &mut self.cell_lattice
    }

    fn update_edges(&mut self, pos: Pos<usize>) -> EdgesUpdate {
        self.update_edges(pos)
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

#[cfg(test)]
pub mod tests {
    use crate::cellular::BasicCell;
    use crate::positional::boundary::UnsafePeriodicBoundary;
    use crate::positional::neighbourhood::MooreNeighbourhood;
    use super::*;
    use crate::positional::pos::Pos;
    use crate::positional::rect::Rect;
    use crate::space::Space;

    fn make_test_env() -> Environment<BasicCell, MooreNeighbourhood, Space<UnsafePeriodicBoundary<f32>>> {
        let env = Environment::new(
            CellContainer::default(),
            MooreNeighbourhood::new(1),
            Space::new(UnsafePeriodicBoundary::new(Rect::new(
                (0., 0.,).into(),
                (10., 10.).into()
            ))).expect("failed to make test `Space`")
        );
        env
    }

    fn add_cell(
        positions: &[Pos<usize>],
        env: &mut Environment<BasicCell, MooreNeighbourhood, Space<UnsafePeriodicBoundary<f32>>>
    ) -> RelCell<BasicCell> {
        let mut cell = RelCell::mock(BasicCell::new_empty(2));
        for &pos in positions {
            cell.shift_position(pos, true, env.space.boundary());
            env.space.cell_lattice_mut()[pos] = cell.spin;
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
        let area = env.spawn_solid(positions.into_iter());
        assert_eq!(area, 3); // One was a duplicate
        for pos in &[
            Pos::new(1, 1),
            Pos::new(2, 2),
            Pos::new(3, 3),
        ] {
            assert_eq!(env.space.cell_lattice()[*pos], Solid.discriminant());
        }
    }

    #[test]
    fn test_update_edges_adds_and_removes() {
        let mut env = make_test_env();
        let spin = LatticeEntity::first_cell_spin();
        env.space.cell_lattice_mut()[Pos::new(5, 5)] = spin;
        let mut edges_update = env.update_edges(Pos::new(5, 5));
        assert_eq!(edges_update.removed, 0);
        assert_eq!(edges_update.added, 8);

        env.space.cell_lattice_mut()[Pos::new(6, 5)] = spin;
        edges_update = env.update_edges(Pos::new(5, 5));
        assert_eq!(edges_update.removed, 1);
        assert_eq!(edges_update.added, 0);

        env.space.cell_lattice_mut()[Pos::new(6, 5)] = spin + 1;
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
        let boxed_positions = env.search_cell_box(&cell, 2.0);
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

        let neighbour_spins = [cell.spin + 1, cell.spin + 2];
        env.space.cell_lattice_mut()[Pos::new(1, 2)] = neighbour_spins[0];
        env.space.cell_lattice_mut()[Pos::new(2, 0)] = neighbour_spins[1];

        let neighs = env.cell_neighbours(&cell, 1.0);

        assert!(neighs.contains(&neighbour_spins[0]));
        assert!(neighs.contains(&neighbour_spins[1]));
        assert!(!neighs.contains(&cell.spin));
    }
}

