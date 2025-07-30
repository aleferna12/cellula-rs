// TODO!: Revise all trait bounds of methods of Environment 
use crate::cell::{CanDivide, Cell, CellLike, ChemSniffer, RelCell};
use crate::cell_container::CellContainer;
use crate::constants::{BoundaryType, NeighbourhoodType, Spin};
use crate::environment::LatticeEntity::*;
use crate::genome::MockGenome;
use crate::positional::boundary::{AsLatticeBoundary, Boundary};
use crate::positional::edge::Edge;
use crate::positional::edge_book::EdgeBook;
use crate::positional::neighbourhood::Neighbourhood;
use crate::positional::pos::Pos;
use crate::positional::rect::Rect;
use crate::space::Space;
use std::fmt::Debug;

pub struct Environment<C, N, B: AsLatticeBoundary> {
    pub space: Space<B>,
    pub cells: CellContainer<C>,
    pub edge_book: EdgeBook,
    pub neighbourhood: N,
    pub update_period: u32,
    pub cell_search_radius: f32,
    pub population_exploded: bool,
    pub max_cells: Spin
}

impl<C, N, B: AsLatticeBoundary> Environment<C, N, B> {
    pub fn new(
        update_period: u32,
        cell_search_radius: f32,
        max_cells: Spin,
        enclose: bool,
        cells: CellContainer<C>,
        space: Space<B>,
        neighbourhood: N
    ) -> Self {
        let mut env = Self {
            space,
            cells,
            neighbourhood,
            edge_book: EdgeBook::new(),
            max_cells,
            update_period,
            cell_search_radius,
            population_exploded: false
        };

        if enclose {
            env.make_border();
        }
    
        env.make_chem_gradient();
        env
    }

    pub fn time_to_update(&self, time_step: u32) -> bool {
        time_step % self.update_period == 0
    }
    
    pub fn width(&self) -> usize {
        self.space.cell_lattice.width()
    }

    pub fn height(&self) -> usize {
        self.space.cell_lattice.height()
    }
    
    pub fn spawn_solid(&mut self, positions: impl Iterator<Item = Pos<usize>>) -> usize {
        let mut area = 0;
        for pos in positions {
            if self.space.cell_lattice[pos] != Medium.discriminant() {
                continue
            }
            self.space.cell_lattice[pos] = Solid.discriminant();
            area += 1;
        }
        area
    }
    
    pub fn make_border(&mut self) {
        let mut border_positions = Vec::<Pos<usize>>::new();
        for x in 0..self.width() {
            border_positions.push((x, 0).into());
        }
        for y in 1..self.height() {
            border_positions.push((self.width() - 1, y).into());
        }
        if self.width() > 1 {
            for y in (1..self.height() - 1).rev() {
                border_positions.push((0, y).into());
            }
        }
        if self.height() > 1 {
            for x in (0..self.width() - 1).rev() {
                border_positions.push((x, self.height() - 1).into());
            }
        }
        
        self.spawn_solid(border_positions.into_iter());
    }
    
    pub fn make_chem_gradient(&mut self) {
        for row in 0..self.height() {
            for col in 0..self.width() {
                self.space.chem_lattice[(col, row).into()] = row.try_into().expect("lattice is too big");
            }
        }
    }
}

impl<C, N, B> Environment<C, N, B> 
where 
    C: CellLike
        + CanDivide
        + ChemSniffer
        + Clone,
    B: AsLatticeBoundary<Coord = f32> {
    // With some unsafe code we can return Vec<&RelCell> from this function, but it would
    // require that self.divide_cell never invalidates any references to self.cells
    // we need thorough testing of self.divide_cells to make this change, and the performance
    // gain is minimal (although the ergonomic gains are significant)
    pub fn reproduce(&mut self) -> Vec<Spin> {
        let mut divide = vec![];
        for cell in &self.cells {
            // Currently cells don't need to express the dividing type to divide, they just need to be big enough
            if cell.area() >= cell.divide_area() {
                divide.push(cell.spin);
            }
        }
        divide.into_iter().filter_map(|spin| {
            if self.cells.n_cells() >= self.max_cells {
                if !self.population_exploded {
                    log::warn!(
                        "Population exceeded maximum threshold `max-cells={}` during cell division", 
                        {self.max_cells}
                    );
                    log::warn!("This warning will be suppressed from now on");
                    self.population_exploded = true;
                }
                return None;
            }
            match self.divide_cell(spin) {
                Err(e) => {
                    log::warn!("Failed to divide cell with spin {spin} with error `{e:?}`");
                    None
                },
                Ok(cell) => Some(cell.spin)
            }
        }).collect()
    }

    // We take spin here because this operation is not safe with &Cell (pushing to vec can cause reallocation)
    pub fn divide_cell(&mut self, spin: Spin) -> Result<&RelCell<C>, DivisionError> {
        let new_spin = self.cells.next_spin();
        let cell_target_area = self.cells.target_area;
        let mut mom_clone = self
            .cells
            .get_entity_mut(spin)
            .expect_cell(&format!("passed non-cell with spin {spin} to `divide_cel()`"))
            .clone();

        let mut new_cell = mom_clone.birth();
        let new_positions: Vec<_> = self
            .space
            .box_cell_positions(&mom_clone, self.cell_search_radius)
            .into_iter()
            .filter(|pos| {
                // TODO!: use principal component to determine division axis
                //  current algorithm hands out all x positions to the right of the cell centre to the new cell
                self.space.bound.displacement(Pos::new(pos.x as f32, pos.y as f32), mom_clone.center()).0 > 0.
            })
            .collect();
        for pos in new_positions {
            if mom_clone.area() == 1 {
                return Err(DivisionError::NewCellTooBig);
            }
            // TODO: this is basically the same as executing a lattice copy, unify the APIs
            self.space.cell_lattice[pos] = new_spin;
            let chem_at = self.space.chem_lattice[pos] as f32;
            new_cell.shift_position(
                pos,
                true,
                &self.space.bound
            );
            new_cell.shift_chem(
                pos,
                chem_at,
                true,
                &self.space.bound
            );
            mom_clone.shift_position(
                pos,
                false,
                &self.space.bound
            );
        }
        if new_cell.area() > 0 {
            let mom_spin = mom_clone.spin;
            mom_clone.set_target_area(cell_target_area);
            self.cells.replace(mom_clone);
            Ok(self.cells.push(new_cell, Some(mom_spin)))
        } else {
            Err(DivisionError::NewCellTooSmall)
        }
    }
}

impl<C, N: Neighbourhood, B: AsLatticeBoundary<Coord = f32>> Environment<C, N, B> {
    pub fn update_edges(&mut self, pos: Pos<usize>) -> (u16, u16) {
        let mut removed = 0;
        let mut added = 0;
        let spin = self.space.cell_lattice[pos];
        let valid_neighs = self
            .space
            .lat_bound
            .valid_positions(self.neighbourhood.neighbours(pos.into()))
            .map(|neigh| neigh.into());
        for neigh in valid_neighs {
            let edge = Edge::new(pos, neigh);
            let spin_neigh = self.space.cell_lattice[neigh];
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
        (removed, added)
    }

    pub fn spawn_rect_cell(&mut self, rect: Rect<usize>, mut empty_cell: C) -> Option<&RelCell<C>>
    where C: CellLike {
        let spin = self.cells.n_cells() as Spin + LatticeEntity::first_cell_spin();

        for pos in rect.iter_positions() {
            if let Some(valid_pos) = self.space.lat_bound.valid_pos(pos.into()) {
                let lat_pos = valid_pos.into();
                if self.space.cell_lattice[lat_pos] != Medium.discriminant() {
                    continue
                }
                self.space.cell_lattice[lat_pos] = spin;
                self.update_edges(lat_pos);
                empty_cell.shift_position(
                    lat_pos,
                    true,
                    &self.space.bound
                );
            }
        }
        if empty_cell.area() == 0 {
            return None;
        }
        Some(self.cells.push(empty_cell, None))
    }
}

impl Environment<Cell<MockGenome>, NeighbourhoodType, BoundaryType> {
    /// Empty environment with arbitrary cell parameters for testing and benchmarking.
    ///
    /// Do not use this in production, no cells can be added to an environment created through this method.
    pub fn new_empty_test(width:usize, height: usize) -> Self {
        Environment::new(
            0,
            0.,
            0,
            false,
            CellContainer::new(
                0,
                false,
                false
            ),
            Space::new(BoundaryType::new(Rect::new(
                (0., 0.,).into(),
                (width as f32, height as f32).into()
            ))).expect("failed to make test `Space`"),
            NeighbourhoodType::new(1)
        )
    }
}

#[derive(Debug)]
pub enum DivisionError {
    NewCellTooSmall,
    NewCellTooBig
}

/// This enum represents anything that can be on the cell lattice.
#[derive(Debug, Copy, Clone)]
pub enum LatticeEntity<C> {
    Solid,
    Medium,
    SomeCell(C),
}

impl<C> LatticeEntity<C> {
    pub fn map<D, F: FnOnce(C) -> D>(self, f: F) -> LatticeEntity<D> {
        match self {
            SomeCell(c) => SomeCell(f(c)),
            Medium => Medium,
            Solid => Solid,
        }
    }
}

impl<C> LatticeEntity<&RelCell<C>> {
    /// Maps the `LatticeEntity` to a unique spin value.
    pub fn spin(&self) -> Spin {
        match self {
            SomeCell(cell) => cell.spin,
            Medium => Medium.discriminant(),
            Solid => Solid.discriminant()
        }
    }
}

impl<C> LatticeEntity<C> {
    pub fn unwrap_cell(self) -> C
    where C: Debug {
        match self {
            SomeCell(cell) => cell,
            _ => panic!("called `LatticeEntity::unwrap_cell()` on a `{self:?}` value")
        }
    }

    pub fn expect_cell(self, message: &str) -> C {
        match self {
            SomeCell(cell) => cell,
            _ => panic!("{}", message)
        }
    }
}

impl LatticeEntity<()> {
    /// Returns the first spin that corresponds to a cell.
    /// 
    /// This is required to be larger than `Medium::spin()` and `Solid::spin()`.
    pub fn first_cell_spin() -> Spin {
        2
    }
    
    pub fn discriminant(&self) -> Spin {
        match self { 
            SomeCell(_) => 2,
            Medium => 0,
            Solid => 1
        }
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;
    use crate::cell::Cell;
    use crate::genome::MockGenome;
    use crate::positional::pos::Pos;
    use crate::positional::rect::Rect;

    fn make_env_for_division() -> Environment<Cell<MockGenome>, NeighbourhoodType, BoundaryType> {
        let env = Environment::new(
            1,
            2.0,
            10,
            false,
            CellContainer::new(
                4,
                true,
                false
            ),
            Space::new(BoundaryType::new(Rect::new(
                (0., 0.,).into(),
                (100., 100.).into()
            ))).expect("failed to make test `Space`"),
            NeighbourhoodType::new(1)
        );
        env
    }

    #[test]
    fn test_spawn_solid() {
        let mut env = Environment::new_empty_test(10, 10);
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
            assert_eq!(env.space.cell_lattice[*pos], Solid.discriminant());
        }
    }

    #[test]
    fn test_update_edges_adds_and_removes() {
        let mut env = Environment::new_empty_test(10, 10);
        let spin = LatticeEntity::first_cell_spin();
        env.space.cell_lattice[Pos::new(5, 5)] = spin;
        let (removed, added) = env.update_edges(Pos::new(5, 5));
        assert_eq!(removed, 0);
        assert_eq!(added, 8);
        
        env.space.cell_lattice[Pos::new(6, 5)] = spin;
        let (removed, added) = env.update_edges(Pos::new(5, 5));
        assert_eq!(removed, 1);
        assert_eq!(added, 0);

        env.space.cell_lattice[Pos::new(6, 5)] = spin + 1;
        let (removed, added) = env.update_edges(Pos::new(5, 5));
        assert_eq!(removed, 0);
        assert_eq!(added, 1);
    }

    #[test]
    fn test_divide_cell() {
        let mut env = make_env_for_division();

        let rect = Rect::new(Pos::new(20, 20), Pos::new(22, 22));
        env.spawn_rect_cell(
            rect, 
            Cell::new_empty(4, 8, MockGenome::new(0))
        );

        let spin = LatticeEntity::first_cell_spin();
        let result = env.divide_cell(spin);
        assert!(result.is_ok());
        let new_cell = result.unwrap();
        assert_ne!(new_cell.spin, spin);
    }

    #[test]
    fn test_reproduce() {
        let mut env = make_env_for_division();
        env.spawn_rect_cell(
            Rect::new(Pos::new(30, 30), Pos::new(32, 32)),
            Cell::new_empty(4, 8, MockGenome::new(0))
        );

        let divided_spins = env.reproduce();
        assert_eq!(divided_spins.len(), 1);
    }

    #[test]
    fn test_reproduce_limit_population() {
        let mut env = make_env_for_division();
        env.max_cells = 1;
        env.spawn_rect_cell(
            Rect::new(Pos::new(30, 30), Pos::new(32, 32)), 
            Cell::new_empty(4, 8, MockGenome::new(0))
        );

        let divided_spins = env.reproduce();
        assert_eq!(divided_spins.len(), 0);
    }
}

