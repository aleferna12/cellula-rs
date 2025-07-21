use crate::cell::{Cell, RelCell};
use crate::cell_container::CellContainer;
use crate::constants::{Spin, SIZE_SCALE};
use crate::environment::LatticeEntity::*;
use crate::io::parameters::{CellParameters, EnvironmentParameters};
use crate::positional::boundary::Boundary;
use crate::positional::edge::Edge;
use crate::positional::edge_book::EdgeBook;
use crate::positional::neighbourhood::{MooreNeighbourhood, Neighbourhood};
use crate::positional::pos::Pos;
use crate::positional::rect::Rect;
use crate::space::Space;
use rand::Rng;
use crate::genome::SpecialisedGrn;

pub struct Environment {
    pub space: Space,
    pub cells: CellContainer<SpecialisedGrn>,
    pub edge_book: EdgeBook,
    pub neighbourhood: MooreNeighbourhood,
    pub update_period: u32,
    pub cell_search_radius: f32,
    pub population_exploded: bool,
    pub max_cells: Spin
}

impl Environment {
    pub fn new(params: EnvironmentParameters, rng: &mut impl Rng) -> Self {
        let mut env = Self::new_empty(
            params.width,
            params.height,
            params.neigh_r,
            params.update_period,
            params.cell_search_radius,
            params.max_cells,
            params.cell
        );

        if params.enclose {
            log::info!("Enclosing environment with a border");
            env.make_border();
        }

        log::info!("Initialising light gradient");
        for row in 0..env.height() {
            for col in 0..env.width() {
                env.space.light_lattice[(col, row).into()] = row.try_into().expect("Lattice is too big");
            }
        }

        log::info!("Creating cells");
        let mut cell_count = 0;
        let cell_side = (params.cell_start_area as f32).sqrt() as usize;
        for _ in 0..params.starting_cells {
            let pos = env.space.cell_lattice.random_pos(rng);
            let cell = env.spawn_rect_cell(
                Rect::new(
                    pos,
                    (pos.x + cell_side, pos.y + cell_side).into()
                ),
            );
            if cell.is_some() {
                cell_count += 1;
            }
        }
        log::info!("Created {} out of the {} cells requested", cell_count, params.starting_cells);
        
        env
    }

    pub fn new_empty(
        width: usize,
        height: usize,
        neigh_r: u8,
        update_period: u32,
        cell_search_radius: f32,
        max_cells: Spin,
        cell_parameters: CellParameters
    ) -> Self {
        Self {
            space: Space::new(width, height),
            cells: CellContainer::from(cell_parameters),
            edge_book: EdgeBook::new(),
            neighbourhood: MooreNeighbourhood::new(neigh_r),
            max_cells,
            update_period,
            cell_search_radius,
            population_exploded: false
        }
    }

    /// Empty environment with arbitrary cell parameters for testing purposes.
    /// 
    /// Do not use this in production, no cells can be added to an environment created through this method.
    pub fn new_empty_test(width:usize, height: usize) -> Self {
        Environment::new_empty(
            width,
            height,
            1,
            0,
            0.,
            0,
            CellParameters {
                target_area: 0,
                div_area: 0,
                divide: false,
                migrate: false
            },
        )
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

    pub fn spawn_rect_cell(&mut self, rect: Rect<usize>) -> Option<&RelCell<SpecialisedGrn>> {
        let spin = self.cells.n_cells() as Spin + LatticeEntity::first_cell_spin();
        let mut cell = Cell::new(
            self.cells.target_area,
            SpecialisedGrn::new(1. / self.height() as f32, SIZE_SCALE)
        );
        
        for pos in rect.iter_positions() {
            if let Some(valid_pos) = self.space.lat_bound.valid_pos(pos.into()) {
                let lat_pos = valid_pos.into();
                if self.space.cell_lattice[lat_pos] != Medium.discriminant() {
                    continue
                }
                self.space.cell_lattice[lat_pos] = spin;
                self.update_edges(lat_pos);
                cell.shift_position(
                    lat_pos, 
                    self.space.light_lattice[lat_pos], 
                    true,
                    &self.space.bound
                );
            }
        }
        if cell.area == 0 { 
            return None;
        }
        self.cells.push(cell, None);
        Some(self.cells.get_entity(spin).unwrap_cell())
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
    
    // With some unsafe code we can return Vec<&RelCell> from this function, but it would
    // require that self.divide_cell never invalidates any references to self.cells
    // we need thorough testing of self.divide_cells to make this change, and the performance
    // gain is minimal (although the ergonomic gains are significant)
    pub fn reproduce(&mut self) -> Vec<Spin> {
        let mut divide = vec![];
        for cell in &self.cells {
            if self.cells.divide && cell.area >= self.cells.div_area {
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
    pub fn divide_cell(&mut self, spin: Spin) -> Result<&RelCell<SpecialisedGrn>, DivisionError> {
        let light_scale = 1. / self.height() as f32;
        let new_spin = self.cells.next_spin();
        let cell_target_area = self.cells.target_area;
        let mom_cell = self
            .cells
            .get_entity_mut(spin)
            .expect_cell(&format!("passed non-cell with spin {spin} to `divide_cel()`"));
        // We modify this mock cell to allow the division to be cancelled in the case of an error
        let mut mom_clone = mom_cell.clone();
        
        let mut new_cell = Cell::new(
            cell_target_area,
            SpecialisedGrn::new(light_scale, SIZE_SCALE)
        );
        let new_positions: Vec<_> = self
            .space
            .box_cell_positions(mom_cell, self.cell_search_radius)
            .into_iter()
            .filter(|pos| { 
                // TODO!: use principal component to determine division axis
                //  current algorithm hands out all x positions to the right of the cell centre to the new cell
                //  it also might be worth writing a faster implementation for FixedBoundary (LatticeBoundary)
                self.space.bound.displacement(Pos::new(pos.x as f32, pos.y as f32), mom_cell.center).0 > 0.
            })
            .collect();
        for pos in new_positions {
            if mom_cell.area == 1 {
                return Err(DivisionError::NewCellTooBig);
            }
            self.space.cell_lattice[pos] = new_spin;
            new_cell.shift_position(
                pos, 
                self.space.light_lattice[pos],
                true,
                &self.space.bound
            );
            mom_clone.shift_position(
                pos,
                self.space.light_lattice[pos],
                false,
                &self.space.bound
            );
        }
        if new_cell.area > 0 {
            let mom_spin = mom_cell.spin;
            mom_clone.target_area = cell_target_area;
            self.cells.replace(mom_clone);
            Ok(self.cells.push(new_cell, Some(mom_spin)))
        } else {
            Err(DivisionError::NewCellTooSmall)
        }
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

impl<G> LatticeEntity<&RelCell<G>> {
    /// Maps the `LatticeEntity` to a unique spin value.
    pub fn spin(&self) -> Spin {
        match self {
            SomeCell(cell) => cell.spin,
            Medium => Medium.discriminant(),
            Solid => Solid.discriminant()
        }
    }
}

impl<C: std::fmt::Debug> LatticeEntity<C> {
    pub fn unwrap_cell(self) -> C {
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
    use crate::positional::rect::Rect;
    use crate::positional::pos::Pos;

    fn make_env_for_division() -> Environment {
        let env = Environment::new_empty(
            100,
            100,
            1,
            1,
            2.0,
            10,
            CellParameters {
                target_area: 4,
                div_area: 4,
                divide: true,
                migrate: false,
            },
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
    fn test_make_border() {
        let mut env = Environment::new_empty_test(10, 5);
        env.make_border();

        for x in 0..10 {
            assert_eq!(env.space.cell_lattice[Pos::new(x, 0)], Solid.discriminant());
            assert_eq!(env.space.cell_lattice[Pos::new(x, 4)], Solid.discriminant());
        }

        for y in 1..5 {
            assert_eq!(env.space.cell_lattice[Pos::new(9, y)], Solid.discriminant());
        }

        for y in 1..4 {
            assert_eq!(env.space.cell_lattice[Pos::new(0, y)], Solid.discriminant());
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
        env.spawn_rect_cell(rect);

        let spin = LatticeEntity::first_cell_spin();
        let result = env.divide_cell(spin);
        assert!(result.is_ok());
        let new_cell = result.unwrap();
        assert_ne!(new_cell.spin, spin);
    }

    #[test]
    fn test_reproduce() {
        let mut env = make_env_for_division();
        env.spawn_rect_cell(Rect::new(Pos::new(30, 30), Pos::new(32, 32)));

        let divided_spins = env.reproduce();
        assert_eq!(divided_spins.len(), 1);
    }

    #[test]
    fn test_reproduce_limit_population() {
        let mut env = make_env_for_division();
        env.max_cells = 1;
        env.spawn_rect_cell(Rect::new(Pos::new(30, 30), Pos::new(32, 32)));

        let divided_spins = env.reproduce();
        assert_eq!(divided_spins.len(), 0);
    }
}

