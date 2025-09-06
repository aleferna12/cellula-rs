use crate::cell::Cell;
use crate::cellular_automata::CellularAutomata;
use crate::chem_environment::ChemEnvironment;
use crate::genetics::genome::Genome;
use cellulars_lib::basic_cell::{Alive, Cellular, RelCell};
use cellulars_lib::constants::Spin;
use cellulars_lib::evolution::selector::Fit;
use cellulars_lib::lattice_entity::LatticeEntity;
use cellulars_lib::lattice_entity::LatticeEntity::Medium;
use cellulars_lib::positional::boundaries::{Boundary, PosValidator};
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::positional::rect::Rect;
use rand::Rng;
use rand_xoshiro::Xoshiro256StarStar;
use cellulars_lib::environment::Habitable;
use crate::clonal_adhesion::ClonalAdhesion;

// TODO: this struct can be made general if CellularAutomata is also general
pub struct Pond {
    pub env: ChemEnvironment,
    pub ca: CellularAutomata<ClonalAdhesion>,
    pub rng: Xoshiro256StarStar,
    pub update_period: u32,
    pub cell_target_area: u32,
    pub division_enabled: bool,
    pub max_cells: u32,
    pub cell_search_scaler: f32,
    population_exploded: bool,
    time_step: u32,
}

impl Pond {
    pub fn new(
        env: ChemEnvironment,
        ca: CellularAutomata<ClonalAdhesion>,
        rng: Xoshiro256StarStar,
        update_period: u32,
        cell_target_area: u32,
        cell_search_scaler: f32,
        division_enabled: bool,
        max_cells: u32
    ) -> Self {
        Self {
            env,
            ca,
            rng,
            update_period,
            cell_target_area,
            cell_search_scaler,
            division_enabled,
            max_cells,
            population_exploded: false,
            time_step: 0
        }
    }
    
    pub fn step(&mut self) {
        self.ca.step(&mut self.env, &mut self.rng);
        if self.time_step % self.update_period == 0 {
            self.env.cells().iter_mut().for_each(|cell| cell.update());
            let new_spins = self.reproduce();
            for spin in new_spins {
                self.ca.adhesion.update_clones(spin, &self.env);
                // We could also instead choose to mutate at a fix rate throughout the cell's life cycle
                if let LatticeEntity::SomeCell(cell) = self.env.cells().get_entity_mut(spin) {
                    cell.genome.attempt_mutate(&mut self.rng);
                } else { 
                    panic!("newborn is not a cell")
                }
            }
        }
        self.time_step += 1;
    }

    pub fn spawn_cell_random(
        &mut self,
        cell_area: u32,
        empty_cell: Cell
    ) -> Option<&RelCell<Cell>> {
        let cell_side = (cell_area as f32).sqrt() as usize;
        let pos = self.env.cell_lattice().random_pos(&mut self.rng);
        self.spawn_cell_rect(
            &Rect::new(
                pos,
                (pos.x + cell_side, pos.y + cell_side).into()
            ),
            empty_cell
        )
    }

    pub fn spawn_cell_rect(
        &mut self,
        rect: &Rect<usize>,
        empty_cell: Cell
    ) -> Option<&RelCell<Cell>> {
        let valid_pos = rect
            .iter_positions()
            .filter_map(|pos| {
                match self.env.lattice_boundary().valid_pos(pos.to_isize()) { 
                    Some(pos) => Some(pos.to_usize()),
                    None => None
                }
            });
        self.env.spawn_cell(empty_cell, valid_pos, self.env.boundaries.boundary())
    }

    fn spawn_cell(
        &mut self,
        empty_cell: Cell,
        positions: impl IntoIterator<Item = Pos<usize>>
    ) -> Option<&RelCell<Cell>> {
        let valid_positions = positions
            .into_iter()
            .filter(|&pos| {
                self.env.cell_lattice()[pos] == Medium.discriminant()
            })
            .collect::<Vec<_>>();

        if valid_positions.is_empty() {
            return None;
        }

        let new_spin = self.env.cells_mut().push(empty_cell, None).spin;
        for pos in valid_positions {
            self.env.grant_position(pos, new_spin, self.env.boundary());
        }
        Some(self.env.cells().get_entity(new_spin).expect_cell("retrieved non-cell while spawning cell"))
    }

    pub fn divide_cell(&mut self, mother: &RelCell<C>) -> Result<&RelCell<C>, cellulars_lib::environment::DivisionError> {
        let new_positions: Vec<_> = self
            .search_cell_box(mother)
            .into_iter()
            .filter(|pos| {
                // TODO!: use principal component to determine division axis
                //  current algorithm hands out all x positions to the right of the cell centre to the new cell
                self.boundaries.boundary.displacement(Pos::new(pos.x as f32, pos.y as f32), mother.center()).0 > 0.
            })
            .collect();
        if new_positions.is_empty() {
            return Err(cellulars_lib::environment::DivisionError::NewCellTooSmall);
        }
        if new_positions.len() >= mother.area() as usize {
            return Err(cellulars_lib::environment::DivisionError::NewCellTooBig);
        }

        let new_spin = self.cells.push(mother.birth(), Some(mother.spin)).spin;
        for pos in new_positions {
            self.grant_position(
                pos,
                new_spin,
            );
        }
        Ok(self.cells.get_entity(new_spin).expect_cell("retrieved non-cell during cell division"))
    }

    pub fn kill_cell(&mut self, cell: &mut RelCell<C>)
    where C: Alive
    {
        for pos in self.search_cell_box(cell) {
            // TODO!: Parameterize chance of medium
            if self.rng.random::<f32>() < 0.1 {
                self.cell_lattice[pos] = Medium.discriminant();
            }
        }
        cell.apoptosis();
    }

    // With some unsafe code we can return Vec<&RelCell> from this function, but it would
    // require that self.divide_cell never invalidates any references to self.cells
    // we need thorough testing of self.divide_cells to make this change, and the performance
    // gain is minimal (although the ergonomic gains are significant)
    pub fn reproduce(&mut self) -> Vec<Spin> {
        let mut divide = vec![];
        if !self.division_enabled {
            return divide;
        }
        
        for cell in self.env.cells().iter() {
            if !cell.is_alive() {
                continue;
            }
            // Currently cells don't need to express the dividing type to divide, they just need to be big enough
            if cell.area() >= cell.divide_area() {
                divide.push(cell.spin);
            }
        }
        divide.into_iter().filter_map(|spin| {
            if !self.can_add_cell() {
                return None;
            }
            let mom = self.env
                .cells()
                .get_entity(spin)
                .expect_cell("retrieved non-cell during reproduction");
            match self.env.divide_cell(mom) {
                Err(e) => {
                    log::warn!("Failed to divide cell with spin {spin} with error `{e:?}`");
                    None
                },
                Ok(cell) => Some(cell.spin)
            }
        }).collect()
    }
}

impl Fit for Pond {
    fn fitness(&self) -> f32 {
        let tot_fit: f32 = self
            .env
            .cells()
            .iter()
            .filter(|cell| cell.is_valid())
            .map(|c| { c.fitness() })
            .sum();
        tot_fit / self.env.cells().n_valid() as f32
    }
}

#[derive(Debug)]
pub enum DivisionError {
    NewCellTooSmall,
    NewCellTooBig
}