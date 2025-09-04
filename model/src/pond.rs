use crate::cell::Cell;
use crate::cellular_automata::CellularAutomata;
use crate::chem_space::ChemEnvironment;
use crate::genetics::genome::Genome;
use cellulars_lib::adhesion::ClonalAdhesion;
use cellulars_lib::cellular::{Cellular, RelCell};
use cellulars_lib::constants::Spin;
use cellulars_lib::environment::DivisionError;
use cellulars_lib::evolution::selector::Fit;
use cellulars_lib::lattice_entity::LatticeEntity;
use cellulars_lib::lattice_entity::LatticeEntity::Medium;
use cellulars_lib::positional::boundary::Boundary;
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::positional::rect::Rect;
use rand::Rng;
use rand_xoshiro::Xoshiro256StarStar;

// TODO: this struct can be made general if CellularAutomata is also general
pub struct Pond {
    pub env: ChemEnvironment,
    pub ca: CellularAutomata<ClonalAdhesion>,
    pub rng: Xoshiro256StarStar,
    pub update_period: u32,
    pub cell_search_radius: f32,
    pub cell_target_area: u32,
    pub division_enabled: bool,
    pub max_cells: u32,
    population_exploded: bool,
    time_step: u32,
}

impl Pond {
    pub fn new(
        env: ChemEnvironment,
        ca: CellularAutomata<ClonalAdhesion>,
        rng: Xoshiro256StarStar,
        update_period: u32,
        cell_search_radius: f32,
        cell_target_area: u32,
        division_enabled: bool,
        max_cells: u32
    ) -> Self {
        Self {
            env,
            ca,
            rng,
            update_period,
            cell_search_radius,
            cell_target_area,
            division_enabled,
            max_cells,
            population_exploded: false,
            time_step: 0
        }
    }
    
    pub fn step(&mut self) {
        self.ca.step(&mut self.env, &mut self.rng);
        if self.time_step % self.update_period == 0 {
            self.env.cells.iter_mut().for_each(|cell| cell.update());
            let new_spins = self.reproduce();
            for spin in new_spins {
                self.ca.adhesion.update_clones(spin, self.cell_search_radius, &self.env);
                // We could also instead choose to mutate at a fix rate throughout the cell's life cycle
                if let LatticeEntity::SomeCell(cell) = self.env.cells.get_entity_mut(spin) {
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
        empty_cell: Cell,
        rng: &mut impl Rng,
    ) -> Option<&RelCell<Cell>> {
        let cell_side = (cell_area as f32).sqrt() as usize;
        let pos = self.env.space.cell_lattice.random_pos(rng);
        self.spawn_rect_cell(
            Rect::new(
                pos,
                (pos.x + cell_side, pos.y + cell_side).into()
            ),
            empty_cell.clone()
        )
    }

    pub fn spawn_rect_cell(&mut self, rect: Rect<usize>, mut empty_cell: Cell) -> Option<&RelCell<Cell>> {
        if !self.can_add_cell() {
            return None;
        }

        let spin = self.env.cells.next_spin();
        for pos in rect.iter_positions() {
            if let Some(valid_pos) = self.env.space.lat_bound.valid_pos(pos.to_isize()) {
                let lat_pos = valid_pos.to_usize();
                if self.env.space.cell_lattice[lat_pos] != Medium.discriminant() {
                    continue
                }
                self.env.space.cell_lattice[lat_pos] = spin;
                self.env.update_edges(lat_pos);
                empty_cell.shift_position(
                    lat_pos,
                    true,
                    &self.env.space.bound
                );
                empty_cell.shift_chem(
                    lat_pos,
                    self.env.space.chem_lattice[lat_pos] as f32,
                    true,
                    &self.env.space.bound
                )
            }
        }
        if empty_cell.area() == 0 {
            return None;
        }
        Some(self.env.cells.push(empty_cell, None))
    }

    pub fn kill_cell(&mut self, cell: &mut RelCell<Cell>) {
        for pos in self.env.search_cell_box(cell, self.cell_search_radius) {
            // TODO!: Parameterize chance of medium
            if self.rng.random::<f32>() < 0.1 {
                self.env.space.cell_lattice[pos] = Medium.discriminant();
            }
        }
        for i in 0..self.ca.adhesion.clone_pairs.length() {
            self.ca.adhesion.clone_pairs[(cell.spin as usize, i)] = false
        }
        cell.apoptosis();
    }
    
    pub fn wipe_out(&mut self) {
        self.env.cells.wipe_out();
        self.env.space.cell_lattice.iter_values_mut().for_each(|value| {
            if *value >= LatticeEntity::first_cell_spin() {
                *value = Medium.discriminant();
            }
        });
        self.ca.adhesion.clone_pairs.clear();
    }

    // We take spin here because this operation is not safe with &Cell (pushing to vec can cause reallocation)
    pub fn divide_cell(&mut self, spin: Spin) -> Result<&RelCell<Cell>, DivisionError> {
        let new_spin = self.env.cells.next_spin();
        let cell_target_area = self.cell_target_area;
        let mut mom_clone = self
            .env
            .cells
            .get_entity_mut(spin)
            .expect_cell(&format!("passed non-cell with spin {spin} to `divide_cel()`"))
            .clone();

        let mut new_cell = mom_clone.birth();
        new_cell.set_target_area(self.cell_target_area);
        let new_positions: Vec<_> = self
            .env
            .search_cell_box(&mom_clone, self.cell_search_radius)
            .into_iter()
            .filter(|pos| {
                // TODO!: use principal component to determine division axis
                //  current algorithm hands out all x positions to the right of the cell centre to the new cell
                self.env.space.bound.displacement(Pos::new(pos.x as f32, pos.y as f32), mom_clone.center()).0 > 0.
            })
            .collect();
        for pos in new_positions {
            if mom_clone.area() == 1 {
                return Err(DivisionError::NewCellTooBig);
            }
            // TODO: this is basically the same as executing a lattice copy, unify the APIs
            //  This can happen when we move functions that change the env to Pond and CA
            self.env.space.cell_lattice[pos] = new_spin;
            let chem_at = self.env.space.chem_lattice[pos] as f32;
            new_cell.shift_position(
                pos,
                true,
                &self.env.space.bound
            );
            new_cell.shift_chem(
                pos,
                chem_at,
                true,
                &self.env.space.bound
            );
            mom_clone.shift_position(
                pos,
                false,
                &self.env.space.bound
            );
            mom_clone.shift_chem(
                pos,
                chem_at,
                false,
                &self.env.space.bound
            );
        }
        if new_cell.area() > 0 {
            let mom_spin = mom_clone.spin;
            mom_clone.set_target_area(cell_target_area);
            self.env.cells.replace(mom_clone);
            Ok(self.env.cells.push(new_cell, Some(mom_spin)))
        } else {
            Err(DivisionError::NewCellTooSmall)
        }
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
        
        for cell in self.env.cells.iter() {
            if cell.is_dying() {
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
            match self.divide_cell(spin) {
                Err(e) => {
                    log::warn!("Failed to divide cell with spin {spin} with error `{e:?}`");
                    None
                },
                Ok(cell) => Some(cell.spin)
            }
        }).collect()
    }

    pub fn make_chem_gradient(&mut self) {
        for row in 0..self.env.height() {
            for col in 0..self.env.width() {
                self.env.space.chem_lattice[(col, row).into()] = row.try_into().expect("lattice is too big");
            }
        }
    }

    pub fn can_add_cell(&mut self) -> bool {
        if self.env.cells.n_alive() < self.max_cells {
            return true;
        }
        if !self.population_exploded {
            log::warn!(
                        "Population exceeded maximum threshold `max-cells={}` during cell division",
                        {self.max_cells}
                    );
            log::warn!("This warning will be suppressed from now on");
            self.population_exploded = true;
        }
        false
    }
}

impl Fit for Pond {
    fn fitness(&self) -> f32 {
        let tot_fit: f32 = self
            .env
            .cells
            .iter()
            .filter(|cell| cell.is_alive())
            .map(|c| { c.fitness() })
            .sum();
        tot_fit / self.env.cells.n_alive() as f32
    }
}