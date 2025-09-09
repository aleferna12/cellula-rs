use bon::Builder;
use crate::cell::Cell;
use crate::cellular_automata::CellularAutomata;
use crate::chem_environment::ChemEnvironment;
use crate::clonal_adhesion::ClonalAdhesion;
use crate::genetics::genome::Genome;
use cellulars_lib::basic_cell::{Alive, Cellular, RelCell};
use cellulars_lib::constants::Spin;
use cellulars_lib::environment::Habitable;
use cellulars_lib::evolution::selector::Fit;
use cellulars_lib::lattice_entity::LatticeEntity;
use cellulars_lib::lattice_entity::LatticeEntity::Medium;
use cellulars_lib::positional::boundaries::Boundary;
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::positional::rect::Rect;
use rand::Rng;
use rand_xoshiro::Xoshiro256StarStar;

// TODO: this struct can be made general if CellularAutomata is also general
#[derive(Clone, Builder)]
pub struct Pond {
    pub env: ChemEnvironment,
    pub ca: CellularAutomata<ClonalAdhesion>,
    pub rng: Xoshiro256StarStar,
    pub update_period: u32,
    pub cell_target_area: u32,
    pub division_enabled: bool,
    pub max_cells: u32,
    pub cell_search_scaler: f32,
    #[builder(default = false)]
    population_exploded: bool,
    #[builder(default = 0)]
    time_step: u32,
}

impl Pond {
    pub fn step(&mut self) {
        self.ca.step(&mut self.env, &mut self.rng);
        if self.time_step % self.update_period == 0 {
            self.env.cells.iter_mut().for_each(|cell| cell.update());
            let new_spins = self.reproduce();
            for spin in new_spins {
                self.ca.adhesion.update_clones(spin, &self.env);
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

    // TODO: make spawn as a circle with center at pos
    pub fn spawn_cell_random(
        &mut self,
        empty_cell: Cell,
        cell_area: u32,
    ) -> &RelCell<Cell> {
        let pos_isize = self.env.cell_lattice.random_pos(&mut self.rng).to_isize();
        let cell_side = ((cell_area as f32).sqrt() / 2.) as isize;
        let rect = Rect::new(
            Pos::new(pos_isize.x - cell_side, pos_isize.y - cell_side),
            Pos::new(pos_isize.x + cell_side, pos_isize.y + cell_side)
        );
        let positions = rect
            .iter_positions()
            .filter_map(|pos| self.env.bounds.lattice_boundary.valid_pos(pos))
            .map(|pos| pos.to_usize())
            .collect::<Vec<_>>();
        self.env.spawn_cell(
            empty_cell,
            positions
        )
    }

    pub fn divide_cell(&mut self, mom_spin: Spin) -> &RelCell<Cell> {
        let mom = self
            .env
            .cells
            .get_entity(mom_spin)
            .expect_cell("retrieved non-cell during cell division");
        let new_positions: Vec<_> = self
            .env
            .search_cell_box(mom, self.cell_search_scaler)
            .into_iter()
            .filter(|pos| {
                // TODO!: use principal component to determine division axis
                //  current algorithm hands out all x positions to the right of the cell centre to the new cell
                self.env.bounds.boundary.displacement(Pos::new(pos.x as f32, pos.y as f32), mom.center()).0 > 0.
            })
            .collect();

        let newborn = mom.birth();
        let new_spin = self.env.cells.push(newborn, Some(mom_spin)).spin;
        for pos in new_positions {
            self.env.grant_position(
                pos,
                new_spin,
            );
        }
        self.env.cells.get_entity(new_spin).expect_cell("retrieved non-cell during cell division")
    }

    pub fn kill_cell(&mut self, cell: &mut RelCell<Cell>) {
        for pos in self.env.search_cell_box(cell, self.cell_search_scaler) {
            // TODO!: Parameterize chance of medium
            if self.rng.random::<f32>() < 0.1 {
                self.env.cell_lattice[pos] = Medium.discriminant();
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
                .cells
                .get_entity(spin)
                .expect_cell("retrieved non-cell during reproduction");
            Some(self.divide_cell(mom.spin).spin)
        }).collect()
    }

    pub fn wipe_out(&mut self) {
        self.env.wipe_out();
        self.ca.adhesion.clone_pairs.clear();
    }

    pub fn can_add_cell(&mut self) -> bool {
        if self.env.cells.n_valid() < self.max_cells {
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
            .cells()
            .iter()
            .filter(|cell| cell.is_valid())
            .map(|c| { c.fitness() })
            .sum();
        tot_fit / self.env.cells().n_valid() as f32
    }
}