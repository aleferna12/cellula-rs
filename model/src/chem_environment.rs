use crate::cell::Cell;
use crate::constants::BoundaryType;
use cellulars_lib::basic_cell::Cellular;
use cellulars_lib::cell_container::CellContainer;
use cellulars_lib::constants::Spin;
use cellulars_lib::environment::{EdgesUpdate, Environment, Habitable};
use cellulars_lib::lattice::Lattice;
use cellulars_lib::lattice_entity::LatticeEntity::SomeCell;
use cellulars_lib::positional::neighbourhood::MooreNeighbourhood;
use cellulars_lib::positional::pos::Pos;
use std::ops::{Deref, DerefMut};

#[derive(Clone)]
pub struct ChemEnvironment {
    env: Environment<Cell, MooreNeighbourhood, BoundaryType>,
    pub(crate) chem_lattice: Lattice<u32>
}

impl ChemEnvironment {
    pub fn new(env: Environment<Cell, MooreNeighbourhood, BoundaryType>) -> Self {
        let mut env_ = Self {
            chem_lattice: env.cell_lattice.clone(),
            env
        };
        env_.make_chem_gradient();
        env_
    }

    pub fn make_chem_gradient(&mut self) {
        for row in 0..self.height() {
            for col in 0..self.width() {
                self.chem_lattice[(col, row).into()] = row.try_into().expect("lattice is too big");
            }
        }
    }
}

impl Deref for ChemEnvironment {
    type Target = Environment<Cell, MooreNeighbourhood, BoundaryType>;
    fn deref(&self) -> &Self::Target {
        &self.env
    }
}

impl DerefMut for ChemEnvironment {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.env
    }
}

impl Habitable for ChemEnvironment {
    type Cell = Cell;

    fn cells(&self) -> &CellContainer<Self::Cell> {
        self.env.cells()
    }

    fn cells_mut(&mut self) -> &mut CellContainer<Self::Cell> {
        self.env.cells_mut()
    }

    fn grant_position(
        &mut self,
        pos: Pos<usize>,
        to: Spin
    ) -> EdgesUpdate {
        // TODO! chem should always be u32
        let chem_at_pos = self.chem_lattice[pos] as f32;
        if let SomeCell(to_cell) = self.env.cells.get_entity_mut(to) {
            to_cell.shift_position(pos, true, &self.env.bounds.boundary);
            to_cell.shift_chem(pos, chem_at_pos, true, &self.env.bounds.boundary);
        }
        let from = self.cell_lattice[pos];
        if let SomeCell(from_cell) = self.env.cells.get_entity_mut(from) {
            from_cell.shift_position(pos, false, &self.env.bounds.boundary);
            from_cell.shift_chem(pos, chem_at_pos, false, &self.env.bounds.boundary);
        }
        // Executes the copy
        self.cell_lattice[pos] = to;
        self.update_edges(pos)
    }
}