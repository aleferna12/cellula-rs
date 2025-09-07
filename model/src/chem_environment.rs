use std::ops::{Deref, DerefMut};
use cellulars_lib::basic_cell::Cellular;
use cellulars_lib::cell_container::CellContainer;
use cellulars_lib::constants::Spin;
use crate::cell::Cell;
use crate::constants::BoundaryType;
use cellulars_lib::environment::{EdgesUpdate, Environment, Habitable};
use cellulars_lib::lattice::Lattice;
use cellulars_lib::lattice_entity::LatticeEntity::SomeCell;
use cellulars_lib::positional::boundaries::Boundary;
use cellulars_lib::positional::neighbourhood::MooreNeighbourhood;
use cellulars_lib::positional::pos::Pos;

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
        let valid_pos = match self.bounds.lattice_boundary.valid_pos(pos.to_isize()) {
            None => return EdgesUpdate { added: 0, removed: 0 },
            Some(pos_isize) => { pos_isize.to_usize() }
        };
        // TODO! chem should always be u32
        let chem_at_pos = self.chem_lattice[valid_pos] as f32;
        if let SomeCell(to_cell) = self.env.cells.get_entity_mut(to) {
            to_cell.shift_position(valid_pos, true, &self.env.bounds.boundary);
            to_cell.shift_chem(valid_pos, chem_at_pos, true, &self.env.bounds.boundary);
        }
        let from = self.cell_lattice[valid_pos];
        if let SomeCell(from_cell) = self.env.cells.get_entity_mut(from) {
            from_cell.shift_position(valid_pos, false, &self.env.bounds.boundary);
            from_cell.shift_chem(valid_pos, chem_at_pos, false, &self.env.bounds.boundary);
        }
        // Executes the copy
        self.cell_lattice[valid_pos] = to;
        self.update_edges(valid_pos)
    }
}