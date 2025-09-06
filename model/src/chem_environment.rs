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
    chem_lattice: Lattice<u32>
}

impl ChemEnvironment {
    pub fn new(env: Environment<Cell, MooreNeighbourhood, BoundaryType>) -> Self {
        Self {
            chem_lattice: env.cell_lattice().clone(),
            env
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

    fn cell_lattice(&self) -> &Lattice<Spin> {
        self.env.cell_lattice()
    }

    fn cell_lattice_mut(&mut self) -> &mut Lattice<Spin> {
        self.env.cell_lattice_mut()
    }

    fn update_edges(&mut self, pos: Pos<usize>) -> EdgesUpdate {
        self.env.update_edges(pos)
    }

    fn grant_position(
        &mut self,
        pos: Pos<usize>,
        to: Spin,
        boundary: &impl Boundary<Coord = f32>
    ) -> EdgesUpdate {
        // TODO: make chem a u32
        let chem_at = self.chem_lattice[pos] as f32;
        if let SomeCell(to_cell) = self.cells_mut().get_entity_mut(to) {
            to_cell.shift_position(pos, true, boundary);
            to_cell.shift_chem(pos, chem_at, true, boundary);
        }
        let from = self.cell_lattice()[pos];
        if let SomeCell(from_cell) = self.cells_mut().get_entity_mut(from) {
            from_cell.shift_position(pos, false, boundary);
            from_cell.shift_chem(pos, chem_at, false, boundary);
        }
        // Executes the copy
        self.cell_lattice_mut()[pos] = to;
        self.update_edges(pos)
    }
}