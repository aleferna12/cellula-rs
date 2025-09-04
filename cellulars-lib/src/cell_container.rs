use crate::cellular::{Cellular, RelCell};
use crate::constants::Spin;
use crate::lattice_entity::LatticeEntity;
use crate::lattice_entity::LatticeEntity::{Medium, Solid, SomeCell};

pub struct CellContainer<C> {
    vec: Vec<RelCell<C>>
}

impl<C> CellContainer<C> {
    pub fn new() -> Self {
        Self {
            vec: vec![],
        }
    }
    
    pub fn n_cells(&self) -> Spin {
        self.vec.len().try_into().expect("there are more cells than supported by the type `Spin`")
    }
    
    /// Replaces the cell at `cell.spin`.
    pub fn replace(&mut self, cell: RelCell<C>) -> RelCell<C> {
        let index = cell.spin - LatticeEntity::first_cell_spin();
        std::mem::replace(&mut self.vec[index as usize], cell)
    }

    pub fn get_entity(&self, spin: Spin) -> LatticeEntity<&RelCell<C>> {
        if spin == Medium.discriminant() {
            return Medium;
        }
        if spin == Solid.discriminant() {
            return Solid;
        }
        SomeCell(&self.vec[(spin - LatticeEntity::first_cell_spin()) as usize])
    }

    pub fn get_entity_mut(&mut self, spin: Spin) -> LatticeEntity<&mut RelCell<C>> {
        if spin == Medium.discriminant() {
            return Medium;
        }
        if spin == Solid.discriminant() {
            return Solid;
        }
        SomeCell(&mut self.vec[(spin - LatticeEntity::first_cell_spin()) as usize])
    }

    pub fn wipe_out(&mut self) {
        self.vec.clear()
    }

    pub fn iter(&self) -> impl Iterator<Item=&RelCell<C>> {
        self.vec.iter()
    }

    pub fn iter_mut(&mut self) -> impl Iterator<Item=&mut RelCell<C>> {
        self.vec.iter_mut()
    }
}

impl<C: Cellular> CellContainer<C> {
    pub fn n_alive(&self) -> Spin {
        self.vec
            .iter()
            .filter(|cell| cell.is_alive())
            .count() as Spin
    }

    pub fn next_spin(&self) -> Spin {
        self.vec
            .iter()
            .find(|cell| !cell.is_alive())
            .map(|cell| cell.spin)
            .unwrap_or(self.n_cells() + LatticeEntity::first_cell_spin())
    }

    pub fn push(&mut self, cell: C, mom_spin: Option<Spin>) -> &RelCell<C> {
        let new_spin = self.next_spin();
        let index = new_spin - LatticeEntity::first_cell_spin();
        let rel_cell = RelCell {
            spin: new_spin,
            mom: mom_spin.unwrap_or(new_spin),
            cell
        };

        if index == self.n_cells() {
            self.vec.push(rel_cell);
        } else {
            self.replace(rel_cell);
        }
        &self.vec[index as usize]
    }
}

impl<C> Default for CellContainer<C> {
    fn default() -> Self {
        Self::new()
    }
}