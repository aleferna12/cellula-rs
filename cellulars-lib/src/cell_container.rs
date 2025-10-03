use crate::basic_cell::{Alive, Cellular, RelCell};
use crate::constants::CellIndex;

#[derive(Clone, Debug)]
pub struct CellContainer<C> {
    vec: Vec<RelCell<C>>
}

impl<C> CellContainer<C> {
    pub fn new() -> Self {
        Self {
            vec: vec![],
        }
    }
    
    pub fn n_cells(&self) -> CellIndex {
        self.vec.len().try_into().expect("there are more cells than supported by the type `Spin`")
    }
    
    /// Replaces the cell at `cell.index`.
    pub fn replace(&mut self, cell: RelCell<C>) -> RelCell<C> {
        std::mem::replace(&mut self.vec[cell.index as usize], cell)
    }

    pub fn get_cell(&self, index: CellIndex) -> &RelCell<C> {
        &self.vec[index as usize]
    }

    pub fn get_cell_mut(&mut self, index: CellIndex) -> &mut RelCell<C> {
        &mut self.vec[index as usize]
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
    pub fn n_valid(&self) -> CellIndex {
        self.vec
            .iter()
            .filter(|cell| cell.is_valid())
            .count() as CellIndex
    }
    
    pub fn n_alive(&self) -> CellIndex
    where C: Alive {
        self.vec
            .iter()
            .filter(|cell| cell.is_alive())
            .count() as CellIndex
    }

    pub fn next_index(&self) -> CellIndex {
        self.vec
            .iter()
            .find(|cell| !cell.is_valid())
            .map(|cell| cell.index)
            .unwrap_or(self.n_cells())
    }

    pub fn add(&mut self, cell: C) -> &mut RelCell<C> {
        let new_index = self.next_index();
        let rel_cell = RelCell {
            index: new_index,
            cell
        };

        if new_index == self.n_cells() {
            self.vec.push(rel_cell);
        } else {
            self.replace(rel_cell);
        }
        &mut self.vec[new_index as usize]
    }

    pub fn push(&mut self, cell: C) -> &mut RelCell<C> {
        let new_index = self.n_cells();
        let rel_cell = RelCell {
            index: new_index,
            cell
        };
        self.vec.push(rel_cell);
        &mut self.vec[new_index as usize]
    }
}

impl<C> Default for CellContainer<C> {
    fn default() -> Self {
        Self::new()
    }
}