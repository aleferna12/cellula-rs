use crate::cell::Cell;
use crate::constants::Spin;
use crate::environment::LatticeEntity;
use crate::environment::LatticeEntity::{Medium, Solid, SomeCell};

#[derive(Default)]
pub struct CellContainer {
    vec: Vec<Cell>
}

impl CellContainer {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn n_cells(&self) -> usize {
        self.vec.len()
    }
    
    pub fn next_spin(&self) -> Spin {
        self.n_cells() as Spin + LatticeEntity::first_cell_spin()
    }

    pub(crate) fn push(&mut self, cell: Cell) -> &Cell {
        assert_eq!(
            cell.spin, 
            self.next_spin(),
            "tried to add cell with incorrect spin {} (correct spin is {})", 
            cell.spin, 
            self.next_spin()
        );
        self.vec.push(cell);
        self.vec.last().unwrap()
    }
    
    /// Replaces the cell at `cell.spin`.
    pub(crate) fn replace(&mut self, cell: Cell) {
        let index = cell.spin - LatticeEntity::first_cell_spin();
        self.vec[index as usize] = cell
    }

    pub fn get_entity(&self, spin: Spin) -> LatticeEntity<&Cell> {
        if spin == Medium::<&Cell>.spin() {
            return Medium;
        }
        if spin == Solid::<&Cell>.spin() {
            return Solid;
        }
        SomeCell(&self.vec[(spin - LatticeEntity::first_cell_spin()) as usize])
    }

    pub fn get_entity_mut(&mut self, spin: Spin) -> LatticeEntity<&mut Cell> {
        if spin == Medium::<&Cell>.spin() {
            return Medium;
        }
        if spin == Solid::<&Cell>.spin() {
            return Solid;
        }
        SomeCell(&mut self.vec[(spin - LatticeEntity::first_cell_spin()) as usize])
    }

    // TODO: move to Cell?
    pub fn update_cells(&mut self, cell_div_area: u32, cells_grow: bool) {
        for cell in self {
            if cells_grow && cell.target_area < cell_div_area {
                cell.target_area += 1;
            }
        }
    }
}

impl<'a> IntoIterator for &'a CellContainer {
    type Item = &'a Cell;
    type IntoIter = std::slice::Iter<'a, Cell>;

    fn into_iter(self) -> Self::IntoIter {
        self.vec.iter()
    }
}

impl<'a> IntoIterator for &'a mut CellContainer {
    type Item = &'a mut Cell;
    type IntoIter = std::slice::IterMut<'a, Cell>;

    fn into_iter(self) -> Self::IntoIter {
        self.vec.iter_mut()
    }
}