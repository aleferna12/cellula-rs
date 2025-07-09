use crate::cell::{RelCell, Cell};
use crate::constants::Spin;
use crate::environment::LatticeEntity;
use crate::environment::LatticeEntity::{Medium, Solid, SomeCell};
use crate::parameters::CellParameters;

pub struct CellContainer {
    pub target_area: u32,
    pub div_area: u32,
    pub divide: bool,
    pub migrate: bool,
    vec: Vec<RelCell>
}

impl CellContainer {
    pub fn n_cells(&self) -> Spin {
        self.vec.len().try_into().expect("there are more cells than supported by the type `Spin`")
    }
    
    pub fn next_spin(&self) -> Spin {
        self.n_cells() as Spin + LatticeEntity::first_cell_spin()
    }

    pub(crate) fn push(&mut self, cell: Cell, mom_spin: Option<Spin>) -> &RelCell {
        let new_spin = self.next_spin();
        self.vec.push(RelCell {
            spin: new_spin,
            mom: mom_spin.unwrap_or(new_spin),
            cell
        });
        self.vec.last().unwrap()
    }
    
    /// Replaces the cell at `cell.spin`.
    pub(crate) fn replace(&mut self, cell: RelCell) {
        let index = cell.spin - LatticeEntity::first_cell_spin();
        self.vec[index as usize] = cell
    }

    pub fn get_entity(&self, spin: Spin) -> LatticeEntity<&RelCell> {
        if spin == Medium::<&RelCell>.spin() {
            return Medium;
        }
        if spin == Solid::<&RelCell>.spin() {
            return Solid;
        }
        SomeCell(&self.vec[(spin - LatticeEntity::first_cell_spin()) as usize])
    }

    pub fn get_entity_mut(&mut self, spin: Spin) -> LatticeEntity<&mut RelCell> {
        if spin == Medium::<&RelCell>.spin() {
            return Medium;
        }
        if spin == Solid::<&RelCell>.spin() {
            return Solid;
        }
        SomeCell(&mut self.vec[(spin - LatticeEntity::first_cell_spin()) as usize])
    }

    pub fn update_cells(&mut self) {
        for cell in &mut self.vec {
            if cell.target_area < self.div_area {
                cell.target_area += 1;
            }
        }
    }
}

impl From<CellParameters> for CellContainer {
    fn from(params: CellParameters) -> Self {
        Self {
            target_area: params.target_area,
            div_area: params.div_area,
            divide: params.divide,
            migrate: params.migrate,
            vec: vec![],
        }
    }
}

impl<'a> IntoIterator for &'a CellContainer {
    type Item = &'a RelCell;
    type IntoIter = std::slice::Iter<'a, RelCell>;

    fn into_iter(self) -> Self::IntoIter {
        self.vec.iter()
    }
}

impl<'a> IntoIterator for &'a mut CellContainer {
    type Item = &'a mut RelCell;
    type IntoIter = std::slice::IterMut<'a, RelCell>;

    fn into_iter(self) -> Self::IntoIter {
        self.vec.iter_mut()
    }
}