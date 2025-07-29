use crate::cell::{CellLike, RelCell};
use crate::constants::Spin;
use crate::environment::LatticeEntity;
use crate::environment::LatticeEntity::{Medium, Solid, SomeCell};

pub struct CellContainer<C> {
    pub target_area: u32,
    pub div_area: u32,
    pub divide: bool,
    pub migrate: bool,
    vec: Vec<RelCell<C>>
}

impl<C> CellContainer<C> {
    pub fn new(
        target_area: u32,
        div_area: u32,
        divide: bool,
        migrate: bool
    ) -> Self {
        Self {
            target_area,
            div_area,
            divide,
            migrate,
            vec: vec![],
        }
    }
    
    pub fn n_cells(&self) -> Spin {
        self.vec.len().try_into().expect("there are more cells than supported by the type `Spin`")
    }
    
    pub fn next_spin(&self) -> Spin {
        self.n_cells() as Spin + LatticeEntity::first_cell_spin()
    }

    pub fn push(&mut self, cell: C, mom_spin: Option<Spin>) -> &RelCell<C> {
        let new_spin = self.next_spin();
        self.vec.push(RelCell {
            spin: new_spin,
            mom: mom_spin.unwrap_or(new_spin),
            cell
        });
        self.vec.last().unwrap()
    }
    
    /// Replaces the cell at `cell.spin`.
    pub fn replace(&mut self, cell: RelCell<C>) {
        let index = cell.spin - LatticeEntity::first_cell_spin();
        self.vec[index as usize] = cell
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

    pub fn update_cells(&mut self)
    where C: CellLike {
        for cell in &mut self.vec {
            cell.update();
        }
    }
}

impl<'a, G> IntoIterator for &'a CellContainer<G> {
    type Item = &'a RelCell<G>;
    type IntoIter = std::slice::Iter<'a, RelCell<G>>;

    fn into_iter(self) -> Self::IntoIter {
        self.vec.iter()
    }
}

impl<'a, G> IntoIterator for &'a mut CellContainer<G> {
    type Item = &'a mut RelCell<G>;
    type IntoIter = std::slice::IterMut<'a, RelCell<G>>;

    fn into_iter(self) -> Self::IntoIter {
        self.vec.iter_mut()
    }
}