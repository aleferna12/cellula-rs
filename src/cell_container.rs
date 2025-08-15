use crate::cell::{Cellular, RelCell};
use crate::constants::Spin;
use crate::environment::LatticeEntity;
use crate::environment::LatticeEntity::{Medium, Solid, SomeCell};

pub struct CellContainer<C> {
    pub target_area: u32,
    pub divide: bool,
    pub migrate: bool,
    vec: Vec<RelCell<C>>
}

impl<C> CellContainer<C> {
    pub fn new(
        target_area: u32,
        divide: bool,
        migrate: bool
    ) -> Self {
        Self {
            target_area,
            divide,
            migrate,
            vec: vec![],
        }
    }
    
    pub fn n_cells(&self) -> Spin {
        self.vec.len().try_into().expect("there are more cells than supported by the type `Spin`")
    }
    
    // TODO!: Reuse first free spin
    pub fn next_spin(&self) -> Spin {
        self.n_cells() + LatticeEntity::first_cell_spin()
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
    pub fn update_cells(&mut self) {
        for cell in &mut self.vec {
            cell.update();
        }
    }

    pub fn n_valid(&self) -> Spin {
        self.vec
            .iter()
            .filter(|cell| cell.is_valid())
            .count() as Spin
    }
}