use crate::cell::{Cell, RelCell};
use crate::constants::Spin;
use crate::environment::LatticeEntity;
use crate::environment::LatticeEntity::{Medium, Solid, SomeCell};
use crate::genome::{CellType, Genome};

pub struct CellContainer<G> {
    pub target_area: u32,
    pub div_area: u32,
    pub divide: bool,
    pub migrate: bool,
    vec: Vec<RelCell<G>>
}

impl<G> CellContainer<G> {
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

    pub fn push(&mut self, cell: Cell<G>, mom_spin: Option<Spin>) -> &RelCell<G> {
        let new_spin = self.next_spin();
        self.vec.push(RelCell {
            spin: new_spin,
            mom: mom_spin.unwrap_or(new_spin),
            cell
        });
        self.vec.last().unwrap()
    }
    
    /// Replaces the cell at `cell.spin`.
    pub fn replace(&mut self, cell: RelCell<G>) {
        let index = cell.spin - LatticeEntity::first_cell_spin();
        self.vec[index as usize] = cell
    }

    pub fn get_entity(&self, spin: Spin) -> LatticeEntity<&RelCell<G>> {
        if spin == Medium.discriminant() {
            return Medium;
        }
        if spin == Solid.discriminant() {
            return Solid;
        }
        SomeCell(&self.vec[(spin - LatticeEntity::first_cell_spin()) as usize])
    }

    pub fn get_entity_mut(&mut self, spin: Spin) -> LatticeEntity<&mut RelCell<G>> {
        if spin == Medium.discriminant() {
            return Medium;
        }
        if spin == Solid.discriminant() {
            return Solid;
        }
        SomeCell(&mut self.vec[(spin - LatticeEntity::first_cell_spin()) as usize])
    }

    pub fn update_cells(&mut self)
    where G: Genome {
        for cell in &mut self.vec {
            if let CellType::Divide = cell.cell_type && cell.target_area < self.div_area {
                cell.target_area += 1;
            }
            let chem_signal = cell.chem_mass;
            cell.genome.update_expression(chem_signal);
            cell.cell_type = cell.genome.get_cell_type();
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