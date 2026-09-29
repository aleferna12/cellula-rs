use rand::Rng;
use cellulars_lib::constants::FloatType;
use cellulars_lib::prelude::{AdhesionSystem, CellIndex, Spin, SymmetricTable};

#[derive(Clone)]
pub struct PairwiseAdhesion {
    medium_energy: FloatType,
    solid_energy: FloatType,
    adh_table: SymmetricTable<FloatType>,
}

impl PairwiseAdhesion {
    pub fn new(medium_energy: FloatType, solid_energy: FloatType, max_cells: CellIndex) -> Self {
        Self {
            medium_energy,
            solid_energy,
            adh_table: SymmetricTable::new(max_cells as usize)
        }
    }

    pub fn randomize_cell_energies(&mut self, min_energy: FloatType, max_energy: FloatType, rng: &mut impl Rng) {
        for (i, j) in self.adh_table.iter_index_pairs(None, None) {
            self.adh_table[(i, j)] = rng.random_range(min_energy..=max_energy);
        }
    }
}

impl<C> AdhesionSystem<C> for PairwiseAdhesion {
    fn adhesion_energy(&self, spin1: Spin, spin2: Spin, _context: &C) -> FloatType {
        match (spin1, spin2) {
            (Spin::Some(c1), Spin::Some(c2)) => {
                if c1 == c2 {
                    0.
                } else {
                    2. * self.adh_table[(c1 as usize, c2 as usize)]
                }
            },
            (Spin::Some(_), Spin::Medium) | (Spin::Medium, Spin::Some(_)) => self.medium_energy,
            (Spin::Some(_), Spin::Solid) | (Spin::Solid, Spin::Some(_)) => self.solid_energy,
            _ => 0.
        }
    }
}