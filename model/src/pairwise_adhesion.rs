use rand::Rng;
use cellulars_lib::constants::FloatType;
use cellulars_lib::prelude::{AdhesionSystem, CellIndex, Spin, SymmetricTable};

pub struct PairwiseAdhesion {
    medium_energy: FloatType,
    solid_energy: FloatType,
    adh_table: SymmetricTable<FloatType>,
}

impl PairwiseAdhesion {
    fn new(medium_energy: FloatType, solid_energy: FloatType, max_cells: CellIndex) -> Self {
        Self {
            medium_energy,
            solid_energy,
            adh_table: SymmetricTable::new(max_cells as usize)
        }
    }

    fn randomize_cell_energies(&mut self, min_energy: FloatType, max_energy: FloatType, rng: &mut impl Rng) {
        for (i, j) in self.adh_table.iter_index_pairs(None, None) {
            self.adh_table[(i, j)] = rng.random_range(min_energy..=max_energy);
        }
    }
}

impl<C> AdhesionSystem<C> for PairwiseAdhesion {
    fn adhesion_energy(&self, spin1: Spin, spin2: Spin, context: &C) -> FloatType {
        match (spin1, spin2) {
            (Spin::Some(c1), Spin::Some(c2)) => 2. * self.adh_table[(c1 as usize, c2 as usize)],
            (Spin::Medium, _) | (_, Spin::Medium) => self.medium_energy,
            (Spin::Solid, _)  | (_, Spin::Solid) => self.solid_energy,
        }
    }
}