use cellulars_lib::adhesion::{AdhesionSystem, StaticAdhesion};
use cellulars_lib::cell_container::CellContainer;
use cellulars_lib::spin::Spin;
use crate::cell::Cell;

/// Adhesion using cell bit ligands and receptors.
///
/// Adhesion strength is defined differently to the first paper.
/// When `static_adhesion.cell_energy == static_adhesion.medium_energy`,
/// ligands and receptors that are 50% complimentary give gamma = 0.
/// 0% complimentary give `-max_bit_energy / 2`.
/// 100% complimentary give `max_bit_energy / 2`.
#[derive(Clone)]
pub struct SimpleAdhesion {
    pub static_adhesion: StaticAdhesion
}

impl AdhesionSystem for SimpleAdhesion {
    type Context = CellContainer<Cell>;

    fn adhesion_energy(&self, spin1: Spin, spin2: Spin, context: &Self::Context) -> f32 {
        // This is a feature because it quite heavily affects performance
        if cfg!(feature = "static-adhesion") {
            return self.static_adhesion.adhesion_energy(spin1, spin2, &());
        }

        match (spin1, spin2) {
            (Spin::Some(c1), Spin::Some(c2)) => {
                if c1 == c2 {
                    0.
                } else {
                    let gen1 = &context.get_cell(c1).genome.val;
                    let gen2 = &context.get_cell(c2).genome.val;
                    2. * self.static_adhesion.cell_energy + gen1 + gen2
                }
            }
            (Spin::Some(_), Spin::Medium) | (Spin::Medium, Spin::Some(_)) => self.static_adhesion.medium_energy,
            (Spin::Some(_), Spin::Solid) | (Spin::Solid, Spin::Some(_)) => self.static_adhesion.solid_energy,
            _ => 0.
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::evolution::bit_genome::BitGenome;
    use rand::SeedableRng;
    use crate::cell::SimpleGenome;
    use super::*;
    
    fn bit_adhesion() -> SimpleAdhesion {
        SimpleAdhesion {
            static_adhesion: StaticAdhesion {
                cell_energy: 16.,
                medium_energy: 16.,
                solid_energy: 16.,
            },
        }
    }

    #[test]
    fn test_adhesion_energy() {
        let bit_adh = bit_adhesion();

        let mut cells = CellContainer::new();
        let cell = Cell::new_empty(
            0,
            0,
            SimpleGenome::new(0., 0.1, 0.1)
        );
        let index1 = cells.push(cell.clone()).index;
        let index2 = cells.push(cell).index;
        let calc_gamma = bit_adh.static_adhesion.medium_energy - bit_adh.adhesion_energy(
            Spin::Some(index1),
            Spin::Some(index2),
            &cells
        ) / 2.;
        assert_eq!(0., calc_gamma);
    }
}
