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
pub struct BitAdhesion {
    pub gene_energy: f32,
    pub static_adhesion: StaticAdhesion
}

impl BitAdhesion {
    pub fn complementarity(lig: u64, rec: u64) -> u8 {
        (lig ^ rec).count_ones() as u8
    }
}

impl AdhesionSystem for BitAdhesion {
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
                    let gen1 = &context.get_cell(c1).genome;
                    let gen2 = &context.get_cell(c2).genome;
                    let complementarity = Self::complementarity(gen1.ligands(), gen2.receptors())
                        + Self::complementarity(gen2.ligands(), gen1.receptors());
                    2. * (
                        self.static_adhesion.cell_energy
                            + self.gene_energy
                            * (0.5 - complementarity as f32 / (gen1.length + gen2.length) as f32)
                    )
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
    use super::*;
    
    fn bit_adhesion() -> BitAdhesion {
        BitAdhesion {
            static_adhesion: StaticAdhesion {
                cell_energy: 16.,
                medium_energy: 16.,
                solid_energy: 16.,
            },
            gene_energy: 16.,
        }
    }

    #[test]
    fn test_complementarity() {
        assert_eq!(0, BitAdhesion::complementarity(0, 0));
        assert_eq!(0, BitAdhesion::complementarity(255, 255));
        assert_eq!(8, BitAdhesion::complementarity(255, 0));
        assert_eq!(5, BitAdhesion::complementarity(250, 103));
    }

    #[test]
    fn test_adhesion_energy() {
        let bit_adh = bit_adhesion();

        let recep_gamma = [
            (0b11110000, 0.),
            (0b00000000, 8.),
            (0b11111111, -8.),
        ];

        for (recep, gamma) in recep_gamma {
            let mut cells = CellContainer::new();
            let cell = Cell::new_empty(
                0,
                0,
                BitGenome::new(
                    0b11111111,
                    recep,
                    0.,
                    8
                ).unwrap(),
                false
            );
            let index1 = cells.push(cell.clone()).index;
            let index2 = cells.push(cell).index;
            let calc_gamma = bit_adh.static_adhesion.medium_energy - bit_adh.adhesion_energy(
                Spin::Some(index1),
                Spin::Some(index2),
                &cells
            ) / 2.;
            assert_eq!(gamma, calc_gamma);
        }
    }
    
    #[test]
    fn test_adhesion_random_genomes() {
        let bit_adh = bit_adhesion();
      
        let mut gammas = vec![];
        let mut rng = rand::rngs::StdRng::seed_from_u64(132415);
        for _ in 0..100_000 {
            let mut cells = CellContainer::new();
            let cell1 = Cell::new_empty(
                0,
                0,
                BitGenome::new_random(0., 8, &mut rng).unwrap(),
                false
            );
            let cell2 = Cell::new_empty(
                0,
                0,
                BitGenome::new_random(0., 8, &mut rng).unwrap(),
                false
            );
            let index1 = cells.push(cell1).index;
            let index2 = cells.push(cell2).index;
            let gamma = bit_adh.static_adhesion.medium_energy - bit_adh.adhesion_energy(
                Spin::Some(index1),
                Spin::Some(index2),
                &cells
            ) / 2.;
            gammas.push(gamma);
        }
        let sum: f32 = gammas.iter().sum();
        let mean = sum / gammas.len() as f32;
        assert!(mean < 0.05 && mean > -0.05);
        assert_eq!(8., gammas.iter().copied().fold(f32::NEG_INFINITY, f32::max));
        assert_eq!(-8., gammas.iter().copied().fold(f32::INFINITY, f32::min));
    }
}
