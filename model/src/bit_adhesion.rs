use crate::my_environment::MyEnvironment;
use cellulars_lib::adhesion::{AdhesionSystem, StaticAdhesion};
use cellulars_lib::spin::Spin;

#[derive(Clone)]
pub struct BitAdhesion {
    pub static_adhesion: StaticAdhesion
}

impl BitAdhesion {
    pub fn complementarity(lig: u64, rec: u64) -> u8 {
        (lig ^ rec).count_ones() as u8
    }
}

impl AdhesionSystem for BitAdhesion {
    type Context = MyEnvironment;

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
                    let gen1 = &context.cells.get_cell(c1).genome;
                    let gen2 = &context.cells.get_cell(c2).genome;
                    let avg_energy = (gen1.length + gen2.length) as f32 / 2.;
                    let complementarity = Self::complementarity(gen1.ligands(), gen2.receptors())
                        + Self::complementarity(gen2.ligands(), gen1.receptors());
                    // TODO!: Test
                    2. * self.static_adhesion.cell_energy + avg_energy - complementarity as f32
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
    use super::*;

    #[test]
    fn test_complementarity() {
        assert_eq!(0, BitAdhesion::complementarity(0, 0));
        assert_eq!(0, BitAdhesion::complementarity(255, 255));
        assert_eq!(8, BitAdhesion::complementarity(255, 0));
        assert_eq!(5, BitAdhesion::complementarity(250, 103));
    }
}
