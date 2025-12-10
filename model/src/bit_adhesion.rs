use crate::my_environment::MyEnvironment;
use cellulars_lib::adhesion::{AdhesionSystem, StaticAdhesion};
use cellulars_lib::spin::Spin;

#[derive(Clone)]
pub struct BitAdhesion {
    pub static_adhesion: StaticAdhesion
}

impl BitAdhesion {
    pub fn contact_energy(lig: u64, rec: u64) -> u8 {
        (lig ^ rec).count_ones() as u8
    }
}

impl AdhesionSystem for BitAdhesion {
    type Context = MyEnvironment;

    fn adhesion_energy(&self, spin1: Spin, spin2: Spin, context: &Self::Context) -> f32 {
        if cfg!(feature = "static_adhesion") {
            return self.static_adhesion.adhesion_energy(spin1, spin2, &());
        }
        match (spin1, spin2) {
            (Spin::Some(c1), Spin::Some(c2)) => {
                if c1 == c2 {
                    0.
                } else {
                    let gen1 = &context.cells.get_cell(c1).genome;
                    let gen2 = &context.cells.get_cell(c2).genome;
                    let energy = Self::contact_energy(gen1.ligands(), gen2.receptors()) as f32
                        + Self::contact_energy(gen2.ligands(), gen1.receptors()) as f32;
                    2. * self.static_adhesion.cell_energy + energy
                }
            }
            (Spin::Some(_), Spin::Medium) | (Spin::Medium, Spin::Some(_)) => self.static_adhesion.medium_energy,
            (Spin::Some(_), Spin::Solid) | (Spin::Solid, Spin::Some(_)) => self.static_adhesion.solid_energy,
            _ => 0.
        }
    }
}