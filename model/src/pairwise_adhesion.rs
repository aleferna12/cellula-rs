use crate::my_environment::MyEnvironment;
use cellulars_lib::constants::FloatType;
use cellulars_lib::prelude::{AdhesionSystem, Spin};

#[derive(Clone)]
pub struct PairwiseAdhesion {
    pub medium_energy: FloatType,
    pub solid_energy: FloatType,
    pub adh_energy: FloatType,
    pub non_adh_energy: FloatType,
}

impl AdhesionSystem<MyEnvironment> for PairwiseAdhesion {
    fn adhesion_energy(&self, spin1: Spin, spin2: Spin, context: &MyEnvironment) -> FloatType {
        match (spin1, spin2) {
            (Spin::Some(c1), Spin::Some(c2)) => {
                if c1 == c2 {
                    0.
                } else {
                    let adh_id1 = context.env.cells[c1].cell.adh_id;
                    let adh_id2 = context.env.cells[c2].cell.adh_id;
                    2. * if adh_id1 == adh_id2 {
                        self.adh_energy
                    } else {
                        self.non_adh_energy
                    }
                }
            },
            (Spin::Some(_), Spin::Medium) | (Spin::Medium, Spin::Some(_)) => self.medium_energy,
            (Spin::Some(_), Spin::Solid) | (Spin::Solid, Spin::Some(_)) => self.solid_energy,
            _ => 0.
        }
    }
}