//! Contains logic associated with adhesion systems.

use crate::spin::Spin;

/// Trait defining a way to calculate adhesion at the interface between two spins.
pub trait AdhesionSystem {
    /// Context required to calculate the adhesion energy.
    type Context;
    /// Returns the energy at the interface between `spin1` and `spin2`, given a context.
    fn adhesion_energy(&self, spin1: Spin, spin2: Spin, context: &Self::Context) -> f32;
}

/// An adhesion system that only considers [Spin]s to determine adhesion energies.
#[derive(Clone)]
pub struct StaticAdhesion {
    /// Energy at a cell-cell interface.
    pub cell_energy: f32,
    /// Energy at a cell-medium interface.
    pub medium_energy: f32,
    /// Energy at a cell-solid interface.
    pub solid_energy: f32
}

impl AdhesionSystem for StaticAdhesion {
    type Context = ();
    fn adhesion_energy(
        &self,
        spin1: Spin,
        spin2: Spin,
        _: &Self::Context,
    ) -> f32 {
        match (spin1, spin2) {
            (Spin::Some(c1), Spin::Some(c2)) => {
                if c1 == c2 {
                    0.
                } else {
                    2. * self.cell_energy
                }
            }
            (Spin::Some(_), Spin::Medium) | (Spin::Medium, Spin::Some(_)) => self.medium_energy,
            (Spin::Some(_), Spin::Solid) | (Spin::Solid, Spin::Some(_)) => self.solid_energy,
            _ => 0.
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_static_adhesion() -> StaticAdhesion {
        StaticAdhesion {
            cell_energy: 3.,
            medium_energy: 1.5,
            solid_energy: 2.,
        }
    }

    #[test]
    fn test_static_adhesion() {
        let static_adhesion = make_static_adhesion();

        assert_eq!(
            static_adhesion.adhesion_energy(Spin::Some(1), Spin::Some(1), &()),
            0.
        );
        assert_eq!(
            static_adhesion.adhesion_energy(Spin::Some(1), Spin::Some(2), &()),
            2. * static_adhesion.cell_energy
        );
        assert_eq!(
            static_adhesion.adhesion_energy(Spin::Some(1), Spin::Medium, &()),
            static_adhesion.medium_energy
        );
        assert_eq!(
            static_adhesion.adhesion_energy(Spin::Solid, Spin::Some(1), &()),
            static_adhesion.solid_energy
        );
    }
}
