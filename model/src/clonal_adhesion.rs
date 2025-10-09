use cellulars_lib::adhesion::{AdhesionSystem, StaticAdhesion};
use cellulars_lib::entity::{Entity, Spin};
use cellulars_lib::symmetric_table::SymmetricTable;

#[derive(Clone)]
pub struct ClonalAdhesion {
    pub clone_energy: f32,
    pub static_adhesion: StaticAdhesion
}

impl ClonalAdhesion {
    pub fn new(clone_energy: f32, static_adhesion: StaticAdhesion) -> Self {
        Self {
            clone_energy,
            static_adhesion
        }
    }
}

impl AdhesionSystem for ClonalAdhesion {
    type Context = SymmetricTable<bool>;

    fn adhesion_energy(&self, entity1: Spin, entity2: Spin, clones_table: &Self::Context) -> f32 {
        if let (Entity::Some(c1), Entity::Some(c2)) = (entity1, entity2) {
            if c1 == c2 {
                return 0.
            }
            if clones_table[(c1 as usize, c2 as usize)] {
                return 2. * self.clone_energy;
            }
        }
        // Handle all other cases
        self.static_adhesion.adhesion_energy(entity1, entity2, &())
    }
}

#[cfg(test)]
mod tests {
    use crate::clonal_adhesion::ClonalAdhesion;
    use cellulars_lib::adhesion::{AdhesionSystem, StaticAdhesion};
    use cellulars_lib::entity::{Entity, Spin};
    use cellulars_lib::symmetric_table::SymmetricTable;

    fn make_static_adhesion() -> StaticAdhesion {
        StaticAdhesion {
            cell_energy: 2.,
            medium_energy: 3.,
            solid_energy: 4.,
        }
    }

    #[test]
    fn test_clonal_adhesion() {
        let clonal_adhesion = ClonalAdhesion::new(
            1.,
            make_static_adhesion()
        );

        let mut clones = SymmetricTable::new(2);
        // Initially clone_pairs empty
        assert_eq!(
            clonal_adhesion.adhesion_energy(Spin::Some(0), Spin::Some(0), &clones),
            0.
        );
        assert_eq!(
            clonal_adhesion.adhesion_energy(Spin::Some(0), Entity::Some(1), &clones),
            2. * clonal_adhesion.static_adhesion.cell_energy
        );

        clones[(0, 1)] = true;
        assert_eq!(
            clonal_adhesion.adhesion_energy(Spin::Some(0), Spin::Some(1), &clones),
            2. * clonal_adhesion.clone_energy
        );
        // ClonalAdhesion falls back to StaticAdhesion for Medium and Solid
        assert_eq!(
            clonal_adhesion.adhesion_energy(Spin::Some(0), Spin::Medium, &clones),
            clonal_adhesion.static_adhesion.medium_energy
        );
        assert_eq!(
            clonal_adhesion.adhesion_energy(Spin::Solid, Spin::Some(0), &clones),
            clonal_adhesion.static_adhesion.solid_energy
        );
    }
}