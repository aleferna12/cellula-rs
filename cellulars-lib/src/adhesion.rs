use crate::basic_cell::RelCell;
use crate::entity::Entity;

pub trait AdhesionSystem {
    fn adhesion_energy<C>(&self, entity1: Entity<&RelCell<C>>, entity2: Entity<&RelCell<C>>) -> f32;
}

#[derive(Clone)]
pub struct StaticAdhesion {
    pub cell_energy: f32,
    pub medium_energy: f32,
    pub solid_energy: f32
}

impl AdhesionSystem for StaticAdhesion {
    fn adhesion_energy<C>(
        &self,
        entity1: Entity<&RelCell<C>>,
        entity2: Entity<&RelCell<C>>
    ) -> f32 {
        match (entity1, entity2) {
            (Entity::Some(c1), Entity::Some(c2)) => {
                if c1.index == c2.index {
                    0.
                } else {
                    2. * self.cell_energy
                }
            }
            (Entity::Some(_), Entity::Medium) | (Entity::Medium, Entity::Some(_)) => self.medium_energy,
            (Entity::Some(_), Entity::Solid) | (Entity::Solid, Entity::Some(_)) => self.solid_energy,
            _ => 0.
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::basic_cell::BasicCell;
    use crate::constants::CellIndex;

    // Helper to create a RelCell with given cell index and mom by mocking and overriding
    fn make_rel_with_index(index: CellIndex) -> RelCell<BasicCell> {
        RelCell {
            index,
            cell: BasicCell::new_empty(10)
        }
    }
    
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

        let cell1 = make_rel_with_index(1);
        let cell2 = make_rel_with_index(2);

        assert_eq!(
            static_adhesion.adhesion_energy(Entity::Some(&cell1), Entity::Some(&cell1)),
            0.
        );
        assert_eq!(
            static_adhesion.adhesion_energy(Entity::Some(&cell1), Entity::Some(&cell2)),
            2. * static_adhesion.cell_energy
        );
        assert_eq!(
            static_adhesion.adhesion_energy(Entity::Some(&cell1), Entity::Medium),
            static_adhesion.medium_energy
        );
        assert_eq!(
            static_adhesion.adhesion_energy(Entity::Solid, Entity::Some(&cell1)),
            static_adhesion.solid_energy
        );
    }
}
