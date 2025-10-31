use crate::my_environment::MyEnvironment;
use cellulars_lib::adhesion::AdhesionSystem;
use cellulars_lib::spin::Spin;
use crate::cell::Cell;

#[derive(Clone)]
pub struct SpeciesAdhesion {
    pub amoeba_energy: f32,
    pub bacteria_energy: f32,
    pub interspecies_energy: f32,
    pub medium_energy: f32,
    pub solid_energy: f32,
}

impl SpeciesAdhesion {
    fn cell_cell_energy(&self, cell1: &Cell, cell2: &Cell) -> f32 {
        match (cell1, cell2) {
            (Cell::Amoeba(_), Cell::Amoeba(_)) => self.amoeba_energy,
            (Cell::Bacterium(_), Cell::Bacterium(_)) => self.bacteria_energy,
            _ => self.interspecies_energy,
        }
    }
}

impl AdhesionSystem for SpeciesAdhesion {
    type Context = MyEnvironment;

    fn adhesion_energy(&self, spin1: Spin, spin2: Spin, env: &Self::Context) -> f32 {
        match (spin1, spin2) {
            (Spin::Some(c1), Spin::Some(c2)) => {
                if c1 == c2 {
                    0.
                } else {
                    let cell1 = &env.cells.get_cell(c1).cell;
                    let cell2 = &env.cells.get_cell(c2).cell;
                    2. * self.cell_cell_energy(cell1, cell2)
                }
            }
            (Spin::Some(_), Spin::Medium) | (Spin::Medium, Spin::Some(_)) => self.medium_energy,
            (Spin::Some(_), Spin::Solid) | (Spin::Solid, Spin::Some(_)) => self.solid_energy,
            _ => 0.
        }
    }
}