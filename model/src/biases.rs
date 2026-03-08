use cellulars::prelude::*;
use crate::my_cell::CellType;
use crate::my_environment::MyEnvironment;

#[derive(Debug, Clone)]
pub struct Biases {
    pub chem_bias: ChemotaxisBias
}

impl CopyBias<MyEnvironment> for Biases {
    fn bias(&self, pos_source: Pos<usize>, pos_target: Pos<usize>, context: &MyEnvironment) -> FloatType {
        let spin = context.env.cell_lattice[pos_source];
        if let Spin::Some(cell_index) = spin
            && let CellType::Migrating =  context.env.cells[cell_index].cell.cell_type {
            self.chem_bias.bias(pos_source, pos_target, &context.chem_lattice)
        } else {
            0.
        }
    }
}