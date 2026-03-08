use cellulars::prelude::*;
use crate::my_cell::CellType;
use crate::my_environment::MyEnvironment;

#[derive(Clone)]
pub struct ChemotaxisBias {
    pub chemotaxis_mu: FloatType
}

impl CopyBias<MyEnvironment> for ChemotaxisBias {
    fn bias(&self, pos_source: Pos<usize>, pos_target: Pos<usize>, env: &MyEnvironment) -> FloatType {
        let Spin::Some(cell_index) = env.env.cell_lattice[pos_source] else {
            return 0.;
        };
        let rel_cell = &env.env.cells[cell_index];
        if let CellType::Dividing = rel_cell.cell.cell_type {
            return 0.;
        }

        let (dx1, dy1) = env.env.bounds.boundary.displacement(
            rel_cell.cell.center(),
            Pos::new(pos_target.x as FloatType, pos_target.y as FloatType)
        );
        let (dx2, dy2) = env.env.bounds.boundary.displacement(
            rel_cell.cell.center(),
            rel_cell.cell.chem_center()
        );

        let dot = dx1 * dx2 + dy1 * dy2;
        let norm1_sq = dx1 * dx1 + dy1 * dy1;
        let norm2_sq = dx2 * dx2 + dy2 * dy2;
        let denom = (norm1_sq * norm2_sq).sqrt();

        if denom <= 0. {
            0.
        } else {
            -self.chemotaxis_mu * (dot / denom)
        }
    }
}