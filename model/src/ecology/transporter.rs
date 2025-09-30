use crate::pond::Pond;
use cellulars_lib::basic_cell::Alive;
use cellulars_lib::constants::CellIndex;
use cellulars_lib::environment::Habitable;
use cellulars_lib::lattice_entity::LatticeEntity;

pub trait Transporter {
    fn transport(&mut self, from: &mut Pond, to: &mut Pond, spins: Vec<CellIndex>);
}

pub struct WipeOut;

impl Transporter for WipeOut {
    fn transport(&mut self, from: &mut Pond, to: &mut Pond, spins: Vec<CellIndex>) {
        to.wipe_out();
        for spin in spins {
            let cell = from.env
                .cells
                .get_entity(spin)
                .expect_cell("tried to transport non-cell");
            let spin_to = to.env.cells.add(cell.birth(), None).index;
            for pos in from.env.search_cell_box(cell, from.cell_search_scaler) {
                from.env.grant_position(
                    pos,
                    LatticeEntity::Medium.discriminant(),
                );
                to.env.grant_position(
                    pos,
                    spin_to
                );
            }
        }
    }
}

