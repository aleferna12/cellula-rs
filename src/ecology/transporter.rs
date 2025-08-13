use crate::cell::Cellular;
use crate::constants::Spin;
use crate::environment::LatticeEntity;
use crate::pond::Pond;

pub trait Transporter {
    fn transport(&mut self, from: &mut Pond, to: &mut Pond, spins: Vec<Spin>);
}

pub struct WipeOut;

impl Transporter for WipeOut {
    fn transport(&mut self, from: &mut Pond, to: &mut Pond, spins: Vec<Spin>) {
        to.wipe_out();
        for spin in spins {
            let cell = from.env
                .cells
                .get_entity(spin)
                .expect_cell("tried to transport non-cell");
            let spin_from = cell.spin;
            let spin_to = to.env.cells.push(cell.birth(), None).spin;

            for pos in from.env.space.search_cell_box(&cell, from.env.cell_search_radius) {
                from.ca.shift_position(
                    pos,
                    spin_from,
                    LatticeEntity::Medium.discriminant(),
                    &mut from.env
                );
                to.ca.shift_position(
                    pos,
                    LatticeEntity::Medium.discriminant(),
                    spin_to,
                    &mut to.env
                );
            }
        }
    }
}

