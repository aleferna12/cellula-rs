use crate::pond::Pond;
use cellulars_lib::cellular::Cellular;
use cellulars_lib::constants::Spin;
use cellulars_lib::lattice_entity::LatticeEntity;

pub trait Transporter {
    fn transport(&mut self, from: &mut Pond, to: &mut Pond, spins: Vec<Spin>);
}

pub struct WipeOut {
    cell_search_radius: f32
}

impl Transporter for WipeOut {
    fn transport(&mut self, from: &mut Pond, to: &mut Pond, spins: Vec<Spin>) {
        to.wipe_out();
        for spin in spins {
            let cell = from.env
                .cells
                .get_entity(spin)
                .expect_cell("tried to transport non-cell");
            let spin_to = to.env.cells.push(cell.birth(), None).spin;
            for pos in from.env.space.search_cell_box(cell, self.cell_search_radius) {
                from.ca.grant_position(
                    pos,
                    LatticeEntity::Medium.discriminant(),
                    &mut from.env
                );
                to.ca.grant_position(
                    pos,
                    spin_to,
                    &mut to.env
                );
            }
        }
    }
}

