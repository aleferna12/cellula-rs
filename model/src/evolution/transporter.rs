use crate::pond::Pond;
use cellulars_lib::basic_cell::Alive;
use cellulars_lib::constants::CellIndex;
use cellulars_lib::habitable::Habitable;
use cellulars_lib::spin::Spin;

pub trait Transporter {
    fn transport(&mut self, from: &mut Pond, to: &mut Pond, cell_indexes: Vec<CellIndex>);
}

pub struct WipeOut;

impl Transporter for WipeOut {
    fn transport(&mut self, from: &mut Pond, to: &mut Pond, cell_indexes: Vec<CellIndex>) {
        to.wipe_out();
        for cell_index in cell_indexes {
            let cell = from.env.cells.get_cell(cell_index);
            let index_to = to.env.cells.add(cell.birth()).index;
            for pos in from.env.search_cell_box(cell, from.env.cell_search_scaler) {
                from.env.grant_position(
                    pos,
                    Spin::Medium,
                );
                to.env.grant_position(
                    pos,
                    Spin::Some(index_to)
                );
            }
        }
    }
}

