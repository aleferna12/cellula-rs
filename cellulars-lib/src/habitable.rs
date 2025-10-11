use crate::basic_cell::{Cellular, RelCell};
use crate::environment::{EdgesUpdate, Environment};
use crate::positional::boundaries::ToLatticeBoundary;
use crate::positional::neighbourhood::Neighbourhood;
use crate::positional::pos::Pos;
use crate::spin::Spin;

pub trait Habitable {
    type Cell: Cellular;

    fn env(&self) -> &Environment<Self::Cell, impl Neighbourhood, impl ToLatticeBoundary>;

    fn env_mut(&mut self) -> &mut Environment<Self::Cell, impl Neighbourhood, impl ToLatticeBoundary>;

    fn grant_position(&mut self, pos: Pos<usize>, to: Spin) -> EdgesUpdate;

    fn spawn_cell(
        &mut self,
        empty_cell: Self::Cell,
        positions: impl IntoIterator<Item = Pos<usize>>
    ) -> &RelCell<Self::Cell> {
        let cell_index = self.env_mut().cells.add(empty_cell).index;
        let new_spin = Spin::Some(cell_index);
        for pos in positions {
            self.grant_position(pos, new_spin);
        }
        self.env().cells.get_cell(cell_index)
    }

    fn spawn_solid(&mut self, positions: impl Iterator<Item = Pos<usize>>) {
        for pos in positions {
            self.grant_position(pos, Spin::Solid);
        }
    }
}