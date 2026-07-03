use std::fmt::{Debug, Formatter};
use crate::constants::{CellIndex, FloatType};
use crate::empty_cell::Empty;
use crate::environment::EdgesUpdate;
use crate::prelude::{Alive, AsEnv, Cellular, Environment, FastPeriodicBoundary, HasCenter, MooreNeighborhood, Neighborhood, Pos, Spin, SymmetricTable, ToLatticeBoundary, TransferPosition};
use crate::traits::track_neighbors::TrackNeighbors;

pub struct NeighEnvironment<C, N = MooreNeighborhood, B: ToLatticeBoundary = FastPeriodicBoundary<FloatType>> {
    pub env: Environment<C, N, B>,
    neigh_tracker: NeighTracker
}

impl<C, N, B: ToLatticeBoundary> NeighEnvironment<C, N, B> {
    pub fn new(env: Environment<C, N, B>, max_cells: CellIndex) -> NeighEnvironment<C, N, B> {
        let array =  vec![0; max_cells as usize].into_boxed_slice();
        NeighEnvironment {
            neigh_tracker: NeighTracker {
                neigh_table: SymmetricTable::new(max_cells as usize),
                med_array: array.clone(),
                solid_array: array,
            },
            env
        }
    }

    pub fn max_cells(&self) -> CellIndex {
        self.neigh_tracker.neigh_table.length() as u32
    }
}

impl<C, N: Neighborhood, B: ToLatticeBoundary> NeighEnvironment<C, N, B> {
    pub fn update_neighbors(&mut self, pos: Pos<usize>, to: Spin) {
        let spin = self.env.cell_lattice[pos];
        let valid_neighs = self.env.valid_neighbors(pos);
        for neigh in valid_neighs {
            let neigh_spin = self.env.cell_lattice[neigh];
            if  neigh_spin != to {
                self.neigh_tracker.shift_neighbor_contact(to, neigh_spin, true);
            } 
            if neigh_spin != spin {
                self.neigh_tracker.shift_neighbor_contact(spin, neigh_spin, false);
            }
        }
    }
}

impl<C, N, B> TransferPosition for NeighEnvironment<C, N, B>
where
    C: Cellular + Alive + HasCenter + Empty,
    N: Neighborhood,
    B: ToLatticeBoundary<Coord = FloatType> {
    fn transfer_position(&mut self, pos: Pos<usize>, to: Spin) -> EdgesUpdate {
        self.update_neighbors(pos, to);
        self.env.transfer_position(pos, to)
    }
}

impl<C, N: Neighborhood, B: ToLatticeBoundary<Coord = FloatType>> AsEnv for NeighEnvironment<C, N, B> {
    type Cell = C;
    type Coord = FloatType;

    fn env(&self) -> &Environment<Self::Cell, impl Neighborhood, impl ToLatticeBoundary<Coord = Self::Coord>> {
        &self.env
    }

    fn env_mut(&mut self) -> &mut Environment<Self::Cell, impl Neighborhood, impl ToLatticeBoundary<Coord = Self::Coord>> {
        &mut self.env
    }
}

impl<C, N, B: ToLatticeBoundary> TrackNeighbors for NeighEnvironment<C, N, B> {
    fn neighbor_contacts(&self, spin1: Spin, spin2: Spin) -> Option<u32> {
        match (spin1, spin2) {
            (Spin::Some(ci1), Spin::Some(ci2)) => {
                self.neigh_tracker.neigh_table.get((ci1 as usize, ci2 as usize))
            },
            (Spin::Some(ci), Spin::Medium) | (Spin::Medium, Spin::Some(ci)) => {
                self.neigh_tracker.med_array.get(ci as usize)
            },
            (Spin::Some(ci), Spin::Solid) | (Spin::Solid, Spin::Some(ci)) => {
                self.neigh_tracker.solid_array.get(ci as usize)
            },
            _ => None
        }.copied()
    }
}

impl<C, N, B: ToLatticeBoundary> Clone for NeighEnvironment<C, N, B>
where
    C: Clone,
    N: Clone,
    B: Clone,
    B::LatticeBoundary: Clone {
    fn clone(&self) -> Self {
        Self {
            env: self.env.clone(),
            neigh_tracker: self.neigh_tracker.clone()
        }
    }
}

impl<C, N, B: ToLatticeBoundary> Debug for NeighEnvironment<C, N, B>
where
    C: Debug,
    N: Debug,
    B: Debug,
    B::LatticeBoundary: Debug {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("NeighEnvironment")
            .field("env", &self.env)
            .field("neigh_tracker", &self.neigh_tracker)
            .finish()
    }
}

impl<C, N, B:ToLatticeBoundary> PartialEq for NeighEnvironment<C, N, B>
where
    C: PartialEq,
    N: PartialEq,
    B: PartialEq,
    B::LatticeBoundary: PartialEq {
    fn eq(&self, other: &Self) -> bool {
        // Assume correctness: no need to check neigh_tracker
        self.env == other.env
    }
}

impl<C, N, B:ToLatticeBoundary> Eq for NeighEnvironment<C, N, B>
where C: Eq, N: Eq, B: Eq, B::LatticeBoundary: Eq {}

#[derive(Clone, Debug)]
struct NeighTracker {
    neigh_table: SymmetricTable<u32>,
    med_array: Box<[u32]>,
    solid_array: Box<[u32]>
}

impl NeighTracker {
    fn shift_neighbor_contact(&mut self, spin1: Spin, spin2: Spin, adding: bool) -> bool {
        let diff = if adding { 1 } else { -1 };
        let to_modify = match (spin1, spin2) {
            (Spin::Some(ci1), Spin::Some(ci2)) => {
                match self.neigh_table.get_mut((ci1 as usize, ci2 as usize)) {
                    Some(val) => val,
                    _ => return false
                }
            },
            (Spin::Some(ci), Spin::Medium) | (Spin::Medium, Spin::Some(ci)) => {
                match self.med_array.get_mut(ci as usize) {
                    Some(val) => val,
                    _ => return false
                }
            },
            (Spin::Some(ci), Spin::Solid) | (Spin::Solid, Spin::Some(ci)) => {
                match self.solid_array.get_mut(ci as usize) {
                    Some(val) => val,
                    _ => return false
                }
            },
            _ => return false
        };
        *to_modify = to_modify.saturating_add_signed(diff);
        true
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::{Boundaries, Cell, Rect};
    use super::*;
    
    #[test]
    fn add_neigh() {
        let mut env = NeighEnvironment::<Cell, _, _>::new(Environment::new_empty(
            MooreNeighborhood::new(1),
            Boundaries::new(FastPeriodicBoundary::new(Rect::new(Pos::new(0., 0.), Pos::new(10., 10.))))
        ), 10);
        env.env.cells.push(Cell::new_empty(1));
        env.env.cells.push(Cell::new_empty(1));
        env.transfer_position(Pos::new(0, 0), Spin::Some(0));
        env.transfer_position(Pos::new(0, 1), Spin::Some(1));
        assert_eq!(env.neighbor_contacts(Spin::Some(0), Spin::Some(1)), Some(1));
    }
}