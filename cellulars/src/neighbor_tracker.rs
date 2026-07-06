use std::fmt::Debug;
use crate::prelude::{CellIndex, Environment, Neighborhood, Pos, Spin, SymmetricTable, ToLatticeBoundary};
use crate::traits::track_neighbors::TrackNeighbors;

/// Tracks neighboring cells in an [`Environment`].
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct NeighborTracker {
    neigh_table: SymmetricTable<u32>,
    med_array: Box<[u32]>,
    solid_array: Box<[u32]>
}

impl NeighborTracker {
    /// Makes a new [`NeighborTracker`].
    ///
    /// Neighbor contacts are only tracked between cells with lower index than `max_cells`.
    /// The memory footprint of this struct grows with `max_cells^2`.
    pub fn new(max_cells: CellIndex) -> Self {
        let neigh_array = vec![0; max_cells as usize].into_boxed_slice();
        Self {
            neigh_table: SymmetricTable::new(max_cells as usize),
            med_array: neigh_array.clone(),
            solid_array: neigh_array
        }
    }

    /// Maximum number of cells tracked.
    pub fn max_cells(&self) -> CellIndex {
        self.neigh_table.length() as CellIndex
    }

    /// Update the neighborhoods of cells in `env` given that the spin `to` is going to be copied into `pos`.
    ///
    /// If the algorithm encounters any cells with a spin >= [`NeighborTracker::max_cells()`], returns `false`.
    /// Otherwise, returns `true`. Regardless, the neighborhoods of cells with spin smaller than max
    ///
    /// <div class="warning">
    /// This must be called before the cell lattice is updated.
    /// </div>
    pub fn update_neighbors<C, N: Neighborhood, B: ToLatticeBoundary>(
        &mut self,
        env: &Environment<C, N, B>,
        pos: Pos<usize>,
        to: Spin
    ) -> bool {
        let spin = env.cell_lattice[pos];
        let valid_neighs = env.valid_neighbors(pos);
        let mut success = true;
        for neigh in valid_neighs {
            let neigh_spin = env.cell_lattice[neigh];
            if  neigh_spin != to && !self.shift_neighbor_contact(to, neigh_spin, true) {
                success = false;
            }
            if neigh_spin != spin && !self.shift_neighbor_contact(spin, neigh_spin, false) {
                success = false;
            }
        }
        success
    }


    pub fn initialize_from_env<C, N: Neighborhood, B: ToLatticeBoundary>(
        &mut self,
        env: &Environment<C, N, B>
    ) -> bool {
        let mut success = true;
        for pos in env.cell_lattice.iter_positions() {
            let spin = env.cell_lattice[pos];
            for neigh_pos in env.valid_neighbors(pos) {
                let neigh_spin = env.cell_lattice[neigh_pos];
                if spin != neigh_spin && !self.shift_neighbor_contact(spin, neigh_spin, true) {
                    success = false;
                }
            }
        }
        // Every position is counted twice, so we divide by two
        for index in self.neigh_table.iter_index_pairs(None, None).collect::<Box<_>>() {
            self.neigh_table[index] /= 2;
        }
        for index in 0..self.max_cells() {
            self.med_array[index as usize] /= 2;
            self.solid_array[index as usize] /= 2;
        }
        success
    }

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

impl TrackNeighbors for NeighborTracker {
    fn neighbor_contacts(&self, spin1: Spin, spin2: Spin) -> Option<u32> {
        match (spin1, spin2) {
            (Spin::Some(ci1), Spin::Some(ci2)) => {
                self.neigh_table.get((ci1 as usize, ci2 as usize))
            },
            (Spin::Some(ci), Spin::Medium) | (Spin::Medium, Spin::Some(ci)) => {
                self.med_array.get(ci as usize)
            },
            (Spin::Some(ci), Spin::Solid) | (Spin::Solid, Spin::Some(ci)) => {
                self.solid_array.get(ci as usize)
            },
            _ => None
        }.copied()
    }
}

#[cfg(test)]
mod tests {
    use crate::prelude::{Boundaries, Cell, FastPeriodicBoundary, MooreNeighborhood, Rect, TransferPosition};
    use super::*;
    
    #[test]
    fn add_neigh() {
        let mut env = Environment::<Cell, _, _>::new_empty(
            MooreNeighborhood::new(1),
            Boundaries::new(FastPeriodicBoundary::new(Rect::new(Pos::new(0., 0.), Pos::new(10., 10.))))
        );
        let mut nt = NeighborTracker::new(2);
        env.cells.push(Cell::new_empty(1));
        env.cells.push(Cell::new_empty(1));
        nt.update_neighbors(&env, Pos::new(0, 0), Spin::Some(0));
        env.transfer_position(Pos::new(0, 0), Spin::Some(0));
        nt.update_neighbors(&env, Pos::new(0, 1), Spin::Some(1));
        env.transfer_position(Pos::new(0, 1), Spin::Some(1));
        assert_eq!(nt.neighbor_contacts(Spin::Some(0), Spin::Some(1)), Some(1));
    }
}