//! Contains logic associated with [`Lattice`].

use crate::positional::boundaries::Boundary;
use crate::positional::neighborhood::Neighborhood;
use crate::positional::pos::Pos;
use crate::positional::rect::Rect;
use rand::RngExt;
use std::collections::VecDeque;
use std::ops::{Index, IndexMut};

/// A 2D rectangular lattice containing some objects of type `T`.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Lattice<T> {
    array: Box<[T]>,
    /// The lattice dimensions.
    pub rect: Rect<usize>
}

// Since the lattice size is naturally usize, boundary coord should be isize to avoid overflow errors
// Although technically it only has to be slightly larger than its defined size
impl<T> Lattice<T> {
    /// Makes a lattice from a pre-existing buffer of values stored as an array.
    pub fn from_array<const N: usize>(buf: [T; N], width: usize, height: usize) -> Option<Self> {
        if width * height != N {
            return None;
        }
        Some(Self {
            array: Box::new(buf),
            rect: Rect::new(Pos::new(0, 0), Pos::new(width, height))
        })
    }
    
    /// Makes a lattice from a pre-existing buffer of values stored in a slice.
    pub fn from_slice(slice: &[T], width: usize, height: usize) -> Option<Self>
    where 
        T: Clone {
        if width * height != slice.len() {
            return None;
        }
        Some(Self {
            array: Box::from(slice),
            rect: Rect::new(Pos::new(0, 0), Pos::new(width, height))
        })
    }

    /// Returns the width of the lattice.
    pub fn width(&self) -> usize {
        self.rect.width()
    }

    /// Returns the height of the lattice.
    pub fn height(&self) -> usize {
        self.rect.height()
    }

    /// Returns a random position from the lattice using a uniform distribution.
    pub fn random_pos(&self, rng: &mut impl RngExt) -> Pos<usize> {
        Pos::new(
            rng.random_range(0..self.width()),
            rng.random_range(0..self.height())
        )
    }

    /// Iterates over all lattice positions in column-major order.
    pub fn iter_positions(&self) -> impl Iterator<Item = Pos<usize>> + use<T> {
        self.rect.iter_positions()
    }

    /// Iterates over all values in the lattice in column-major order.
    pub fn iter_values(&self) -> impl Iterator<Item = &T> { self.array.iter() }

    /// Mutable version of [`Lattice::iter_values()`].
    pub fn iter_values_mut(&mut self) -> impl Iterator<Item = &mut T> { self.array.iter_mut() }

    /// Returns a slice to the values contained in the lattice.
    pub fn as_slice(&self) -> &[T] {
        &self.array
    }
}

impl<T: Default + Clone> Lattice<T> {
    /// Makes a lattice using `T`.
    pub fn new(width: usize, height: usize) -> Self {
        Self::from_slice(&vec![T::default(); width * height], width, height).unwrap()
    }

    /// Clears the lattice by setting all values to the default of the inner type of the lattice.
    pub fn clear(&mut self) {
        self.array.fill(T::default());
    }
}

impl<T: PartialEq> Lattice<T> {
    /// Searches for `value` by creating a box around `center_pos` and iterating all the positions inside it.
    pub fn search_box(
        &self,
        value: &T,
        center_pos: Pos<usize>,
        box_side: usize,
        bound: &impl Boundary<Coord = isize>
    ) -> impl Iterator<Item = Pos<usize>> {
        let center_isize = center_pos.cast_as::<isize>();
        let radius = (box_side / 2) as isize;
        let rect = Rect::new(
            (center_isize.x - radius, center_isize.y - radius).into(),
            (center_isize.x + radius, center_isize.y + radius).into(),
        );
        bound
            .valid_positions(rect.iter_positions())
            .filter_map(|pos| {
                let lat_pos = pos.cast_as();
                if self[lat_pos].eq(value) {
                    return Some(lat_pos);
                }
                None
            })
    }

    /// Searches for `value` using a BFS algorithm that iterates neighbors.
    ///
    /// Is considerably slower than [Lattice::search_box()].
    pub fn search_contiguous(
        &self,
        value: &T,
        start_pos: Pos<usize>,
        bound: &impl Boundary<Coord = isize>,
        neighborhood: &impl Neighborhood
    ) -> Box<[Pos<usize>]> {
        let mut found = vec![];
        let mut queue = VecDeque::from([start_pos.cast_as()]);
        let mut visited = Lattice::<bool>::from(self.rect.clone());
        visited[start_pos] = true;

        while let Some(pos) = queue.pop_front() {
            let lat_pos = pos.cast_as();
            if !self[lat_pos].eq(value) {
                continue;
            }
            bound
                .valid_positions(neighborhood.neighbors(pos))
                .for_each(|neigh| {
                    let lat_neigh = neigh.cast_as();
                    if !visited[lat_neigh] {
                        visited[lat_neigh] = true;
                        queue.push_back(neigh);
                    }
                });
            found.push(lat_pos);
        }
        found.into()
    }

    /// Returns the outline of a contiguous area containing `value`.
    ///
    /// The first outline position is automatically determined using [Lattice::search_box()] at `center_pos`.
    pub fn search_outline(
        &self,
        value: &T,
        center_pos: Pos<usize>,
        box_side: usize,
        bound: &impl Boundary<Coord = isize>,
        neighborhood: &impl Neighborhood
    ) -> Box<[Pos<usize>]> {
        let mut found = vec![];
        let border_pos = match self.search_box(
            value,
            center_pos,
            box_side,
            bound
        ).find_map(|pos| {
            if let Some(neigh) = bound.valid_pos(Pos::new(
                pos.x as isize - 1,
                pos.y as isize
            ))
                && &self[neigh.cast_as()] != value {
                Some(neigh)
            } else {
                None
            }
        }) {
            Some(neigh) => neigh,
            None => return found.into()
        };

        let mut queue = VecDeque::from([border_pos]);
        let mut visited = Lattice::<bool>::from(self.rect.clone());
        visited[border_pos.cast_as()] = true;

        while let Some(pos) = queue.pop_front() {
            let mut diff_spin_neighs = Vec::with_capacity(neighborhood.n_neighs().into());
            let mut has_value_neighbor = false;
            for neigh in bound.valid_positions(neighborhood.neighbors(pos)) {
                let neigh_pos = neigh.cast_as();
                if visited[neigh_pos] {
                    continue;
                }

                let neigh_spin = &self[neigh_pos];
                if neigh_spin == value {
                    has_value_neighbor = true;
                } else {
                    diff_spin_neighs.push(neigh);
                }
            }

            if has_value_neighbor {
                found.push(pos.cast_as());
                for neigh in diff_spin_neighs {
                    visited[neigh.cast_as()] = true;
                    queue.push_back(neigh)
                }
            }
        }
        found.into()
    }
}

impl<T> Index<Pos<usize>> for Lattice<T> {
    type Output = T;

    // Tested tiled-row-major and z order and normal row/column-major was fastest
    fn index(&self, pos: Pos<usize>) -> &Self::Output {
        &self.array[pos.col_major(self.height())]
    }
}

impl<T> IndexMut<Pos<usize>> for Lattice<T> {
    fn index_mut(&mut self, pos: Pos<usize>) -> &mut Self::Output {
        &mut self.array[pos.col_major(self.height())]
    }
}

impl<T: Default + Clone> From<Rect<usize>> for Lattice<T> {
    fn from(rect: Rect<usize>) -> Self {
        Self {
            array: vec![T::default(); rect.width() * rect.height()].into_boxed_slice(),
            rect,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::positional::boundaries::{ToLatticeBoundary, FastPeriodicBoundary};
    use crate::positional::neighborhood::MooreNeighborhood;
    use crate::positional::pos::{CastCoords, Pos};
    use crate::positional::rect::Rect;
    use rand::{rngs::StdRng, SeedableRng};

    #[test]
    fn test_lattice_indexing_get_and_set() {
        let rect = Rect::new((0, 0).into(), (3, 3).into());
        let mut lattice: Lattice<i32> = Lattice::from(rect);
        let pos = Pos::new(1, 2);
        lattice[pos] = 42;
        assert_eq!(lattice[pos], 42);
    }

    #[test]
    fn test_random_pos_within_bounds() {
        let rect = Rect::new((0, 0).into(), (10, 10).into());
        let lattice: Lattice<u8> = Lattice::from(rect);
        let mut rng = StdRng::seed_from_u64(42);
        for _ in 0..100 {
            let p = lattice.random_pos(&mut rng);
            assert!(p.x < lattice.width());
            assert!(p.y < lattice.height());
        }
    }

    #[test]
    fn test_search_box() {
        let rect = Rect::new((0., 0.).into(), (10., 10.).into());
        let mut lattice: Lattice<u8> = Lattice::from(rect.cast_coords());
        lattice[(5, 5).into()] = 1;
        lattice[(5, 6).into()] = 1;
        lattice[(4, 5).into()] = 1;
        let outline = lattice.search_box(
            &1,
            Pos::new(5, 5),
            5,
            &FastPeriodicBoundary::new(rect).to_lattice_boundary()
        ).collect::<Vec<_>>();
        assert_eq!(outline.len(), 3);
    }

    #[test]
    fn test_search_outline() {
        let rect = Rect::new((0., 0.).into(), (10., 10.).into());
        let mut lattice: Lattice<u8> = Lattice::from(rect.cast_coords());
        lattice[(5, 5).into()] = 1;
        lattice[(5, 6).into()] = 1;
        lattice[(4, 5).into()] = 1;
        let outline = lattice.search_outline(
            &1,
            Pos::new(5, 5),
            5,
            &FastPeriodicBoundary::new(rect).to_lattice_boundary(),
            &MooreNeighborhood::new(1)
        );
        assert_eq!(outline.len(), 12);
    }
}
