use crate::positional::boundary::Boundary;
use crate::positional::neighbourhood::Neighbourhood;
use crate::positional::pos::Pos;
use crate::positional::rect::Rect;
use rand::Rng;
use std::collections::VecDeque;
use std::ops::{Index, IndexMut};

pub struct Lattice<T> {
    array: Box<[T]>,
    pub rect: Rect<usize>
}

// Since the lattice size is naturally usize, boundary coord should be isize to avoid overflow errors
// Although technically it only has to be slightly larger than its defined size
impl<T> Lattice<T> {
    pub fn new(rect: Rect<usize>) -> Self
    where T: Default + Clone {
        Self {
            array: vec![T::default(); rect.width() * rect.height()]
                .into_boxed_slice(),
            rect,
        }
    }

    pub fn width(&self) -> usize {
        self.rect.width()
    }

    pub fn height(&self) -> usize {
        self.rect.height()
    }

    pub fn random_pos(&self, rng: &mut impl Rng) -> Pos<usize> {
        Pos::new(
            rng.random_range(0..self.width()),
            rng.random_range(0..self.height())
        )
    }

    pub fn iter_positions(&self) -> impl Iterator<Item = Pos<usize>> {
        self.rect.iter_positions()
    }

    pub fn iter_values(&self) -> impl Iterator<Item = &T> {
        self.iter_positions()
            .map(|pos| {
                &self[pos]
            })
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
        let center_isize = Pos::new(
            center_pos.x as isize,
            center_pos.y as isize
        );
        let radius = (box_side / 2) as isize;
        let rect = Rect::new(
            (center_isize.x - radius, center_isize.y - radius).into(),
            (center_isize.x + radius, center_isize.y + radius).into(),
        );
        bound
            .valid_positions(rect.iter_positions())
            .filter_map(|pos| {
                let lat_pos = pos.to_usize();
                if self[lat_pos].eq(value) {
                    return Some(lat_pos);
                }
                None
            })
    }

    /// Searches for `value` using a BFS algorithm that iterates neighbours.
    ///
    /// Is considerably slower than `search_box()`.
    pub fn search_contiguous(
        &self,
        value: &T,
        start_pos: Pos<usize>,
        bound: &impl Boundary<Coord = isize>,
        neighbourhood: &impl Neighbourhood
    ) -> Vec<Pos<usize>> {
        let mut found = vec![];
        let mut queue = VecDeque::from([start_pos.to_isize()]);
        let mut visited = Lattice::<bool>::new(self.rect.clone());
        visited[start_pos] = true;

        while let Some(pos) = queue.pop_front() {
            let lat_pos = pos.to_usize();
            if !self[lat_pos].eq(value) {
                continue;
            }
            bound
                .valid_positions(neighbourhood.neighbours(pos))
                .for_each(|neigh| {
                    let lat_neigh = neigh.to_usize();
                    if !visited[lat_neigh] {
                        visited[lat_neigh] = true;
                        queue.push_back(neigh);
                    }
                });
            found.push(lat_pos);
        }
        found
    }

    /// Returns the outline of a contiguous area containing `value`.
    ///
    /// The first outline position is automatically determined using `search_box()` at `center_pos`.
    pub fn search_outline(
        &self,
        value: &T,
        center_pos: Pos<usize>,
        box_side: usize,
        bound: &impl Boundary<Coord = isize>,
        neighbourhood: &impl Neighbourhood
    ) -> Vec<Pos<usize>> {
        let mut found = vec![];
        let border_pos = match self.search_box(
            value,
            center_pos,
            box_side,
            bound
        ).find_map(|pos| {
            if let Some(neigh) = bound.valid_pos(Pos::new(pos.x as isize - 1, pos.y as isize))
                && &self[neigh.to_usize()] != value {
                Some(neigh)
            } else {
                None
            }
        }) {
            Some(neigh) => neigh,
            None => return found
        };

        let mut queue = VecDeque::from([border_pos]);
        let mut visited = Lattice::<bool>::new(self.rect.clone());
        visited[border_pos.to_usize()] = true;

        while let Some(pos) = queue.pop_front() {
            let mut diff_spin_neighs = Vec::with_capacity(neighbourhood.n_neighs() as usize);
            let mut has_value_neighbour = false;
            for neigh in bound.valid_positions(neighbourhood.neighbours(pos)) {
                let neigh_pos = neigh.to_usize();
                if visited[neigh_pos] {
                    continue;
                }

                let neigh_spin = &self[neigh_pos];
                if neigh_spin == value {
                    has_value_neighbour = true;
                } else {
                    diff_spin_neighs.push(neigh);
                }
            }

            if has_value_neighbour {
                found.push(pos.to_usize());
                for neigh in diff_spin_neighs {
                    visited[neigh.to_usize()] = true;
                    queue.push_back(neigh)
                }
            }
        }
        found
    }
}

impl<T> Index<Pos<usize>> for Lattice<T> {
    type Output = T;

    fn index(&self, pos: Pos<usize>) -> &Self::Output {
        &self.array[pos.row_major(self.height())]
    }
}

impl<T> IndexMut<Pos<usize>> for Lattice<T> {
    fn index_mut(&mut self, pos: Pos<usize>) -> &mut Self::Output {
        &mut self.array[pos.row_major(self.height())]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::{BoundaryType, NeighbourhoodType};
    use crate::positional::boundary::AsLatticeBoundary;
    use crate::positional::pos::Pos;
    use crate::positional::rect::Rect;
    use rand::{rngs::StdRng, SeedableRng};

    #[test]
    fn test_lattice_indexing_get_and_set() {
        let rect = Rect::new((0, 0).into(), (3, 3).into());
        let mut lattice: Lattice<i32> = Lattice::new(rect);
        let pos = Pos::new(1, 2);
        lattice[pos] = 42;
        assert_eq!(lattice[pos], 42);
    }

    #[test]
    fn test_random_pos_within_bounds() {
        let rect = Rect::new((0, 0).into(), (10, 10).into());
        let lattice: Lattice<u8> = Lattice::new(rect);
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
        let mut lattice: Lattice<u8> = Lattice::new(rect.clone().try_into().unwrap());
        lattice[(5, 5).into()] = 1;
        lattice[(5, 6).into()] = 1;
        lattice[(4, 5).into()] = 1;
        let outline = lattice.search_box(
            &1,
            Pos::new(5, 5),
            5,
            &BoundaryType::new(rect).as_lattice_boundary().unwrap()
        ).collect::<Vec<_>>();
        assert_eq!(outline.len(), 3);
    }

    #[test]
    fn test_search_outline() {
        let rect = Rect::new((0., 0.).into(), (10., 10.).into());
        let mut lattice: Lattice<u8> = Lattice::new(rect.clone().try_into().unwrap());
        lattice[(5, 5).into()] = 1;
        lattice[(5, 6).into()] = 1;
        lattice[(4, 5).into()] = 1;
        let outline = lattice.search_outline(
            &1,
            Pos::new(5, 5),
            5,
            &BoundaryType::new(rect).as_lattice_boundary().unwrap(),
            &NeighbourhoodType::new(1)
        );
        assert_eq!(outline.len(), 12);
    }
}
