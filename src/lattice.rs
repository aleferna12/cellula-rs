use crate::positional::pos::Pos;
use crate::positional::rect::Rect;
use rand::Rng;
use std::ops::{Index, IndexMut};

pub struct Lattice<T> {
    array: Box<[T]>,
    pub rect: Rect<usize>
}

// Since the lattice size is naturally usize, boundary coord should be isize to avoid overflow errors
// Although technically it only has to be slightly larger than its defined size
impl<T: Default + Copy> Lattice<T> {
    pub fn new(rect: Rect<usize>) -> Self {
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

    pub fn iter_values(&self) -> impl Iterator<Item = T> {
        self.iter_positions()
            .map(|pos| {
                self[pos]
            })
    }
}

impl<T: Copy + Default> Index<Pos<usize>> for Lattice<T> {
    type Output = T;

    fn index(&self, pos: Pos<usize>) -> &Self::Output {
        &self.array[pos.row_major(self.height())]
    }
}

impl<T: Copy + Default> IndexMut<Pos<usize>> for Lattice<T> {
    fn index_mut(&mut self, pos: Pos<usize>) -> &mut Self::Output {
        &mut self.array[pos.row_major(self.height())]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
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
}
