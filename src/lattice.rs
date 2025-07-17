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
        self.rect.iter_positions().map(|p| Pos::new(
            p.x,
            p.y
        ))
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
