use std::ops::{Index, IndexMut};
use rand::Rng;
use crate::boundary::LatticeBoundary;
use crate::pos::Pos2D;

pub struct Lattice<T, B: LatticeBoundary> {
    array: Box<[T]>,
    pub bound: B
}

impl<T: Default + Copy, B: LatticeBoundary> Lattice<T, B> {
    pub fn new(bound: B) -> Self {
        Self {
            array: vec![T::default(); bound.rect().width() * bound.rect().height()].into_boxed_slice(),
            bound,
        }
    }
}
impl<T, B: LatticeBoundary> Lattice<T, B> {
    pub fn width(&self) -> usize {
        self.bound.rect().width()
    }

    pub fn height(&self) -> usize {
        self.bound.rect().height()
    }

    pub fn random_pos(&self, rng: &mut impl Rng) -> Pos2D<usize> {
        Pos2D::new(rng.random_range(0..self.width() - 1), rng.random_range(0..self.height() - 1))
    }

    pub fn iter_positions(&self) -> impl Iterator<Item = Pos2D<usize>> {
        self.bound.rect().iter_positions()
    }
}

impl<T: Copy, B: LatticeBoundary> Lattice<T, B> {
    pub fn iter_values(&self) -> impl Iterator<Item = T> {
        self.bound
            .rect()
            .iter_positions()
            .map(|pos| {
                self[pos]
            })
    }
}

impl<T, B: LatticeBoundary> Index<Pos2D<usize>> for Lattice<T, B> {
    type Output = T;

    fn index(&self, pos: Pos2D<usize>) -> &Self::Output {
        &self.array[pos.row_major(self.height())]
    }
}
impl<T, B: LatticeBoundary> Index<(usize, usize)> for Lattice<T, B> {
    type Output = T;

    fn index(&self, pos: (usize, usize)) -> &Self::Output {
        &self[Pos2D::<usize>::from(pos)]
    }
}

impl<T, B: LatticeBoundary> IndexMut<Pos2D<usize>> for Lattice<T, B> {
    fn index_mut(&mut self, pos: Pos2D<usize>) -> &mut Self::Output {
        &mut self.array[pos.row_major(self.height())]
    }
}
impl<T, B: LatticeBoundary> IndexMut<(usize, usize)> for Lattice<T, B> {
    fn index_mut(&mut self, pos: (usize, usize)) -> &mut Self::Output {
        &mut self[Pos2D::<usize>::from(pos)]
    }
}