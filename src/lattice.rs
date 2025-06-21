use std::ops::{Index, IndexMut};
use rand::Rng;
use crate::boundary::Boundary;
use crate::pos::Pos2D;

pub struct Lattice<T, B> {
    array: Box<[T]>,
    pub bound: B
}
// Since the lattice size is naturally usize, boundary coord should be isize to avoid overflow errors
// Although technically it only has to be slightly larger than its defined size
impl<T: Default + Copy, B: Boundary<Coord = isize>> Lattice<T, B> {
    pub fn new(bound: B) -> Self {
        Self {
            array: vec![T::default(); bound.rect().width() as usize * bound.rect().height() as usize]
                .into_boxed_slice(),
            bound,
        }
    }

    pub fn width(&self) -> usize {
        self.bound.rect().width() as usize
    }

    pub fn height(&self) -> usize {
        self.bound.rect().height() as usize
    }

    pub fn random_pos(&self, rng: &mut impl Rng) -> Pos2D<usize> {
        Pos2D::new(
            rng.random_range(0..self.width()),
            rng.random_range(0..self.height())
        )
    }

    pub fn iter_positions(&self) -> impl Iterator<Item = Pos2D<usize>> {
        self.bound.rect().iter_positions().map(|p| Pos2D::new(
            p.x as usize,
            p.y as usize
        ))
    }

    pub fn iter_values(&self) -> impl Iterator<Item = T> {
        self.iter_positions()
            .map(|pos| {
                self[pos]
            })
    }
}

impl<T: Copy + Default, B: Boundary<Coord = isize>> Index<Pos2D<usize>> for Lattice<T, B> {
    type Output = T;

    fn index(&self, pos: Pos2D<usize>) -> &Self::Output {
        &self.array[pos.row_major(self.height())]
    }
}
impl<T: Copy + Default, B: Boundary<Coord = isize>> Index<(usize, usize)> for Lattice<T, B> {
    type Output = T;

    fn index(&self, pos: (usize, usize)) -> &Self::Output {
        &self[Pos2D::<usize>::from(pos)]
    }
}

impl<T: Copy + Default, B: Boundary<Coord = isize>> IndexMut<Pos2D<usize>> for Lattice<T, B> {
    fn index_mut(&mut self, pos: Pos2D<usize>) -> &mut Self::Output {
        &mut self.array[pos.row_major(self.height())]
    }
}
impl<T: Copy + Default, B: Boundary<Coord = isize>> IndexMut<(usize, usize)> for Lattice<T, B> {
    fn index_mut(&mut self, pos: (usize, usize)) -> &mut Self::Output {
        &mut self[Pos2D::<usize>::from(pos)]
    }
}