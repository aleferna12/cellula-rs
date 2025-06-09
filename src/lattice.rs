use std::ops::Index;
use rand::random_range;
use crate::pos::Pos2D;

pub struct Lattice<T> {
    pub width: usize,
    pub height: usize,
    array: Box<[T]>,
}

impl<T: Default + Clone> Lattice<T> {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            array: vec![T::default(); width * height].into_boxed_slice()
        }
    }

    // TODO: this should take an rng
    pub fn random_pos(&self) -> Pos2D<usize> {
        Pos2D::new(random_range(0..self.width - 1), random_range(0..self.height - 1))
    }
}

impl<T> Index<Pos2D<usize>> for Lattice<T> {
    type Output = T;

    fn index(&self, pos: Pos2D<usize>) -> &Self::Output {
        &self.array[pos.row_major(self.height)]
    }
}

impl<T> Index<(usize, usize)> for Lattice<T> {
    type Output = T;

    fn index(&self, pos: (usize, usize)) -> &Self::Output {
        &self[Pos2D::<usize>::from(pos)]
    }
}