use std::ops::{Index, IndexMut};
use rand::Rng;
use crate::pos::{Pos2D, Rect};

pub struct Lattice<T> {
    rect: Rect<usize>,
    array: Box<[T]>,
}

impl<T: Default + Clone> Lattice<T> {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            rect: Rect::new((0, 0).into(), (width, height).into()),
            array: vec![T::default(); width * height].into_boxed_slice()
        }
    }
}
impl<T> Lattice<T> {
    pub fn width(&self) -> usize {
        self.rect.width()
    }

    pub fn height(&self) -> usize {
        self.rect.height()
    }

    pub fn random_pos(&self, rng: &mut impl Rng) -> Pos2D<usize> {
        Pos2D::new(rng.random_range(0..self.width() - 1), rng.random_range(0..self.height() - 1))
    }

    pub fn filter_inbounds<S: Iterator<Item = Pos2D<usize>>>(
        &self,
        pos_it: S
    ) -> impl Iterator<Item = Pos2D<usize>> + use<S, T> {
    let rect = self.rect;
    pos_it.filter(move |pos| { rect.inbounds(pos) })
    }
}

impl<T> Index<Pos2D<usize>> for Lattice<T> {
    type Output = T;

    fn index(&self, pos: Pos2D<usize>) -> &Self::Output {
        &self.array[pos.row_major(self.height())]
    }
}
impl<T> Index<(usize, usize)> for Lattice<T> {
    type Output = T;

    fn index(&self, pos: (usize, usize)) -> &Self::Output {
        &self[Pos2D::<usize>::from(pos)]
    }
}

impl<T> IndexMut<Pos2D<usize>> for Lattice<T> {
    fn index_mut(&mut self, pos: Pos2D<usize>) -> &mut Self::Output {
        &mut self.array[pos.row_major(self.height())]
    }
}
impl<T> IndexMut<(usize, usize)> for Lattice<T> {
    fn index_mut(&mut self, pos: (usize, usize)) -> &mut Self::Output {
        &mut self[Pos2D::<usize>::from(pos)]
    }
}