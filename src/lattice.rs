use std::ops::{Index, IndexMut};
use rand::Rng;
use crate::pos::{Pos2D, Rect};

#[derive(Debug)]
pub enum LatticeEntity<C> {
    Medium,
    SomeCell(C)
}
impl<C: std::fmt::Debug> LatticeEntity<C> {
    pub fn unwrap(self) -> C {
        match self {
            LatticeEntity::SomeCell(cell) => cell,
            _ => panic!("called `LatticeEntity::unwrap()` on a `{:?}` value", self)
        }
    }
}

pub struct Lattice<T> {
    rect: Rect<usize>,
    array: Box<[T]>,
}

impl<T: Default + Copy> Lattice<T> {
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
    
    /// Validates that positions are in bounds.
    /// 
    /// With fixed boundary conditions, that means filtering invalid positions.
    pub fn validate<S: Iterator<Item = Pos2D<usize>>>(
        &self,
        pos_it: S
    ) -> impl Iterator<Item = Pos2D<usize>> + use<S, T> {
        let rect = self.rect.clone();
        pos_it.filter(move |pos| { rect.inbounds(pos.clone()) })
    }

    pub fn iter_positions(&self) -> impl Iterator<Item = Pos2D<usize>> {
        self.rect.iter_positions()
    }
}

impl<T: Copy> Lattice<T> {
    pub fn iter_values(&self) -> impl Iterator<Item = T> {
        self.rect
            .iter_positions()
            .map(|pos| {
                self[pos]
            })
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