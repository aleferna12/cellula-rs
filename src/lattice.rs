use std::ops::{Index, IndexMut};
use rand::Rng;
use crate::boundary::Boundary;
use crate::pos::{GeneralCoord, LatticeCoord, Pos2D};

pub struct Lattice<T, B: Boundary<Coord = GeneralCoord>> {
    array: Box<[T]>,
    pub bound: B
}

impl<T: Default + Copy, B: Boundary<Coord = GeneralCoord>> Lattice<T, B> {
    pub fn new(bound: B) -> Self {
        Self {
            array: vec![T::default(); bound.rect().width() as usize * bound.rect().height() as usize]
                .into_boxed_slice(),
            bound,
        }
    }

    pub fn width(&self) -> LatticeCoord {
        self.bound.rect().width() as LatticeCoord
    }

    pub fn height(&self) -> LatticeCoord {
        self.bound.rect().height() as LatticeCoord
    }

    pub fn random_pos(&self, rng: &mut impl Rng) -> Pos2D<LatticeCoord> {
        Pos2D::new(
            rng.random_range(0..self.width() - 1) as LatticeCoord,
            rng.random_range(0..self.height() - 1) as LatticeCoord
        )
    }

    pub fn iter_positions(&self) -> impl Iterator<Item = Pos2D<LatticeCoord>> {
        self.bound.rect().iter_positions().map(|p| Pos2D::new(
            p.x as LatticeCoord,
            p.y as LatticeCoord
        ))
    }

    pub fn iter_values(&self) -> impl Iterator<Item = T> {
        self.iter_positions()
            .map(|pos| {
                self[pos]
            })
    }
}

impl<T: Default + Copy, B: Boundary<Coord = GeneralCoord>> Index<Pos2D<LatticeCoord>> for Lattice<T, B> {
    type Output = T;

    fn index(&self, pos: Pos2D<LatticeCoord>) -> &Self::Output {
        &self.array[pos.row_major(self.height())]
    }
}
impl<T: Default + Copy, B: Boundary<Coord = GeneralCoord>> Index<(LatticeCoord, LatticeCoord)> for Lattice<T, B> {
    type Output = T;

    fn index(&self, pos: (LatticeCoord, LatticeCoord)) -> &Self::Output {
        &self[Pos2D::<LatticeCoord>::from(pos)]
    }
}

impl<T: Default + Copy, B: Boundary<Coord = GeneralCoord>> IndexMut<Pos2D<LatticeCoord>> for Lattice<T, B> {
    fn index_mut(&mut self, pos: Pos2D<LatticeCoord>) -> &mut Self::Output {
        &mut self.array[pos.row_major(self.height())]
    }
}
impl<T: Default + Copy, B: Boundary<Coord = GeneralCoord>> IndexMut<(LatticeCoord, LatticeCoord)> for Lattice<T, B> {
    fn index_mut(&mut self, pos: (LatticeCoord, LatticeCoord)) -> &mut Self::Output {
        &mut self[Pos2D::<LatticeCoord>::from(pos)]
    }
}