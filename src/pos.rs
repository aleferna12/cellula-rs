use std::hash::{Hash, Hasher};
use crate::pos::EdgeError::{NotNeighbours, SamePosition};

pub enum EdgeError {
    SamePosition,
    NotNeighbours
}

/// 2D position in space.
#[derive(PartialEq, Eq)]
pub struct Pos2D<T> {
    pub x: T,
    pub y: T
}

impl<T> Pos2D<T> {
    pub fn new(x: T, y: T) -> Self {
        Self { x, y }
    }
}

// This is only used for indexing and therefore can be implemented for usize only
impl Pos2D<usize> {
    pub fn row_major(&self, height: usize) -> usize {
        self.x * height + self.y
    }
}

impl<T> From<(T, T)> for Pos2D<T> {
    fn from(value: (T, T)) -> Self {
        Pos2D::<T>::new(value.0, value.1)
    }
}

#[derive(Eq)]
pub struct Edge {
    p1: Pos2D<usize>,
    p2: Pos2D<usize>,
    direction: usize,
    pub neigh_r: u32
}

impl Edge {
    pub fn new(&self, p1: Pos2D<usize>, p2: Pos2D<usize>, neigh_r: u32) -> Result<Self, EdgeError> {
        let dir = (p1.x - p2.x + 1) * 3 + (p1.y - p2.y) + 1;
        match dir {
            0..4 => Ok(Self { p1, p2, direction: dir, neigh_r }),
            4 => Err(SamePosition),
            5..9 => Ok(Self { p1: p2, p2: p1, direction: 8 - dir, neigh_r }),
            _ => Err(NotNeighbours)
        }
    }

    pub fn p1(&self) -> &Pos2D<usize> {
        &self.p1
    }

    pub fn p2(&self) -> &Pos2D<usize> {
        &self.p2
    }
}

impl PartialEq<Self> for Edge {
    fn eq(&self, other: &Self) -> bool {
        (self.p1() == other.p1()) & (self.direction == other.direction)
    }
}

// TODO: test the perfect hash algorithm that steven sent me
impl Hash for Edge {
    fn hash<H: Hasher>(&self, state: &mut H) {
        // TODO: test if saving neigh_r as the matrix height itself makes things faster
        (self.p1().row_major(self.neigh_r as usize * 2 + 1) + self.direction).hash(state)
    }
}