use std::hash::{Hash, Hasher};
use std::mem;
use crate::pos::EdgeError::{NotNeighbours, SamePosition};

#[derive(Debug)]
pub enum EdgeError {
    SamePosition,
    NotNeighbours
}

/// 2D position in space.
#[derive(PartialEq, Eq, Copy, Clone)]
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
    #[inline]
    pub fn row_major(&self, height: usize) -> usize {
        self.x * height + self.y
    }
}

impl<T> From<(T, T)> for Pos2D<T> {
    fn from(value: (T, T)) -> Self {
        Pos2D::<T>::new(value.0, value.1)
    }
}

// This currently only supports a Moore neighbourhood of 1
#[derive(Eq, PartialEq)]
pub struct Edge {
    p1: Pos2D<usize>,
    p2: Pos2D<usize>,
    direction: usize,
}

impl Edge {
    pub fn new(p1: Pos2D<usize>, p2: Pos2D<usize>) -> Result<Self, EdgeError> {
        let dir = ((p1.x as i32 - p2.x as i32 + 1) * 3 + (p1.y as i32 - p2.y as i32) + 1) as usize;
        match dir {
            4 => Err(SamePosition),
            9..usize::MAX => Err(NotNeighbours),
            _ => Ok(Self { p1, p2, direction: dir }),
        }
    }
}

// TODO: test the perfect hash algorithm that steven sent me (mine was not faster)
impl Hash for Edge {
    fn hash<H: Hasher>(&self, state: &mut H) {
        let mut p1 = self.p1;
        let mut p2 = self.p2;
        if self.direction > 3 {
            mem::swap(&mut p1, &mut p2);
        }
        // This is worse than just default implementation
        // let index = p.row_major(20) + dir;
        // index.hash(state);
        (p1.x, p1.y, p2.x, p2.y).hash(state)
    }
}