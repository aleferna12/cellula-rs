use std::hash::{Hash, Hasher};
use std::mem;
use crate::pos::EdgeError::{NotNeighbours, SamePosition};

#[derive(Debug)]
pub enum EdgeError {
    SamePosition,
    NotNeighbours
}

/// 2D position in space.
#[derive(PartialEq, Eq, Copy, Clone, Debug)]
pub struct Pos2D<T> {
    pub x: T,
    pub y: T
}

impl<T> Pos2D<T> {
    pub fn new(x: T, y: T) -> Self {
        Self { x, y }
    }
}

impl Pos2D<usize> {
    fn pack_u32(&self) -> u32 {
        ((self.x as u32) << 16) | self.y as u32
    }

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
#[derive(Eq)]
pub struct Edge {
    p1: Pos2D<usize>,
    p2: Pos2D<usize>
}

impl Edge {
    pub fn new(p1: Pos2D<usize>, p2: Pos2D<usize>, neigh_r: u8) -> Result<Self, EdgeError> {
        let cx = p1.x.abs_diff(p2.x);
        let cy = p1.y.abs_diff(p2.y);
        let sum = cx + cy;
        if sum == 0 {
            return Err(SamePosition);
        }
        if sum > (neigh_r * 2) as usize {
            return Err(NotNeighbours);
        }
        Ok(Self { p1, p2})
    }

    #[inline(always)]
    fn hash_u64(&self) -> u64 {
        let mut u1 = self.p1.pack_u32();
        let mut u2 = self.p2.pack_u32();
        if u1 > u2 {
            mem::swap(&mut u1, &mut u2);
        }
        ((u1 as u64) << 32) | (u2 as u64)
    }
}

impl PartialEq for Edge {
    fn eq(&self, other: &Self) -> bool {
        self.hash_u64() == other.hash_u64()
    }
}

// TODO: test the perfect hash algorithm that steven sent me (nothing else helped)
impl Hash for Edge {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.hash_u64().hash(state);
    }
}

pub struct Rect<T> (pub Pos2D<T>, pub Pos2D<T>);
impl<T> Rect<T> {
    pub fn new(p1: Pos2D<T>, p2: Pos2D<T>) -> Self {
        Self(p1, p2)
    }
}

impl Rect<usize> {
    pub fn iterate_pos(&self) -> RectAreaIt {
        RectAreaIt::new(self)
    }
}

pub struct RectAreaIt<'a> {
    curr: Pos2D<usize>,
    rect: &'a Rect<usize>
}
impl<'a> RectAreaIt<'a> {
    fn new(rect: &'a Rect<usize>) -> Self {
        let mut p = rect.0.clone();
        p.x -= 1;
        Self {
            curr: p,
            rect
        }
    }
}

impl Iterator for RectAreaIt<'_> {
    type Item = Pos2D<usize>;

    fn next(&mut self) -> Option<Self::Item> {
        let same_row = self.curr.x < self.rect.1.x - 1;
        if same_row {
            self.curr.x += 1;
            return Some(self.curr)
        }
        if self.curr.y < self.rect.1.y - 1 {
            self.curr.x = self.rect.0.x;
            self.curr.y += 1;
            return Some(self.curr);
        }
        None
    }
}