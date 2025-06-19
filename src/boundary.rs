use std::ops::{Mul, Sub};
use crate::pos::{Pos2D, Rect};

pub trait Boundary<T> {
    /// Expose the boundary as a `Rect`.
    fn rect(&self) -> &Rect<T>;
    
    fn inbounds(&self, pos: Pos2D<T>) -> bool;

    /// Validates that positions are in bounds.
    ///
    /// With fixed boundary conditions, that means filtering invalid positions.
    fn validate_positions<S>(&self, pos_it: S) -> impl Iterator<Item = Pos2D<T>>
    where
        S: Iterator<Item = Pos2D<T>>;
}

pub struct FixedBoundary<T> {
    rect: Rect<T>
}
impl<T> FixedBoundary<T>
where
    T: Sub<Output = T>
        + Mul<Output = T>
        + PartialOrd
        + Copy {
    pub fn new(rect: Rect<T>) -> Self {
        Self { rect }
    }
}

impl<T: PartialOrd + Copy> Boundary<T> for FixedBoundary<T> {
    fn rect(&self) -> &Rect<T> {
        &self.rect
    }

    fn inbounds(&self, pos: Pos2D<T>) -> bool {
        (self.rect.min.x..self.rect.max.x).contains(&pos.x) && (self.rect.min.y..self.rect.max.y).contains(&pos.y)
    }

    fn validate_positions<S: Iterator<Item = Pos2D<T>>>(
        &self,
        pos_it: S
    ) -> impl Iterator<Item = Pos2D<T>> {
        pos_it.filter(|pos| { 
            self.inbounds(*pos) 
        })
    }
}