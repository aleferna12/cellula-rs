use std::ops::{Add, Sub};
use crate::pos::{Pos2D, Rect};

pub trait Boundary {
    type Coord: Copy;
    
    /// Expose the boundary as a `Rect`.
    fn rect(&self) -> &Rect<Self::Coord>;
    
    /// Validates that positions are in bounds.
    ///
    /// With fixed boundary conditions, that means filtering invalid positions.
    fn valid_pos(&self, pos: Pos2D<Self::Coord>) -> Option<Pos2D<Self::Coord>>;
}

pub struct FixedBoundary<T> {
    rect: Rect<T>
}
impl<T> FixedBoundary<T> {
    pub fn new(rect: Rect<T>) -> Self {
        Self { rect }
    }
}

impl<T: PartialOrd + Copy> Boundary for FixedBoundary<T> {
    type Coord = T;

    fn rect(&self) -> &Rect<T> {
        &self.rect
    }

    fn valid_pos(&self, pos: Pos2D<T>) -> Option<Pos2D<T>> {
        if !(self.rect.min.x..self.rect.max.x).contains(&pos.x) {
            return None;
        } 
        if !(self.rect.min.y..self.rect.max.y).contains(&pos.y) {
            return None
        }
        Some(pos)
    }
}

pub struct PeriodicBoundary<T> {
    rect: Rect<T>
}
impl<T> PeriodicBoundary<T> {
    pub fn new(rect: Rect<T>) -> Self {
        Self { rect }
    }
}
impl<T> PeriodicBoundary<T>
where
    T: Copy 
        + PartialOrd 
        + Add<Output = T> 
        + Sub<Output = T> {
    fn wrap_scalar(val: T, min: T, max: T) -> T {
        if val < min {
            max - (min - val)
        } else if val >= max { 
            min + (val - max)
        } else { 
            val
        }
    }
}

impl<T> Boundary for PeriodicBoundary<T>
where
    T: Copy
        + PartialOrd
        + Add<Output = T>
        + Sub<Output = T> {
    type Coord = T;

    fn rect(&self) -> &Rect<T> {
        &self.rect
    }

    // TODO: this is quite slow compared to FixedBoundary, try to implement wrap_scalar as if else statements 
    fn valid_pos(&self, pos: Pos2D<Self::Coord>) -> Option<Pos2D<Self::Coord>> {
        let p = Pos2D::new(
            Self::wrap_scalar(pos.x, self.rect.min.x, self.rect.max.x),
            Self::wrap_scalar(pos.y, self.rect.min.y, self.rect.max.y)
        );
        Some(p)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_periodic() {
        let per = PeriodicBoundary::new(Rect::new((0, 0).into(), (10, 10).into()));
        assert_eq!(per.valid_pos((1, 1).into()).unwrap(), (1, 1).into());
        assert_eq!(per.valid_pos((-1, -1).into()).unwrap(), (9, 9).into());
        assert_eq!(per.valid_pos((10, 10).into()).unwrap(), (0, 0).into());
        assert_eq!(per.valid_pos((11, 11).into()).unwrap(), (1, 1).into());
    }
}