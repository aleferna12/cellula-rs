use num::Num;
use num::traits::Euclid;
use crate::pos::{Pos2D, Rect};

pub trait Boundary {
    type Coord;
    
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
    T: Copy + Num + Euclid {
    fn wrap_scalar(val: T, min: T, max: T) -> T {
        let range = max - min;
        let mut offset = val - min;
        offset = offset.rem_euclid(&range);
        min + offset
    }
}

impl<T> Boundary for PeriodicBoundary<T>
where
    T: Copy + Num + Euclid {
    type Coord = T;

    fn rect(&self) -> &Rect<T> {
        &self.rect
    }
    
    fn valid_pos(&self, pos: Pos2D<Self::Coord>) -> Option<Pos2D<Self::Coord>> {
        let p = Pos2D::new(
            Self::wrap_scalar(pos.x, self.rect.min.x, self.rect.max.x),
            Self::wrap_scalar(pos.y, self.rect.min.y, self.rect.max.y)
        );
        Some(p)
    }
}

pub struct UnsafePeriodicBoundary<T> {
    bound: PeriodicBoundary<T>
}
impl<T> UnsafePeriodicBoundary<T>
where
    T: Copy
        + Num 
        + PartialOrd {
    pub fn new(rect: Rect<T>) -> Self {
        Self { bound: PeriodicBoundary::new(rect) }
    }
    
    pub fn wrap_scalar(&self, val: T, min: T, max: T) -> T {
        if val < min {
            max - (min - val)
        } else if val >= max {
            min + (val - max)
        } else {
            val
        }
    }
}

impl<T> Boundary for UnsafePeriodicBoundary<T>
where
    T: Copy + Num + Euclid + PartialOrd {
    type Coord = T;

    fn rect(&self) -> &Rect<Self::Coord> {
        self.bound.rect()
    }

    /// This wraps the position inside the boundary ONCE.
    ///
    /// If the position is more than `width()` or `height` away, this will not produce a valid position.
    /// If you need this reassurance, use `PeriodicBoundary`, which is slower.
    fn valid_pos(&self, pos: Pos2D<Self::Coord>) -> Option<Pos2D<Self::Coord>> {
        Some(Pos2D::new(
            self.wrap_scalar(pos.x, self.rect().min.x, self.rect().max.x),
            self.wrap_scalar(pos.y, self.rect().min.y, self.rect().max.y)
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_periodic_boundary() {
        let per = PeriodicBoundary::new(Rect::new((0, 0).into(), (10, 10).into()));
        assert_eq!(per.valid_pos((1, 1).into()).unwrap(), (1, 1).into());
        assert_eq!(per.valid_pos((-1, -1).into()).unwrap(), (9, 9).into());
        assert_eq!(per.valid_pos((10, 10).into()).unwrap(), (0, 0).into());
        assert_eq!(per.valid_pos((11, 11).into()).unwrap(), (1, 1).into());
        assert_eq!(per.valid_pos((30, 30).into()).unwrap(), (0, 0).into())
    }
    
    #[test]
    fn test_unsafe_periodic_boundary() {
        let unsafeper = UnsafePeriodicBoundary::new(Rect::new((0, 0).into(), (10, 10).into()));
        assert_eq!(unsafeper.valid_pos((1, 1).into()).unwrap(), (1, 1).into());
        assert_eq!(unsafeper.valid_pos((-1, -1).into()).unwrap(), (9, 9).into());
        assert_eq!(unsafeper.valid_pos((10, 10).into()).unwrap(), (0, 0).into());
        assert_eq!(unsafeper.valid_pos((11, 11).into()).unwrap(), (1, 1).into());
        assert_ne!(unsafeper.valid_pos((30, 30).into()).unwrap(), (0, 0).into())
    }
}