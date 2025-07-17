use std::ops::Sub;
use crate::positional::pos::Pos;
use crate::positional::rect::Rect;
use num::traits::Euclid;
use num::Num;

pub trait Boundary {
    type Coord;
    
    /// Expose the boundary as a `Rect`.
    fn rect(&self) -> &Rect<Self::Coord>;
    
    /// Validates that positions are in bounds.
    ///
    /// With fixed boundary conditions, that means filtering invalid positions.
    fn valid_pos(&self, pos: Pos<Self::Coord>) -> Option<Pos<Self::Coord>>;

    #[inline]
    fn valid_positions(
        &self,
        positions: impl Iterator<Item = Pos<Self::Coord>>
    ) -> impl Iterator<Item = Pos<Self::Coord>> {
        positions.filter_map(|pos| self.valid_pos(pos))
    }

    fn displacement(&self, pos1: Pos<Self::Coord>, pos2: Pos<Self::Coord>) -> (Self::Coord, Self::Coord);
}

pub trait PeriodicBoundary: Boundary
where
    Self::Coord: From<u8> + Num + Copy + Euclid {
    #[inline]
    fn periodic_displacement(&self, pos1: Pos<Self::Coord>, pos2: Pos<Self::Coord>) -> (Self::Coord, Self::Coord) {
        let two = Self::Coord::from(2);
        let dx = ((pos2.x - pos1.x + self.rect().width() / two)
            .rem_euclid(&self.rect().width())) - self.rect().width() / two;
        let dy = ((pos2.y - pos1.y + self.rect().height() / two)
            .rem_euclid(&self.rect().height())) - self.rect().height() / two;
        (dx, dy)
    }

    #[inline]
    fn periodic_valid_pos(&self, pos: Pos<Self::Coord>) -> Option<Pos<Self::Coord>> {
        let p = Pos::new(
            Self::wrap_scalar(pos.x, self.rect().min.x, self.rect().max.x),
            Self::wrap_scalar(pos.y, self.rect().min.y, self.rect().max.y)
        );
        Some(p)
    }

    #[inline]
    fn wrap_scalar(val: Self::Coord, min: Self::Coord, max: Self::Coord) -> Self::Coord {
        let range = max - min;
        let mut offset = val - min;
        offset = offset.rem_euclid(&range);
        min + offset
    }
}

#[derive(Clone)]
pub struct FixedBoundary<T> {
    rect: Rect<T>
}

impl<T> FixedBoundary<T> {
    pub fn new(rect: Rect<T>) -> Self {
        Self { rect }
    }
}

impl<T> Boundary for FixedBoundary<T>
where 
    T: PartialOrd 
        + Copy 
        + Sub<Output = T> {
    type Coord = T;

    fn rect(&self) -> &Rect<T> {
        &self.rect
    }

    #[inline]
    fn valid_pos(&self, pos: Pos<T>) -> Option<Pos<T>> {
        if !(self.rect.min.x..self.rect.max.x).contains(&pos.x) {
            return None;
        }
        if !(self.rect.min.y..self.rect.max.y).contains(&pos.y) {
            return None
        }
        Some(pos)
    }

    #[inline]
    fn displacement(&self, pos1: Pos<Self::Coord>, pos2: Pos<Self::Coord>) -> (Self::Coord, Self::Coord) {
        (pos2.x - pos1.x, pos2.y - pos1.y)
    }
}

#[derive(Clone)]
pub struct SafePeriodicBoundary<T> {
    rect: Rect<T>
}

impl<T> SafePeriodicBoundary<T>
where
    T: Copy + Num + Euclid {
    pub fn new(rect: Rect<T>) -> Self {
        Self { rect }
    }
}

impl<T> Boundary for SafePeriodicBoundary<T>
where
    T: Num + Copy + Euclid + From<u8> {
        type Coord = T;

        fn rect(&self) -> &Rect<T> {
            &self.rect
        }

    #[inline]
    fn valid_pos(&self, pos: Pos<Self::Coord>) -> Option<Pos<Self::Coord>> {
        self.periodic_valid_pos(pos)
    }

    #[inline]
    fn displacement(&self, pos1: Pos<Self::Coord>, pos2: Pos<Self::Coord>) -> (Self::Coord, Self::Coord) {
        self.periodic_displacement(pos1, pos2)
    }
}

impl<T> PeriodicBoundary for SafePeriodicBoundary<T>
where
    T: Num + Copy + Euclid + From<u8> {}

/// This struct can only validate positions that are at most one `width()` or `height()` away from the boundaries.
///
/// <div class="warning">
///
/// Only use when you are confident that all input positions are close to the boundary.
///
/// </div>
#[derive(Clone)]
pub struct UnsafePeriodicBoundary<T> {
    rect: Rect<T>
}

impl<T> UnsafePeriodicBoundary<T>
where
    T: Copy
        + Num 
        + PartialOrd {
    pub fn new(rect: Rect<T>) -> Self {
        Self { rect }
    }
}

impl<T> Boundary for UnsafePeriodicBoundary<T>
where
    T: Num + Copy + Euclid + From<u8> {
    type Coord = T;

    fn rect(&self) -> &Rect<T> {
        &self.rect
    }
    
    #[inline]
    fn valid_pos(&self, pos: Pos<Self::Coord>) -> Option<Pos<Self::Coord>> {
        self.periodic_valid_pos(pos)
    }

    #[inline]
    fn displacement(&self, pos1: Pos<Self::Coord>, pos2: Pos<Self::Coord>) -> (Self::Coord, Self::Coord) {
        self.periodic_displacement(pos1, pos2)
    }
}

impl<T> PeriodicBoundary for UnsafePeriodicBoundary<T>
where
    T: Num + Copy + Euclid + From<u8> {
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_safe_periodic_boundary() {
        let per = SafePeriodicBoundary::new(Rect::new((0, 0).into(), (10, 10).into()));
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