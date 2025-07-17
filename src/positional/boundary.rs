use std::ops::Sub;
use crate::positional::pos::Pos;
use crate::positional::rect::Rect;
use num::traits::Euclid;
use num::Num;

pub trait Boundary {
    type Coord;

    fn rect(&self) -> &Rect<Self::Coord>;

    fn valid_pos(&self, pos: Pos<Self::Coord>) -> Option<Pos<Self::Coord>>;

    #[inline]
    fn valid_positions(
        &self,
        positions: impl Iterator<Item = Pos<Self::Coord>>,
    ) -> impl Iterator<Item = Pos<Self::Coord>> {
        positions.filter_map(|pos| self.valid_pos(pos))
    }

    fn displacement(
        &self,
        pos1: Pos<Self::Coord>,
        pos2: Pos<Self::Coord>,
    ) -> (Self::Coord, Self::Coord);
}

/// Fixed boundary with clamped positions.
#[derive(Clone)]
pub struct FixedBoundary<T> {
    rect: Rect<T>,
}

impl<T> FixedBoundary<T> {
    pub fn new(rect: Rect<T>) -> Self {
        Self { rect }
    }
}

impl<T> Boundary for FixedBoundary<T>
where
    T: PartialOrd + Copy + Sub<Output = T> {
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
            return None;
        }
        Some(pos)
    }

    #[inline]
    fn displacement(&self, pos1: Pos<T>, pos2: Pos<T>) -> (T, T) {
        (pos2.x - pos1.x, pos2.y - pos1.y)
    }
}

/// Safe periodic boundary: always wraps.
#[derive(Clone)]
pub struct SafePeriodicBoundary<T> {
    rect: Rect<T>,
}

impl<T> SafePeriodicBoundary<T>
where
    T: Copy + Num + Euclid + From<u8> {
    pub fn new(rect: Rect<T>) -> Self {
        Self { rect }
    }

    #[inline]
    fn wrap_scalar(val: T, min: T, max: T) -> T {
        let range = max - min;
        let offset = (val - min).rem_euclid(&range);
        min + offset
    }

    #[inline]
    fn periodic_displacement(&self, pos1: Pos<T>, pos2: Pos<T>) -> (T, T) {
        let two = T::from(2);
        let w = self.rect.width();
        let h = self.rect.height();
        let dx = ((pos2.x - pos1.x + w / two).rem_euclid(&w)) - w / two;
        let dy = ((pos2.y - pos1.y + h / two).rem_euclid(&h)) - h / two;
        (dx, dy)
    }
}

impl<T> Boundary for SafePeriodicBoundary<T>
where
    T: Copy + Num + Euclid + From<u8> {
    type Coord = T;

    fn rect(&self) -> &Rect<T> {
        &self.rect
    }

    #[inline]
    fn valid_pos(&self, pos: Pos<T>) -> Option<Pos<T>> {
        Some(Pos::new(
            Self::wrap_scalar(pos.x, self.rect.min.x, self.rect.max.x),
            Self::wrap_scalar(pos.y, self.rect.min.y, self.rect.max.y),
        ))
    }

    #[inline]
    fn displacement(&self, pos1: Pos<T>, pos2: Pos<T>) -> (T, T) {
        self.periodic_displacement(pos1, pos2)
    }
}

/// Unsafe version of periodic boundary — assumes positions are close.
#[derive(Clone)]
pub struct UnsafePeriodicBoundary<T> {
    rect: Rect<T>,
}

impl<T> UnsafePeriodicBoundary<T>
where
    T: Copy + Num + PartialOrd {
    pub fn new(rect: Rect<T>) -> Self {
        Self { rect }
    }

    #[inline]
    fn wrap_scalar(&self, val: T, min: T, max: T) -> T {
        if val < min {
            max - (min - val)
        } else if val >= max {
            min + (val - max)
        } else {
            val
        }
    }

    #[inline]
    fn periodic_displacement(&self, pos1: Pos<T>, pos2: Pos<T>) -> (T, T)
    where
        T: Euclid + From<u8>,
    {
        let two = T::from(2);
        let w = self.rect.width();
        let h = self.rect.height();
        let dx = ((pos2.x - pos1.x + w / two).rem_euclid(&w)) - w / two;
        let dy = ((pos2.y - pos1.y + h / two).rem_euclid(&h)) - h / two;
        (dx, dy)
    }
}

impl<T> Boundary for UnsafePeriodicBoundary<T>
where
    T: Copy + Num + Euclid + PartialOrd + From<u8> {
    type Coord = T;

    fn rect(&self) -> &Rect<T> {
        &self.rect
    }

    #[inline]
    fn valid_pos(&self, pos: Pos<T>) -> Option<Pos<T>> {
        Some(Pos::new(
            self.wrap_scalar(pos.x, self.rect.min.x, self.rect.max.x),
            self.wrap_scalar(pos.y, self.rect.min.y, self.rect.max.y),
        ))
    }

    #[inline]
    fn displacement(&self, pos1: Pos<T>, pos2: Pos<T>) -> (T, T) {
        self.periodic_displacement(pos1, pos2)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_periodic_boundary() {
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