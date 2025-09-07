// TODO! remove inlines in this mod and benchmark
use crate::positional::pos::Pos;
use crate::positional::rect::{Rect, RectConversionError};
use num::traits::Euclid;
use num::Num;
use std::error::Error;
use std::ops::Sub;

pub struct Boundaries<B: ToLatticeBoundary> {
    pub boundary: B,
    pub lattice_boundary: B::LatticeBoundary,
}

impl<B: ToLatticeBoundary> Boundaries<B> {
    pub fn new(bound: B) -> Result<Self, B::Error>
    where
        B: ToLatticeBoundary<Coord = f32>,
        B::Error: 'static + Error {
        Ok(Self {
            lattice_boundary: bound.to_lattice_boundary()?,
            boundary: bound,
        })
    }
}

pub trait Boundary {
    type Coord;

    /// Expose the boundary as a `Rect`.
    fn rect(&self) -> &Rect<Self::Coord>;

    /// Validates that positions are in bounds.
    ///
    /// With fixed boundary conditions, that means filtering invalid positions.
    fn valid_pos(&self, pos: Pos<Self::Coord>) -> Option<Pos<Self::Coord>>;

    #[inline(always)]
    fn valid_positions(
        &self,
        positions: impl Iterator<Item = Pos<Self::Coord>>
    ) -> impl Iterator<Item = Pos<Self::Coord>> {
        positions.filter_map(|pos| self.valid_pos(pos))
    }

    fn displacement(&self, pos1: Pos<Self::Coord>, pos2: Pos<Self::Coord>) -> (Self::Coord, Self::Coord);
}

/// Private trait that provides default implementations for periodic boundary types.
trait PeriodicBoundary: Boundary
where
    Self::Coord: From<u8> + Num + Copy + Euclid {
    #[inline(always)]
    fn periodic_displacement(&self, pos1: Pos<Self::Coord>, pos2: Pos<Self::Coord>) -> (Self::Coord, Self::Coord) {
        let two = Self::Coord::from(2);
        let dx = ((pos2.x - pos1.x + self.rect().width() / two)
            .rem_euclid(&self.rect().width())) - self.rect().width() / two;
        let dy = ((pos2.y - pos1.y + self.rect().height() / two)
            .rem_euclid(&self.rect().height())) - self.rect().height() / two;
        (dx, dy)
    }

    #[inline(always)]
    fn periodic_valid_pos(&self, pos: Pos<Self::Coord>) -> Option<Pos<Self::Coord>> {
        let p = Pos::new(
            Self::wrap_scalar(pos.x, self.rect().min.x, self.rect().max.x),
            Self::wrap_scalar(pos.y, self.rect().min.y, self.rect().max.y)
        );
        Some(p)
    }

    fn wrap_scalar(val: Self::Coord, min: Self::Coord, max: Self::Coord) -> Self::Coord;
}

pub trait ToLatticeBoundary: Boundary {
    type LatticeBoundary: Boundary<Coord = isize>;
    type Error;
    fn to_lattice_boundary(&self) -> Result<Self::LatticeBoundary, Self::Error>;
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

    #[inline(always)]
    fn rect(&self) -> &Rect<T> {
        &self.rect
    }

    #[inline(always)]
    fn valid_pos(&self, pos: Pos<T>) -> Option<Pos<T>> {
        if !(self.rect.min.x..self.rect.max.x).contains(&pos.x) {
            return None;
        }
        if !(self.rect.min.y..self.rect.max.y).contains(&pos.y) {
            return None
        }
        Some(pos)
    }

    #[inline(always)]
    fn displacement(&self, pos1: Pos<Self::Coord>, pos2: Pos<Self::Coord>) -> (Self::Coord, Self::Coord) {
        (pos2.x - pos1.x, pos2.y - pos1.y)
    }
}

impl ToLatticeBoundary for FixedBoundary<f32> {
    type LatticeBoundary = FixedBoundary<isize>;
    type Error = RectConversionError;

    fn to_lattice_boundary(&self) -> Result<FixedBoundary<isize>, Self::Error> {
        Ok(FixedBoundary::new(Rect::try_from(self.rect.clone())?))
    }
}

#[derive(Clone)]
pub struct SafePeriodicBoundary<T> {
    rect: Rect<T>
}

impl<T> SafePeriodicBoundary<T> {
    pub fn new(rect: Rect<T>) -> Self {
        Self { rect }
    }
}

impl<T> Boundary for SafePeriodicBoundary<T>
where
    T: Num + Copy + Euclid + From<u8> {
    type Coord = T;

    #[inline(always)]
    fn rect(&self) -> &Rect<T> {
        &self.rect
    }

    #[inline(always)]
    fn valid_pos(&self, pos: Pos<Self::Coord>) -> Option<Pos<Self::Coord>> {
        self.periodic_valid_pos(pos)
    }

    #[inline(always)]
    fn displacement(&self, pos1: Pos<Self::Coord>, pos2: Pos<Self::Coord>) -> (Self::Coord, Self::Coord) {
        self.periodic_displacement(pos1, pos2)
    }
}

impl<T> PeriodicBoundary for SafePeriodicBoundary<T>
where
    T: Num + Copy + Euclid + From<u8> {
    #[inline(always)]
    fn wrap_scalar(val: Self::Coord, min: Self::Coord, max: Self::Coord) -> Self::Coord {
        let range = max - min;
        let mut offset = val - min;
        offset = offset.rem_euclid(&range);
        min + offset
    }
}

impl ToLatticeBoundary for SafePeriodicBoundary<f32> {
    type LatticeBoundary = SafePeriodicBoundary<isize>;
    type Error = RectConversionError;

    fn to_lattice_boundary(&self) -> Result<SafePeriodicBoundary<isize>, Self::Error> {
        Ok(SafePeriodicBoundary::new(Rect::try_from(self.rect.clone())?))
    }
}

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

impl<T> UnsafePeriodicBoundary<T> {
    pub fn new(rect: Rect<T>) -> Self {
        Self { rect }
    }
}

impl<T> Boundary for UnsafePeriodicBoundary<T>
where
    T: Num + Copy + Euclid + From<u8> + PartialOrd {
    type Coord = T;

    #[inline(always)]
    fn rect(&self) -> &Rect<T> {
        &self.rect
    }

    #[inline(always)]
    fn valid_pos(&self, pos: Pos<Self::Coord>) -> Option<Pos<Self::Coord>> {
        self.periodic_valid_pos(pos)
    }

    #[inline(always)]
    fn displacement(&self, pos1: Pos<Self::Coord>, pos2: Pos<Self::Coord>) -> (Self::Coord, Self::Coord) {
        self.periodic_displacement(pos1, pos2)
    }
}

impl<T> PeriodicBoundary for UnsafePeriodicBoundary<T>
where
    T: Num + Copy + Euclid + From<u8> + PartialOrd {
    #[inline(always)]
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

impl ToLatticeBoundary for UnsafePeriodicBoundary<f32> {
    type LatticeBoundary = UnsafePeriodicBoundary<isize>;
    type Error = RectConversionError;

    fn to_lattice_boundary(&self) -> Result<UnsafePeriodicBoundary<isize>, Self::Error> {
        Ok(UnsafePeriodicBoundary::new(Rect::try_from(self.rect.clone())?))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::positional::pos::Pos;

    fn rect_10x10() -> Rect<isize> {
        Rect::new(Pos::new(0, 0), Pos::new(10, 10))
    }

    #[test]
    fn test_fixed_valid() {
        let fixed = FixedBoundary::new(rect_10x10());

        // In bounds
        assert_eq!(fixed.valid_pos(Pos::new(5, 5)), Some(Pos::new(5, 5)));

        // Out of bounds
        assert_eq!(fixed.valid_pos(Pos::new(-1, 5)), None);
        assert_eq!(fixed.valid_pos(Pos::new(5, -1)), None);
        assert_eq!(fixed.valid_pos(Pos::new(10, 5)), None);
        assert_eq!(fixed.valid_pos(Pos::new(5, 10)), None);
    }

    #[test]
    fn test_fixed_displacement() {
        let fixed = FixedBoundary::new(rect_10x10());
        let d = fixed.displacement(Pos::new(3, 3), Pos::new(6, 1));
        assert_eq!(d, (3, -2));
    }

    #[test]
    fn test_safe_periodic_valid() {
        let per = SafePeriodicBoundary::new(rect_10x10());

        // Valid values wrap around
        assert_eq!(per.valid_pos(Pos::new(-1, -1)), Some(Pos::new(9, 9)));
        assert_eq!(per.valid_pos(Pos::new(10, 10)), Some(Pos::new(0, 0)));
        assert_eq!(per.valid_pos(Pos::new(30, 30)), Some(Pos::new(0, 0)));
    }

    #[test]
    fn test_safe_periodic_displacement() {
        let per = SafePeriodicBoundary::new(rect_10x10());

        let d = per.displacement(Pos::new(9, 0), Pos::new(1, 0)); // wrapped
        assert_eq!(d, (2, 0)); // shortest dx is -8 → wraps to +2

        let d2 = per.displacement(Pos::new(1, 0), Pos::new(9, 0)); // reverse
        assert_eq!(d2, (-2, 0));
    }

    #[test]
    fn test_unsafe_periodic_valid() {
        let unsafe_per = UnsafePeriodicBoundary::new(rect_10x10());

        // Near bounds: valid
        assert_eq!(unsafe_per.valid_pos(Pos::new(10, 0)), Some(Pos::new(0, 0)));
        assert_eq!(unsafe_per.valid_pos(Pos::new(-1, 0)), Some(Pos::new(9, 0)));

        // Far outside bounds: wrong result, but not panicked
        let pos = unsafe_per.valid_pos(Pos::new(30, 30)).unwrap();
        assert_ne!(pos, Pos::new(0, 0)); // Incorrect but expected under "unsafe"
    }

    #[test]
    fn test_safe_periodic_wrap_scalar() {
        let wrapped = SafePeriodicBoundary::<isize>::wrap_scalar(11, 0, 10);
        assert_eq!(wrapped, 1);

        let wrapped_neg = SafePeriodicBoundary::<isize>::wrap_scalar(-1, 0, 10);
        assert_eq!(wrapped_neg, 9);
    }

    #[test]
    fn test_unsafe_periodic_wrap_scalar() {
        let wrap = UnsafePeriodicBoundary::<isize>::wrap_scalar(11, 0, 10);
        assert_eq!(wrap, 1);

        let wrap_neg = UnsafePeriodicBoundary::<isize>::wrap_scalar(-1, 0, 10);
        assert_eq!(wrap_neg, 9);

        let in_bounds = UnsafePeriodicBoundary::<isize>::wrap_scalar(5, 0, 10);
        assert_eq!(in_bounds, 5);
    }
}
