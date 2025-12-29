//! Contains logic associated with boundary conditions and position validation.

use crate::positional::pos::Pos;
use crate::positional::rect::Rect;
use num::traits::Euclid;
use num::Num;
use std::ops::Sub;

/// A conjugate pair of continuous/discrete boundaries that have the same boundary conditions.
#[derive(Clone, Debug, PartialEq)]
pub struct Boundaries<B: ToLatticeBoundary> {
    /// The conjugate continuous boundary type.
    pub boundary: B,
    /// The conjugate discrete boundary type.
    pub lattice_boundary: B::LatticeBoundary,
}

impl<B: ToLatticeBoundary> Boundaries<B> {
    /// Makes a new boundary pair from a continuous boundary type.
    pub fn new(bound: B) -> Self {
        Self {
            lattice_boundary: bound.to_lattice_boundary(),
            boundary: bound,
        }
    }
}

/// Types that implement boundary conditions to validate positions in space.
pub trait Boundary {
    /// Coordinate system associated with the boundaries.
    type Coord;

    /// Expose the size of the boundary as a [Rect].
    fn rect(&self) -> &Rect<Self::Coord>;

    /// Validates that positions are in bounds.
    ///
    /// With fixed boundary conditions, that means filtering invalid positions.
    fn valid_pos(&self, pos: Pos<Self::Coord>) -> Option<Pos<Self::Coord>>;

    /// Validates a series of positions using [Boundary::valid_pos()].
    fn valid_positions(
        &self,
        positions: impl Iterator<Item = Pos<Self::Coord>>
    ) -> impl Iterator<Item = Pos<Self::Coord>> {
        positions.filter_map(|pos| self.valid_pos(pos))
    }

    /// Calculates the minimum displacement between two positions by taking into account boundary conditions.
    fn displacement(&self, pos1: Pos<Self::Coord>, pos2: Pos<Self::Coord>) -> (Self::Coord, Self::Coord);
}

/// Private trait that provides default implementations for periodic boundary types.
trait PeriodicBoundary: Boundary
where
    Self::Coord: From<u8> + Num + Copy + Euclid {
    fn periodic_displacement(&self, pos1: Pos<Self::Coord>, pos2: Pos<Self::Coord>) -> (Self::Coord, Self::Coord) {
        let two = Self::Coord::from(2);
        let dx = ((pos2.x - pos1.x + self.rect().width() / two)
            .rem_euclid(&self.rect().width())) - self.rect().width() / two;
        let dy = ((pos2.y - pos1.y + self.rect().height() / two)
            .rem_euclid(&self.rect().height())) - self.rect().height() / two;
        (dx, dy)
    }

    fn periodic_valid_pos(&self, pos: Pos<Self::Coord>) -> Option<Pos<Self::Coord>> {
        let p = Pos::new(
            Self::wrap_scalar(pos.x, self.rect().min.x, self.rect().max.x),
            Self::wrap_scalar(pos.y, self.rect().min.y, self.rect().max.y)
        );
        Some(p)
    }

    fn wrap_scalar(val: Self::Coord, min: Self::Coord, max: Self::Coord) -> Self::Coord;
}

/// This trait is used to define conjugate (continuous, discrete) boundary conditions used to make [Boundaries].
pub trait ToLatticeBoundary: Boundary {
    /// Conjugate discrete boundary type of this continuous boundary type.
    type LatticeBoundary: Boundary<Coord = isize>;
    /// Returns an instance of the conjugate [ToLatticeBoundary::LatticeBoundary] type.
    fn to_lattice_boundary(&self) -> Self::LatticeBoundary;
}

/// Implementation of fixed (truncated) boundary conditions.
#[derive(Clone, Debug, PartialEq)]
pub struct FixedBoundary<T> {
    rect: Rect<T>
}

impl<T> FixedBoundary<T> {
    /// Makes a new boundary with size determined by `rect`.
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

    fn valid_pos(&self, pos: Pos<T>) -> Option<Pos<T>> {
        if !(self.rect.min.x..self.rect.max.x).contains(&pos.x) {
            return None;
        }
        if !(self.rect.min.y..self.rect.max.y).contains(&pos.y) {
            return None
        }
        Some(pos)
    }

    fn displacement(&self, pos1: Pos<Self::Coord>, pos2: Pos<Self::Coord>) -> (Self::Coord, Self::Coord) {
        (pos2.x - pos1.x, pos2.y - pos1.y)
    }
}

impl ToLatticeBoundary for FixedBoundary<f32> {
    type LatticeBoundary = FixedBoundary<isize>;

    fn to_lattice_boundary(&self) -> FixedBoundary<isize> {
        FixedBoundary::new(self.rect.to_isize())
    }
}

/// Mathematically sound (albeit slow) implementation of periodic boundary conditions.
#[derive(Clone, Debug, PartialEq)]
pub struct SafePeriodicBoundary<T> {
    rect: Rect<T>
}

impl<T> SafePeriodicBoundary<T> {
    /// Makes a new boundary with size determined by `rect`.
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

    fn valid_pos(&self, pos: Pos<Self::Coord>) -> Option<Pos<Self::Coord>> {
        self.periodic_valid_pos(pos)
    }

    fn displacement(&self, pos1: Pos<Self::Coord>, pos2: Pos<Self::Coord>) -> (Self::Coord, Self::Coord) {
        self.periodic_displacement(pos1, pos2)
    }
}

impl<T> PeriodicBoundary for SafePeriodicBoundary<T>
where
    T: Num + Copy + Euclid + From<u8> {
    fn wrap_scalar(val: Self::Coord, min: Self::Coord, max: Self::Coord) -> Self::Coord {
        let range = max - min;
        let mut offset = val - min;
        offset = offset.rem_euclid(&range);
        min + offset
    }
}

impl ToLatticeBoundary for SafePeriodicBoundary<f32> {
    type LatticeBoundary = SafePeriodicBoundary<isize>;

    fn to_lattice_boundary(&self) -> SafePeriodicBoundary<isize> {
        SafePeriodicBoundary::new(self.rect.to_isize())
    }
}

/// This type can only validate positions that are at most one [UnsafePeriodicBoundary::rect]
/// away from the boundaries (in either x or y direction).
///
/// <div class="warning">
///
/// Only use when you are confident that all input positions are close to the boundary.
///
/// </div>
#[derive(Clone, Debug, PartialEq)]
pub struct UnsafePeriodicBoundary<T> {
    rect: Rect<T>
}

impl<T> UnsafePeriodicBoundary<T> {
    /// Makes a new boundary with size determined by `rect`.
    pub fn new(rect: Rect<T>) -> Self {
        Self { rect }
    }
}

impl<T> Boundary for UnsafePeriodicBoundary<T>
where
    T: Num + Copy + Euclid + From<u8> + PartialOrd {
    type Coord = T;

    fn rect(&self) -> &Rect<T> {
        &self.rect
    }

    fn valid_pos(&self, pos: Pos<Self::Coord>) -> Option<Pos<Self::Coord>> {
        self.periodic_valid_pos(pos)
    }

    fn displacement(&self, pos1: Pos<Self::Coord>, pos2: Pos<Self::Coord>) -> (Self::Coord, Self::Coord) {
        self.periodic_displacement(pos1, pos2)
    }
}

impl<T> PeriodicBoundary for UnsafePeriodicBoundary<T>
where
    T: Num + Copy + Euclid + From<u8> + PartialOrd {
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

    fn to_lattice_boundary(&self) -> UnsafePeriodicBoundary<isize> {
        UnsafePeriodicBoundary::new(self.rect.to_isize())
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
