//! Contains logic associated with [Rect].

use crate::positional::pos::Pos;
use num::{Integer, Num, ToPrimitive};
use std::ops::AddAssign;
use thiserror::Error;

/// A rectangle defined by two corners.
#[derive(Clone, Debug)]
pub struct Rect<T> {
    /// Bottom-left corner.
    pub min: Pos<T>,
    /// Upper-right corner.
    pub max: Pos<T>
}

impl<T> Rect<T>
where
    T: Num + Copy {
    /// Makes a new rectangle spanning from positions `min` to `max`.
    pub fn new(min: Pos<T>, max: Pos<T>) -> Self {
        Self{ min, max }
    }

    /// Returns the width of the rectangle.
    pub fn width(&self) -> T {
        self.max.x - self.min.x
    }

    /// Returns the height of the rectangle.
    pub fn height(&self) -> T {
        self.max.y - self.min.y
    }

    /// Returns the area of the rectangle.
    pub fn area(&self) -> T {
        self.width() * self.height()
    }
}

impl<T> Rect<T>
where
    T: Integer
        + AddAssign
        + Copy {
    /// Iterates over all discrete positions contained in the rectangle.
    ///
    /// The iterator range is inclusive on both ends (position \[width, height\] is included for example).
    pub fn iter_positions(&self) -> RectAreaIt<T> {
        RectAreaIt::new(self.clone())
    }
}

/// Iterator over all discrete positions contained in a [Rect].
///
/// The iterator range is inclusive on both ends (position \[width, height\] is included for example).
pub struct RectAreaIt<T> {
    curr: Pos<T>,
    rect: Rect<T>
}

impl<T: Copy> RectAreaIt<T> {
    fn new(rect: Rect<T>) -> Self {
        Self {
            curr: rect.min,
            rect
        }
    }
}

impl<T> Iterator for RectAreaIt<T>
where
    T: Copy
        + Integer
        + AddAssign {
    type Item = Pos<T>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr.x >= self.rect.max.x {
            return None;
        }
        let ret_pos = self.curr;
        if self.curr.y < self.rect.max.y - T::one() {
            self.curr.y += T::one();
        } else {
            self.curr.y = self.rect.min.y;
            self.curr.x += T::one();
        }
        Some(ret_pos)
    }
}

// TODO: make these explicit (like in Pos)
impl TryFrom<Rect<f32>> for Rect<isize> {
    type Error = RectConversionError;

    fn try_from(value: Rect<f32>) -> Result<Self, Self::Error> {
        Ok(Rect::new(
            Pos::new(
                value.min.x.to_isize().ok_or(Self::Error {})?,
                value.min.y.to_isize().ok_or(Self::Error {})?,
            ),
            Pos::new(
                value.max.x.to_isize().ok_or(Self::Error {})?,
                value.max.y.to_isize().ok_or(Self::Error {})?,
            )
        ))
    }
}

impl TryFrom<Rect<f32>> for Rect<usize> {
    type Error = RectConversionError;

    fn try_from(value: Rect<f32>) -> Result<Self, Self::Error> {
        Ok(Rect::new(
            Pos::new(
                value.min.x.to_usize().ok_or(Self::Error {})?,
                value.min.y.to_usize().ok_or(Self::Error {})?,
            ),
            Pos::new(
                value.max.x.to_usize().ok_or(Self::Error {})?,
                value.max.y.to_usize().ok_or(Self::Error {})?,
            )
        ))
    }
}

/// Error raised when the coordinate type of [Rect] cannot be converted.
#[derive(Error, Debug)]
#[error("failed to convert `Rect`")]
pub struct RectConversionError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_rect_area() {
        let r = Rect::<usize>::new((0, 0).into(), (10, 10).into());
        let v: Vec<_> = r.iter_positions().collect();
        assert_eq!(r.area(), v.len())
    }
    
    #[test]
    fn test_rect_dimensions() {
        let r = Rect::new(Pos::new(2, 3), Pos::new(7, 9));
        assert_eq!(r.width(), 5);  // 7 - 2
        assert_eq!(r.height(), 6); // 9 - 3
        assert_eq!(r.area(), 30);  // 5 * 6
    }

    #[test]
    fn test_iter_positions_order_and_coverage() {
        let rect = Rect::new(Pos::new(0, 0), Pos::new(2, 2));
        let expected = vec![
            Pos::new(0, 0), Pos::new(0, 1),
            Pos::new(1, 0), Pos::new(1, 1),
        ];
        let result: Vec<_> = rect.iter_positions().collect();
        assert_eq!(result, expected);
    }

    #[test]
    fn test_iter_empty_if_min_equals_max() {
        let rect = Rect::new(Pos::new(5, 5), Pos::new(5, 5));
        let result: Vec<_> = rect.iter_positions().collect();
        assert_eq!(result.len(), 0);
    }

    #[test]
    fn test_area_matches_iter_count() {
        let r = Rect::<usize>::new((0, 0).into(), (10, 3).into());
        let iter_count = r.iter_positions().count();
        assert_eq!(r.area(), iter_count);
    }
}