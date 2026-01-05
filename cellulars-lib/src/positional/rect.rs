//! Contains logic associated with[`Rect`].

use crate::positional::pos::Pos;
use num::{Integer, Num};
use std::ops::AddAssign;

/// A rectangle defined by two corners.
#[derive(Clone, Debug, PartialEq)]
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

impl Rect<f32> {
    /// Rounds the rectangle's coordinates to their nearest integer value.
    pub fn round(&self) -> Rect<f32> {
        Rect::new(self.min.round(), self.max.round())
    }

    /// Casts the rect's coordinate type to [isize] (by truncating).
    pub fn to_isize(&self) -> Rect<isize> {
        Rect::new(
            self.min.to_isize(),
            self.max.to_isize(),
        )
    }

    /// Casts the rect's coordinate type to [usize] (by truncating).
    pub fn to_usize(&self) -> Rect<usize> {
        Rect::new(
            self.min.to_usize(),
            self.max.to_usize(),
        )
    }
}

impl Rect<usize> {
    /// Casts the rect's coordinate type to [isize].
    pub fn to_isize(&self) -> Rect<isize> {
        Rect::new(
            self.min.to_isize(),
            self.max.to_isize(),
        )
    }

    /// Casts the rect's coordinate type to [f32].
    pub fn to_f32(&self) -> Rect<f32> {
        Rect::new(
            self.min.to_f32(),
            self.max.to_f32(),
        )
    }
}

impl Rect<isize> {
    /// Casts the rect's coordinate type to [usize].
    pub fn to_usize(&self) -> Rect<usize> {
        Rect::new(
            self.min.to_usize(),
            self.max.to_usize(),
        )
    }

    /// Casts the rect's coordinate type to [f32].
    pub fn to_f32(&self) -> Rect<f32> {
        Rect::new(
            self.min.to_f32(),
            self.max.to_f32(),
        )
    }
}

/// Iterator over all discrete positions contained in a [Rect].
///
/// The iterator range is inclusive on both ends (position \[width, height\] is included for example).
#[derive(Clone, Debug, PartialEq)]
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