use crate::positional::pos::Pos;
use num::{Integer, Num};
use std::ops::AddAssign;

#[derive(Clone, Debug)]
pub struct Rect<T> {
    pub min: Pos<T>,
    pub max: Pos<T>
}

impl<T> Rect<T>
where
    T: Num
    + Copy {
    pub fn new(min: Pos<T>, max: Pos<T>) -> Self {
        Self{ min, max }
    }

    pub fn width(&self) -> T {
        self.max.x - self.min.x
    }

    pub fn height(&self) -> T {
        self.max.y - self.min.y
    }

    pub fn area(&self) -> T {
        self.width() * self.height()
    }
}

impl<T> Rect<T>
where
    T: Integer 
        + AddAssign
        + Copy {
    pub fn iter_positions(&self) -> RectAreaIt<T> {
        RectAreaIt::new(self.clone())
    }
}

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
        if self.curr.y >= self.rect.max.y {
            return None;
        }
        let ret_pos = self.curr;
        if self.curr.x < self.rect.max.x - T::one() {
            self.curr.x += T::one();
        } else {
            self.curr.x = self.rect.min.x;
            self.curr.y += T::one();
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
            Pos::new(0, 0), Pos::new(1, 0),
            Pos::new(0, 1), Pos::new(1, 1),
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