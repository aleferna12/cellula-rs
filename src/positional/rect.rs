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
    + Copy
{
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
}