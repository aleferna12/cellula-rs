use std::ops::AddAssign;
use num::{Integer, Num};

#[derive(Debug)]
pub enum EdgeError {
    SamePosition,
    NotNeighbours
}

/// 2D position in space.
#[derive(PartialEq, Eq, Copy, Clone, Debug)]
pub struct Pos2D<T> {
    pub x: T,
    pub y: T
}

impl<T> Pos2D<T> {
    pub fn new(x: T, y: T) -> Self {
        Self { x, y }
    }
}

impl<T> From<(T, T)> for Pos2D<T> {
    fn from(value: (T, T)) -> Self {
        Pos2D::<T>::new(value.0, value.1)
    }
}

impl Pos2D<usize> {
    pub(crate) fn pack_u32(self) -> u32 {
        ((self.x as u32) << 16) | self.y as u32
    }

    pub fn row_major(self, height: usize) -> usize {
        self.x * height + self.y
    }
}

impl From<Pos2D<usize>> for Pos2D<isize> {
    fn from(value: Pos2D<usize>) -> Self {
        Pos2D::new(value.x as isize, value.y as isize)
    }
}

impl From<Pos2D<isize>> for Pos2D<usize> {
    fn from(value: Pos2D<isize>) -> Self {
        let message = "overflow when translating position from general to lattice coordinates";
        Pos2D::new(
            value.x.try_into().expect(message), 
            value.y.try_into().expect(message)
        )
    }
}

#[derive(Clone, Debug)]
pub struct Rect<T> {
    pub min: Pos2D<T>,
    pub max: Pos2D<T>
}
impl<T> Rect<T>
where
    T: Num
    + Copy
{
    pub fn new(min: Pos2D<T>, max: Pos2D<T>) -> Self {
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
        RectAreaIt::new(self)
    }
}

pub struct RectAreaIt<'a, T> {
    curr: Pos2D<T>,
    rect: &'a Rect<T>
}
impl<'a, T: Copy> RectAreaIt<'a, T> {
    fn new(rect: &'a Rect<T>) -> Self {
        Self {
            curr: rect.min,
            rect
        }
    }
}

impl<T> Iterator for RectAreaIt<'_, T>
where
    T: Copy
    + Integer
    + AddAssign {
    type Item = Pos2D<T>;

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