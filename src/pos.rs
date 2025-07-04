use std::f32::consts::TAU;
use std::ops::AddAssign;
use num::{Integer, Num};

/// 2D position in space.
#[derive(PartialEq, Eq, Copy, Clone, Debug)]
#[derive(Hash)]
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

impl Pos2D<f32> {
    /// Unwraps the `AngularProjection` into a position.
    pub(crate) fn from_projection(proj: &AngularProjection, width: usize, height: usize) -> Self {
        let (angle_x, angle_y) = proj.angles();
        Self {
            x: width as f32 * angle_x / TAU,
            y: height as f32 * angle_y / TAU
        }
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

#[derive(Debug, Clone)]
pub(crate) struct AngularProjection {
    pub(crate) x_sin: f32,
    pub(crate) x_cos: f32,
    pub(crate) y_sin: f32,
    pub(crate) y_cos: f32
}

impl AngularProjection {
    pub(crate) fn from_pos(pos: Pos2D<f32>, width: usize, height: usize) -> Self {
        let cx = TAU * pos.x / width as f32;
        let cy = TAU * pos.y  / height as f32;
        Self {
            x_sin: cx.sin(),
            x_cos: cx.cos(),
            y_sin: cy.sin(),
            y_cos: cy.cos()
        }
    }
    
    /// Returns the angles associated with this projection.
    pub(crate) fn angles(&self) -> (f32, f32) {
        (self.x_sin.atan2(self.x_cos).rem_euclid(TAU), self.y_sin.atan2(self.y_cos).rem_euclid(TAU))
    }
    
    // TODO!: this can be optimised to not require atan2 
    //  (and has a significant impact in performance due to hot call in CA)
    pub(crate) fn delta_angles(&self, other: &AngularProjection) -> (f32, f32) {
        // This avoids unecessary multiple calls to rem_euclid
        let x = (self.x_sin * other.x_cos - self.x_cos * other.x_sin)
            .atan2(self.x_cos * other.x_cos + self.x_sin * other.x_sin);
        let y = (self.y_sin * other.y_cos - self.y_cos * other.y_sin)
            .atan2(self.y_cos * other.y_cos + self.y_sin * other.y_sin);
        (x, y)
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
        RectAreaIt::new(self.clone())
    }
}

pub struct RectAreaIt<T> {
    curr: Pos2D<T>,
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