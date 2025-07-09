use num::traits::Euclid;
use num::Num;
use crate::cell::Cell;
use crate::positional::pos::{AngularProjection, Pos2D};
use crate::positional::rect::Rect;

pub trait Boundary {
    type Coord;
    
    /// Expose the boundary as a `Rect`.
    fn rect(&self) -> &Rect<Self::Coord>;
    
    /// Validates that positions are in bounds.
    ///
    /// With fixed boundary conditions, that means filtering invalid positions.
    fn valid_pos(&self, pos: Pos2D<Self::Coord>) -> Option<Pos2D<Self::Coord>>;
    
    fn valid_positions(
        &self,
        positions: impl Iterator<Item = Pos2D<Self::Coord>>
    ) -> impl Iterator<Item = Pos2D<Self::Coord>> {
        positions.filter_map(|pos| self.valid_pos(pos))
    }
}

pub trait LatticeBoundary: Boundary<Coord = isize> {
    fn shift_cell_center(cell: &mut Cell, pos: Pos2D<usize>, width: usize, height: usize, add: bool);
    // TODO!: this can make FixedBoundary quite a lot faster
    // fn delta_angle();
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

impl<T: PartialOrd + Copy> Boundary for FixedBoundary<T> {
    type Coord = T;

    fn rect(&self) -> &Rect<T> {
        &self.rect
    }

    fn valid_pos(&self, pos: Pos2D<T>) -> Option<Pos2D<T>> {
        if !(self.rect.min.x..self.rect.max.x).contains(&pos.x) {
            return None;
        }
        if !(self.rect.min.y..self.rect.max.y).contains(&pos.y) {
            return None
        }
        Some(pos)
    }
}

impl LatticeBoundary for FixedBoundary<isize> {
    fn shift_cell_center(cell: &mut Cell, pos: Pos2D<usize>, _width: usize, _height: usize, add: bool) {
        let shift = if add { 1. } else { -1. };
        let area = cell.area as f32;
        
        cell.center.pos = Pos2D::new(
            cell.center.pos.x + shift * (pos.x as f32 - cell.center.pos.x) / area,
            cell.center.pos.y + shift * (pos.y as f32 - cell.center.pos.y) / area
        );
    }
}

#[derive(Clone)]
pub struct PeriodicBoundary<T> {
    rect: Rect<T>
}

impl<T> PeriodicBoundary<T> {
    pub fn new(rect: Rect<T>) -> Self {
        Self { rect }
    }
}

impl<T> PeriodicBoundary<T>
where
    T: Copy + Num + Euclid {
    fn wrap_scalar(val: T, min: T, max: T) -> T {
        let range = max - min;
        let mut offset = val - min;
        offset = offset.rem_euclid(&range);
        min + offset
    }
}

impl<T> Boundary for PeriodicBoundary<T>
where
    T: Copy + Num + Euclid {
    type Coord = T;

    fn rect(&self) -> &Rect<T> {
        &self.rect
    }

    fn valid_pos(&self, pos: Pos2D<Self::Coord>) -> Option<Pos2D<Self::Coord>> {
        let p = Pos2D::new(
            Self::wrap_scalar(pos.x, self.rect.min.x, self.rect.max.x),
            Self::wrap_scalar(pos.y, self.rect.min.y, self.rect.max.y)
        );
        Some(p)
    }
}

impl LatticeBoundary for PeriodicBoundary<isize> {
    fn shift_cell_center(cell: &mut Cell, pos: Pos2D<usize>, width: usize, height: usize, add: bool) {
        let shift = if add { 1. } else { -1. };

        let proj = AngularProjection::from_pos(Pos2D::new(pos.x as f32, pos.y as f32), width, height);
        let sum_proj = &mut cell.center.projection;
        sum_proj.x_sin += shift * proj.x_sin;
        sum_proj.x_cos += shift * proj.x_cos;
        sum_proj.y_sin += shift * proj.y_sin;
        sum_proj.y_cos += shift * proj.y_cos;
        
        cell.center.pos = Pos2D::from_projection(&cell.center.projection, width, height);
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

impl<T> UnsafePeriodicBoundary<T>
where
    T: Copy
        + Num 
        + PartialOrd {
    pub fn new(rect: Rect<T>) -> Self {
        Self { rect }
    }

    pub fn wrap_scalar(&self, val: T, min: T, max: T) -> T {
        if val < min {
            max - (min - val)
        } else if val >= max {
            min + (val - max)
        } else {
            val
        }
    }
}

impl LatticeBoundary for UnsafePeriodicBoundary<isize> {
    fn shift_cell_center(cell: &mut Cell, pos: Pos2D<usize>, width: usize, height: usize, add: bool) {
        PeriodicBoundary::shift_cell_center(cell, pos, width, height, add)
    }
}

impl<T> Boundary for UnsafePeriodicBoundary<T>
where
    T: Copy + Num + Euclid + PartialOrd {
    type Coord = T;

    fn rect(&self) -> &Rect<Self::Coord> {
        &self.rect
    }

    fn valid_pos(&self, pos: Pos2D<Self::Coord>) -> Option<Pos2D<Self::Coord>> {
        Some(Pos2D::new(
            self.wrap_scalar(pos.x, self.rect().min.x, self.rect().max.x),
            self.wrap_scalar(pos.y, self.rect().min.y, self.rect().max.y)
        ))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_periodic_boundary() {
        let per = PeriodicBoundary::new(Rect::new((0, 0).into(), (10, 10).into()));
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