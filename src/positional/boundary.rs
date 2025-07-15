use crate::positional::pos::Pos;
use crate::positional::rect::Rect;
use num::traits::Euclid;
use num::Num;

pub trait Boundary {
    type Coord;
    
    /// Expose the boundary as a `Rect`.
    fn rect(&self) -> &Rect<Self::Coord>;
    
    /// Validates that positions are in bounds.
    ///
    /// With fixed boundary conditions, that means filtering invalid positions.
    fn valid_pos(&self, pos: Pos<Self::Coord>) -> Option<Pos<Self::Coord>>;
    
    fn valid_positions(
        &self,
        positions: impl Iterator<Item = Pos<Self::Coord>>
    ) -> impl Iterator<Item = Pos<Self::Coord>> {
        positions.filter_map(|pos| self.valid_pos(pos))
    }
}

pub trait LatticeBoundary: Boundary<Coord = isize> {
    /// Shifts the center of mass (`com`) by `pos` taking into account its `mass`.
    /// 
    /// `add` determined whether to add or remove `pos` from `com`.
    fn shift_center_of_mass(
        &self,
        com: Pos<f32>,
        pos: Pos<f32>,
        mass: f32,
        add: bool
    ) -> Option<Pos<f32>> {
        let shift = if add { 1. } else { -1. };
        let new_mass = mass + shift;
        if new_mass <= 0.0 {
            return None;
        }

        let w = self.rect().width() as f32;
        let h = self.rect().height() as f32;
        // TODO! generalise displacement
        let dx = ((pos.x - com.x + w / 2.).rem_euclid(w)) - w / 2.;
        let dy = ((pos.y - com.y + h / 2.).rem_euclid(h)) - h / 2.;
        let x = (com.x + dx * shift / new_mass).rem_euclid(w);
        let y = (com.y + dy * shift / new_mass).rem_euclid(h);
        Some(Pos::new(
            x, y
        ))
    }

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

    fn valid_pos(&self, pos: Pos<T>) -> Option<Pos<T>> {
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
    fn shift_center_of_mass(
        &self,
        com: Pos<f32>,
        pos: Pos<f32>,
        mass: f32,
        add: bool
    ) -> Option<Pos<f32>> {
        let shift = if add { 1. } else { -1. };
        let new_mass = mass + shift;
        if new_mass <= 0.0 {
            return None;
        }
        Some(Pos::new(
            com.x + shift * (pos.x - com.x) / new_mass,
            com.y + shift * (pos.y - com.y) / new_mass
        ))
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

    fn valid_pos(&self, pos: Pos<Self::Coord>) -> Option<Pos<Self::Coord>> {
        let p = Pos::new(
            Self::wrap_scalar(pos.x, self.rect.min.x, self.rect.max.x),
            Self::wrap_scalar(pos.y, self.rect.min.y, self.rect.max.y)
        );
        Some(p)
    }
}

impl LatticeBoundary for PeriodicBoundary<isize> {}

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

impl<T> Boundary for UnsafePeriodicBoundary<T>
where
    T: Copy + Num + Euclid + PartialOrd {
    type Coord = T;

    fn rect(&self) -> &Rect<Self::Coord> {
        &self.rect
    }

    fn valid_pos(&self, pos: Pos<Self::Coord>) -> Option<Pos<Self::Coord>> {
        Some(Pos::new(
            self.wrap_scalar(pos.x, self.rect().min.x, self.rect().max.x),
            self.wrap_scalar(pos.y, self.rect().min.y, self.rect().max.y)
        ))
    }
}

impl LatticeBoundary for UnsafePeriodicBoundary<isize> {}

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