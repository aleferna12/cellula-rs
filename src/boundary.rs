use crate::pos::{Pos2D, Rect};

pub trait Boundary {
    type Coord: Copy;
    
    /// Expose the boundary as a `Rect`.
    fn rect(&self) -> &Rect<Self::Coord>;
    
    fn inbounds(&self, pos: Pos2D<Self::Coord>) -> bool;

    /// Validates that positions are in bounds.
    ///
    /// With fixed boundary conditions, that means filtering invalid positions.
    fn validate_positions(
        &self,
        pos_it: impl Iterator<Item = Pos2D<Self::Coord>>
    ) -> impl Iterator<Item = Pos2D<Self::Coord>>;
}

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

    fn inbounds(&self, pos: Pos2D<T>) -> bool {
        (self.rect.min.x..self.rect.max.x).contains(&pos.x) 
            && (self.rect.min.y..self.rect.max.y).contains(&pos.y)
    }

    fn validate_positions(
        &self,
        pos_it: impl Iterator<Item = Pos2D<T>>
    ) -> impl Iterator<Item = Pos2D<T>> {
        pos_it.filter(|pos| { 
            self.inbounds(*pos) 
        })
    }
}