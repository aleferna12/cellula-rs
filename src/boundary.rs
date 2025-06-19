use crate::pos::{Pos2D, Rect};

pub trait LatticeBoundary {
    /// Expose the boundary as a `Rect`.
    fn rect(&self) -> &Rect<usize>;

    /// Validates that positions are in bounds.
    ///
    /// With fixed boundary conditions, that means filtering invalid positions.
    fn validate_positions<S>(&self, pos_it: S) -> impl Iterator<Item = Pos2D<usize>>
    where
        S: Iterator<Item = Pos2D<usize>>;
}

pub struct FixedBoundary {
    rect: Rect<usize>
}
impl FixedBoundary {
    pub fn new(width: usize, height: usize) -> Self {
        Self { rect: Rect::new((0, 0).into(), (width, height).into()) }
    }
}

impl LatticeBoundary for FixedBoundary {
    fn rect(&self) -> &Rect<usize> {
        &self.rect
    }

    fn validate_positions<S: Iterator<Item = Pos2D<usize>>>(
        &self,
        pos_it: S
    ) -> impl Iterator<Item = Pos2D<usize>> {
        pos_it.filter(move |pos| { self.rect.inbounds(*pos) })
    }
}