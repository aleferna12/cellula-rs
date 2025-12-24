//! Contains logic associated with [Pos].

use num::ToPrimitive;

/// We use this everytime we expect a position conversion to work.
pub(crate) const CONV_ERROR: &str = "failed to convert between positional coordinates";

/// A 2D position in space.
#[derive(PartialEq, Eq, Copy, Clone, Debug)]
#[derive(Hash)]
pub struct Pos<T> {
    /// X component of the position.
    pub x: T,
    /// Y component of the position.
    pub y: T
}

impl<T> Pos<T> {
    /// Makes a new position with x, y coordinates.
    pub fn new(x: T, y: T) -> Self {
        Self { x, y }
    }
}

impl Pos<usize> {
    pub(crate) fn pack_u32(self) -> u32 {
        ((self.x as u32) << 16) | self.y as u32
    }

    /// Compresses the position to 1D using column-major ordering.
    pub fn col_major(self, height: usize) -> usize {
        self.x * height + self.y
    }

    /// Converts the position's coordinate type to [isize].
    ///
    /// Fails under the same conditions that make [num::ToPrimitive](ToPrimitive) fail.
    pub fn to_isize(self) -> Option<Pos<isize>> {
        Some(Pos::new(self.x.to_isize()?, self.y.to_isize()?))
    }

    /// Converts the position's coordinate type to [f32].
    ///
    /// Fails under the same conditions that make [num::ToPrimitive](ToPrimitive) fail.
    pub fn to_f32(self) -> Option<Pos<f32>> {
        Some(Pos::new(self.x.to_f32()?, self.y.to_f32()?))
    }
}

impl Pos<f32> {
    /// Rounds the position's coordinates to their nearest integer value.
    pub fn round(self) -> Pos<f32> {
        Pos::new(self.x.round(), self.y.round())
    }

    /// Converts the position's coordinate type to [isize] (by truncating).
    ///
    /// Fails under the same conditions that make [num::ToPrimitive](ToPrimitive) fail.
    pub fn to_isize(self) -> Option<Pos<isize>> {
        Some(Pos::new(self.x.to_isize()?, self.y.to_isize()?))
    }

    /// Converts the position's coordinate type to [usize] (by truncating).
    ///
    /// Fails under the same conditions that make [num::ToPrimitive](ToPrimitive) fail.
    pub fn to_usize(self) -> Option<Pos<usize>> {
        Some(Pos::new(self.x.to_usize()?, self.y.to_usize()?))
    }
}

impl Pos<isize> {
    /// Converts the position's coordinate type to [f32].
    ///
    /// Fails under the same conditions that make [num::ToPrimitive](ToPrimitive) fail.
    pub fn to_f32(self) -> Option<Pos<f32>> {
        Some(Pos::new(self.x.to_f32()?, self.y.to_f32()?))
    }

    /// Converts the position's coordinate type to [usize].
    ///
    /// Fails under the same conditions that make [num::ToPrimitive](ToPrimitive) fail.
    pub fn to_usize(self) -> Option<Pos<usize>> {
        Some(Pos::new(self.x.to_usize()?, self.y.to_usize()?))
    }
}

impl<T> From<(T, T)> for Pos<T> {
    fn from(value: (T, T)) -> Self {
        Pos::<T>::new(value.0, value.1)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_col_major() {
        let pos = Pos::new(10, 10);
        assert_eq!(pos.col_major(10), 110);
    }
}