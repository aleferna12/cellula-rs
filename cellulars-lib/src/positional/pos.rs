//! Contains logic associated with [`Pos`].

use num::cast::AsPrimitive;

/// A 2D position in space.
#[derive(PartialEq, Eq, Copy, Clone, Debug, Hash)]
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

    /// Casts the position to another coordinate type using the `as` operator.
    pub fn cast_as<U>(self) -> Pos<U>
    where
        T: AsPrimitive<U>,
        U: Copy + 'static {
        Pos::new(
            self.x.as_(),
            self.y.as_(),
        )
    }
}

impl<F: AsPrimitive<T>, T: Copy + 'static> CastCoords<T> for Pos<F> {
    type Outer<U> = Pos<U>;

    fn cast_coords(&self) -> Self::Outer<T> {
        Pos::new(
            self.x.as_(),
            self.y.as_(),
        )
    }
}

impl<T> From<(T, T)> for Pos<T> {
    fn from(value: (T, T)) -> Self {
        Pos::<T>::new(value.0, value.1)
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
}

impl Pos<f32> {
    /// Rounds the position's coordinates to their nearest integer value.
    pub fn round(self) -> Pos<f32> {
        Pos::new(self.x.round(), self.y.round())
    }
}

impl Pos<f64> {
    /// Rounds the position's coordinates to their nearest integer value.
    pub fn round(self) -> Pos<f64> {
        Pos::new(self.x.round(), self.y.round())
    }
}

/// Indicates that the coordinate type of a structure can be converted to `T`. 
pub trait CastCoords<T: Copy + 'static> {
    /// Outer type, which should match the implementor.
    type Outer<U>;

    /// Converts the coordinates system of `self` to type `T`.
    ///
    /// This conversion is performed using the `as` operator, which can lead to narrowing and loss of precision
    /// (see [`num::traits::AsPrimitive`](AsPrimitive)).
    fn cast_coords(&self) -> Self::Outer<T>;
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