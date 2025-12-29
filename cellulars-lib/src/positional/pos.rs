//! Contains logic associated with [Pos].

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
}

impl Pos<usize> {
    pub(crate) fn pack_u32(self) -> u32 {
        ((self.x as u32) << 16) | self.y as u32
    }

    /// Compresses the position to 1D using column-major ordering.
    pub fn col_major(self, height: usize) -> usize {
        self.x * height + self.y
    }

    /// Casts the position's coordinate type to [isize].
    pub fn to_isize(self) -> Pos<isize> {
        Pos::new(self.x as isize, self.y as isize)
    }

    /// Casts the position's coordinate type to [f32].
    pub fn to_f32(self) -> Pos<f32> {
        Pos::new(self.x as f32, self.y as f32)
    }
}

impl Pos<f32> {
    /// Rounds the position's coordinates to their nearest integer value.
    pub fn round(self) -> Pos<f32> {
        Pos::new(self.x.round(), self.y.round())
    }

    /// Casts the position's coordinate type to [isize] (by truncating).
    pub fn to_isize(self) -> Pos<isize> {
        Pos::new(self.x as isize, self.y as isize)
    }

    /// Casts the position's coordinate type to [usize] (by truncating).
    pub fn to_usize(self) -> Pos<usize> {
        Pos::new(self.x as usize, self.y as usize)
    }
}

impl Pos<isize> {
    /// Casts the position's coordinate type to [f32].
    pub fn to_f32(self) -> Pos<f32> {
        Pos::new(self.x as f32, self.y as f32)
    }

    /// Casts the position's coordinate type to [usize].
    pub fn to_usize(self) -> Pos<usize> {
        Pos::new(self.x as usize, self.y as usize)
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