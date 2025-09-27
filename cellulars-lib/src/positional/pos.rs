/// 2D position in space.
#[derive(PartialEq, Eq, Copy, Clone, Debug)]
#[derive(Hash)]
pub struct Pos<T> {
    pub x: T,
    pub y: T
}

impl<T> Pos<T> {
    pub fn new(x: T, y: T) -> Self {
        Self { x, y }
    }
}

impl Pos<usize> {
    pub(crate) fn pack_u32(self) -> u32 {
        ((self.x as u32) << 16) | self.y as u32
    }

    pub fn row_major(self, width: usize) -> usize {
        self.x * width + self.y
    }

    pub fn to_isize(self) -> Pos<isize> {
        Pos::new(self.x as isize, self.y as isize)
    }

    pub fn to_f32(self) -> Pos<f32> {
        Pos::new(self.x as f32, self.y as f32)
    }
}

impl Pos<f32> {
    pub fn to_isize(self) -> Pos<isize> {
        Pos::new(self.x as isize, self.y as isize)
    }

    pub fn to_usize(self) -> Pos<usize> {
        Pos::new(self.x as usize, self.y as usize)
    }
}

impl Pos<isize> {
    pub fn to_f32(self) -> Pos<f32> {
        Pos::new(self.x as f32, self.y as f32)
    }

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
    fn test_row_major() {
        let pos = Pos::new(10, 10);
        assert_eq!(pos.row_major(10), 110);
    }
}