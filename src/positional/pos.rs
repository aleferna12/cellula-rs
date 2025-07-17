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

impl<T> From<(T, T)> for Pos<T> {
    fn from(value: (T, T)) -> Self {
        Pos::<T>::new(value.0, value.1)
    }
}

impl Pos<usize> {
    pub(crate) fn pack_u32(self) -> u32 {
        ((self.x as u32) << 16) | self.y as u32
    }

    pub fn row_major(self, height: usize) -> usize {
        self.x * height + self.y
    }
}

impl From<Pos<usize>> for Pos<isize> {
    fn from(value: Pos<usize>) -> Self {
        Pos::new(value.x as isize, value.y as isize)
    }
}

impl From<Pos<isize>> for Pos<usize> {
    fn from(value: Pos<isize>) -> Self {
        let message = "overflow when translating position from general to lattice coordinates";
        Pos::new(
            value.x.try_into().expect(message), 
            value.y.try_into().expect(message)
        )
    }
}