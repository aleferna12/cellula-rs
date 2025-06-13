use std::ops::{Mul, Sub};

const MAX_NEIGH_R: u8 = 16;
const NEIGHBOURHOOD_SIZE: usize = 4 * MAX_NEIGH_R as usize * (MAX_NEIGH_R as usize + 1);
const MOORE_NEIGHS: [(i16, i16); NEIGHBOURHOOD_SIZE] = {
    let mut ret = [(0i16, 0i16); NEIGHBOURHOOD_SIZE];
    let mut r = 1;
    let mut flat_index = 0usize;
    while r <= MAX_NEIGH_R as i16 {
        let mut i = -r;
        while i <= r {
            let mut j = -r;
            while j <= r {
                let max_abs = if i.abs() > j.abs() { i.abs() } else { j.abs() };
                if max_abs == r {
                    ret[flat_index] = (i, j);
                    flat_index += 1;
                }
                j += 1;
            }
            i += 1;
        }
        r += 1;
    }
    ret
};

#[derive(Debug)]
pub enum EdgeError {
    SamePosition,
    NotNeighbours
}

/// 2D position in space.
#[derive(PartialEq, Eq, Copy, Clone, Debug)]
pub struct Pos2D<T> {
    pub x: T,
    pub y: T
}

impl<T> Pos2D<T> {
    pub fn new(x: T, y: T) -> Self {
        Self { x, y }
    }
}

impl Pos2D<usize> {
    pub(crate) fn pack_u32(&self) -> u32 {
        ((self.x as u32) << 16) | self.y as u32
    }

    pub fn row_major(&self, height: usize) -> usize {
        self.x * height + self.y
    }

    pub fn moore_neighs(&self, neigh_r: u8) -> impl Iterator<Item = Pos2D<usize>> {
        let vec_size = 4 * neigh_r * (neigh_r + 1);
        MOORE_NEIGHS[..vec_size as usize]
            .iter()
            .map(|(i, j)| {
                Pos2D::<usize>::new(
                    (self.x as i16 + i) as usize,
                    (self.y as i16 + j) as usize,
                )
            })
    }
}

impl<T> From<(T, T)> for Pos2D<T> {
    fn from(value: (T, T)) -> Self {
        Pos2D::<T>::new(value.0, value.1)
    }
}

#[derive(Clone)]
pub struct Rect<T> {
    pub min: Pos2D<T>,
    pub max: Pos2D<T>
}
impl<T> Rect<T>
where
    T: Sub<Output = T>
        + Mul<Output = T>
        + PartialOrd
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

    pub fn inbounds(&self, pos: Pos2D<T>) -> bool {
        (self.min.x..self.max.x).contains(&pos.x) && (self.min.y..self.max.y).contains(&pos.y)
    }
}

impl Rect<usize> {
    pub fn iterate_pos(&self) -> RectAreaIt {
        RectAreaIt::new(self)
    }
}

pub struct RectAreaIt<'a> {
    curr: Pos2D<usize>,
    rect: &'a Rect<usize>
}
impl<'a> RectAreaIt<'a> {
    fn new(rect: &'a Rect<usize>) -> Self {
        Self {
            curr: rect.min,
            rect
        }
    }
}

impl Iterator for RectAreaIt<'_> {
    type Item = Pos2D<usize>;

    fn next(&mut self) -> Option<Self::Item> {
        if self.curr.y >= self.rect.max.y {
            return None;
        }
        let ret_pos = self.curr;
        if self.curr.x < self.rect.max.x - 1 {
            self.curr.x += 1;
        } else {
            self.curr.x = self.rect.min.x;
            self.curr.y += 1;
        }
        Some(ret_pos)
    }
}

#[cfg(test)]
mod tests {
    use crate::pos::{Rect, MOORE_NEIGHS};

    #[test]
    fn test_rect_area() {
        let r = Rect::<usize>::new((0, 0).into(), (10, 10).into());
        let v: Vec<_> = r.iterate_pos().collect();
        assert_eq!(r.area(), v.len())
    }
    
    #[test]
    fn test_moore() {
        let first_8 = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)];
        assert_eq!(first_8, MOORE_NEIGHS[..8]);
    }
}