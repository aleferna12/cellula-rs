use std::ops::{Index, IndexMut};
use rand::Rng;
use crate::moore::MOORE_NEIGHS;
use crate::pos::Pos2D;

pub struct Lattice<T> {
    pub width: usize,
    pub height: usize,
    array: Box<[T]>,
}

impl<T: Default + Clone> Lattice<T> {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            array: vec![T::default(); width * height].into_boxed_slice()
        }
    }

    pub fn random_pos(&self, rng: &mut impl Rng) -> Pos2D<usize> {
        Pos2D::new(rng.random_range(0..self.width - 1), rng.random_range(0..self.height - 1))
    }
    
    pub fn inbounds(&self, pos: &Pos2D<usize>) -> bool {
        (0..self.width).contains(&pos.x) && (0..self.height).contains(&pos.y)
    }
    
    // TODO: implement this as an iterator to avoid heap allocation? Dont know if that is faster though
    #[inline]
    pub fn moore_neighs(&self, pos: &Pos2D<usize>, neigh_r: u8) -> Vec<Pos2D<usize>> {
        let vec_size = 4 * neigh_r * (neigh_r + 1);
        MOORE_NEIGHS[..vec_size as usize]
            .iter()
            .filter_map(|(i, j)| {
                let p = Pos2D::<usize>::new(
                    (pos.x as i32 + i) as usize,
                    (pos.y as i32 + j) as usize,
                );
                if self.inbounds(&p) { 
                    return Some(p) 
                }
                None
            })
            .collect()
    }
}

impl<T> Index<Pos2D<usize>> for Lattice<T> {
    type Output = T;

    fn index(&self, pos: Pos2D<usize>) -> &Self::Output {
        &self.array[pos.row_major(self.height)]
    }
}
impl<T> Index<(usize, usize)> for Lattice<T> {
    type Output = T;

    fn index(&self, pos: (usize, usize)) -> &Self::Output {
        &self[Pos2D::<usize>::from(pos)]
    }
}

impl<T> IndexMut<Pos2D<usize>> for Lattice<T> {
    fn index_mut(&mut self, pos: Pos2D<usize>) -> &mut Self::Output {
        &mut self.array[pos.row_major(self.height)]
    }
}
impl<T> IndexMut<(usize, usize)> for Lattice<T> {
    fn index_mut(&mut self, pos: (usize, usize)) -> &mut Self::Output {
        &mut self[Pos2D::<usize>::from(pos)]
    }
}