use crate::constants::Spin;
use std::ops::{Index, IndexMut};

pub struct SpinTable<T> {
    array: Box<[T]>,
    max_spin: Spin,
}

impl<T: Default + Clone> SpinTable<T> {
    pub fn new(max_spin: Spin) -> Self {
        let size = max_spin * (max_spin + 1) / 2;
        Self {
            array: vec![T::default(); size as usize].into_boxed_slice(),
            max_spin
        }
    }
    
    pub fn iter_pairs(&self) -> impl Iterator<Item = (Spin, Spin)> {
        (0..=self.max_spin).flat_map(|i| (0..=self.max_spin).map(move |j| (i, j)))
    }

    fn flat_index(&self, i: Spin, j: Spin) -> usize {
        let (i, j) = if i > j { (j, i) } else { (i, j) };
        let ind = i * (2 * self.max_spin - i - 1) / 2 + j - i;
        ind as usize
    }
}

impl<T: Default + Clone> Index<(Spin, Spin)> for SpinTable<T> {
    type Output = T;

    fn index(&self, index: (Spin, Spin)) -> &Self::Output {
        &self.array[self.flat_index(index.0, index.1)]
    }
}

impl<T: Default + Clone> IndexMut<(Spin, Spin)> for SpinTable<T> {
    fn index_mut(&mut self, index: (Spin, Spin)) -> &mut Self::Output {
        &mut self.array[self.flat_index(index.0, index.1)]
    }
}