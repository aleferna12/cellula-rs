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
    
    pub fn iter_pairs(&self, start: Spin, end: Spin) -> impl Iterator<Item = (Spin, Spin)> {
        (start..end).flat_map(move |i| (i..end).map(move |j| (i, j)))
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::constants::Spin;

    #[test]
    fn test_index_access() {
        let mut table: SpinTable<u32> = SpinTable::new(4);
        table[(1, 3)] = 42;

        // Check symmetry
        assert_eq!(table[(1, 3)], 42);
        assert_eq!(table[(3, 1)], 42);

        // Overwrite via reversed index
        table[(3, 1)] = 99;
        assert_eq!(table[(1, 3)], 99);
    }

    #[test]
    fn test_iter_pairs_produces_correct_pairs() {
        let table: SpinTable<u8> = SpinTable::new(4);
        let pairs: Vec<(Spin, Spin)> = table.iter_pairs(1, 4).collect();
        let expected = vec![
            (1, 1), (1, 2), (1, 3),
            (2, 2), (2, 3),
            (3, 3),
        ];
        assert_eq!(pairs, expected);
    }
}