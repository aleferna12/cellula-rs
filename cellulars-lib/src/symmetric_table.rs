//! Contains logic associated with [SymmetricTable].

use std::ops::{Index, IndexMut};

/// A symmetric table where indexes (x, y) and (y, x) map to the same value.
#[derive(Clone, Debug, PartialEq)]
pub struct SymmetricTable<T> {
    array: Box<[T]>,
    length: usize,
}

impl<T> SymmetricTable<T> {
    /// Returns the length of the sides of the table.
    pub fn length(&self) -> usize {
        self.length
    }

    /// Iterates over all unique pairs of indexes that can be used for the table.
    pub fn iter_index_pairs(&self, start: Option<usize>, end: Option<usize>) -> impl Iterator<Item = (usize, usize)> + use<T> {
        let start = start.unwrap_or(0);
        let end = end.unwrap_or(self.length);
        (start..end).flat_map(move |i| (i..end).map(move |j| (i, j)))
    }

    fn flat_index(&self, i: usize, j: usize) -> usize {
        assert!(i < self.length && j < self.length, "index out of bounds");
        let (i, j) = if i > j { (j, i) } else { (i, j) };
        // The upper triangle is stored row by row, diagonal included, so row `i` starts at
        // `i * length - i * (i - 1) / 2` (the full rows above it, minus the triangle they skip)
        i * (2 * self.length - i + 1) / 2 + j - i
    }
}

impl<T: Default + Clone> SymmetricTable<T> {
    /// Makes a new `length`x`length` table
    pub fn new(length: usize) -> Self {
        let size = length * (length + 1) / 2;
        Self {
            array: vec![T::default(); size].into_boxed_slice(),
            length
        }
    }

    /// Clears the table by setting all values to the default of the table's inner type.
    pub fn clear(&mut self) {
        self.array.fill(T::default());
    }
}

impl<T> Index<(usize, usize)> for SymmetricTable<T> {
    type Output = T;

    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.array[self.flat_index(index.0, index.1)]
    }
}

impl<T> IndexMut<(usize, usize)> for SymmetricTable<T> {
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        &mut self.array[self.flat_index(index.0, index.1)]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_index_access() {
        let mut table: SymmetricTable<u32> = SymmetricTable::new(4);
        table[(1, 3)] = 42;

        // Check symmetry
        assert_eq!(table[(1, 3)], 42);
        assert_eq!(table[(3, 1)], 42);

        // Overwrite via reversed index
        table[(3, 1)] = 99;
        assert_eq!(table[(1, 3)], 99);
    }

    #[test]
    fn test_every_pair_maps_to_its_own_slot() {
        let length = 5;
        let mut table: SymmetricTable<u32> = SymmetricTable::new(length);
        let pairs: Vec<_> = table.iter_index_pairs(None, None).collect();
        // Writing a unique value per pair must not disturb any of the other pairs
        for (n, &(i, j)) in pairs.iter().enumerate() {
            table[(i, j)] = n as u32 + 1;
        }
        for (n, &(i, j)) in pairs.iter().enumerate() {
            assert_eq!(table[(i, j)], n as u32 + 1, "index ({i}, {j}) aliases another pair");
            assert_eq!(table[(j, i)], n as u32 + 1, "index ({j}, {i}) aliases another pair");
        }
    }

    #[test]
    fn test_flat_indexes_fill_the_array() {
        let length = 6;
        let table: SymmetricTable<u8> = SymmetricTable::new(length);
        let mut flat: Vec<_> = table
            .iter_index_pairs(None, None)
            .map(|(i, j)| table.flat_index(i, j))
            .collect();
        flat.sort_unstable();
        // Every slot of the stored upper triangle is used exactly once
        assert_eq!(flat, (0..length * (length + 1) / 2).collect::<Vec<_>>());
    }

    #[test]
    #[should_panic(expected = "index out of bounds")]
    fn test_index_past_the_side_length_panics() {
        let table: SymmetricTable<u8> = SymmetricTable::new(4);
        // Without the bounds check this maps onto a slot that belongs to a valid pair
        let _ = table[(0, 4)];
    }

    #[test]
    fn test_iter_pairs_produces_correct_pairs() {
        let table: SymmetricTable<u8> = SymmetricTable::new(4);
        let pairs: Vec<(usize, usize)> = table.iter_index_pairs(Some(1), Some(4)).collect();
        let expected = vec![
            (1, 1), (1, 2), (1, 3),
            (2, 2), (2, 3),
            (3, 3),
        ];
        assert_eq!(pairs, expected);
    }
}