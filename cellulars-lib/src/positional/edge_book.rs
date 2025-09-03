use crate::positional::edge::Edge;
use indexmap::IndexSet;
use rand::Rng;
use std::ops::Index;

/// This struct exists to book-keep the edges in a lattice.
pub struct EdgeBook {
    // TODO: profile using this crate, I have no clue of whether it's fast enough
    edge_set: IndexSet<Edge>,
}

impl EdgeBook {
    pub fn new() -> Self {
        Self { edge_set: IndexSet::new() }
    }

    pub fn len(&self) -> usize { self.edge_set.len() }

    pub fn is_empty(&self) -> bool { self.len() == 0 }

    pub fn at(&self, index: usize) -> &Edge {
        self.edge_set.index(index)
    }

    pub fn remove_at(&mut self, index: usize) -> Option<Edge> {
        self.edge_set.swap_remove_index(index)
    }

    pub fn insert(&mut self, edge: Edge) -> bool {
        self.edge_set.insert(edge)
    }

    pub fn remove(&mut self, edge: &Edge) -> bool {
        self.edge_set.swap_remove(edge)
    }

    pub fn random_index(&self, rng: &mut impl Rng) -> usize {
        rng.random_range(0..self.edge_set.len())
    }
}

impl Default for EdgeBook {
    fn default() -> Self {
        EdgeBook::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::positional::edge::Edge;
    use crate::positional::pos::Pos;
    use rand::{SeedableRng, rngs::StdRng};
    use std::collections::HashSet;

    fn make_edge(a: (usize, usize), b: (usize, usize)) -> Edge {
        Edge::new(Pos::from(a), Pos::from(b))
    }

    #[test]
    fn test_insert_unique_and_duplicate() {
        let mut book = EdgeBook::new();
        let edge = make_edge((0, 0), (0, 1));

        assert!(book.insert(edge.clone()));
        assert_eq!(book.len(), 1);

        // Insert same edge again: should not increase len
        assert!(!book.insert(edge));
        assert_eq!(book.len(), 1);
    }

    #[test]
    fn test_remove_by_value_and_index() {
        let mut book = EdgeBook::new();
        let e1 = make_edge((0, 0), (0, 1));
        let e2 = make_edge((1, 1), (1, 2));
        book.insert(e1.clone());
        book.insert(e2.clone());

        assert!(book.remove(&e1));
        assert_eq!(book.len(), 1);
        assert!(!book.remove(&e1)); // Already removed

        // remove_at
        let removed = book.remove_at(0);
        assert_eq!(removed.unwrap(), e2);
        assert_eq!(book.len(), 0);
    }

    #[test]
    fn test_random_index_within_bounds() {
        let mut rng = StdRng::seed_from_u64(1234);
        let mut book = EdgeBook::new();

        for i in 0..10 {
            book.insert(make_edge((i, 0), (i, 1)));
        }

        for _ in 0..100 {
            let idx = book.random_index(&mut rng);
            assert!(idx < book.len());
        }
    }

    #[test]
    fn test_preserves_uniqueness() {
        let mut book = EdgeBook::new();
        let mut inserted = HashSet::new();
        for i in 0..5 {
            let e = make_edge((i, 0), (i, 1));
            assert!(book.insert(e.clone()));
            inserted.insert(e);
        }

        for e in &inserted {
            assert!(book.remove(e));
        }

        assert!(book.is_empty());
    }
}
