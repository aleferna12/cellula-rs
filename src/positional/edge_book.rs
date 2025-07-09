use std::ops::Index;
use indexmap::IndexSet;
use rand::Rng;
use crate::positional::edge::Edge;

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