use std::hash::{Hash, Hasher};
use std::mem;
use std::ops::Index;
use indexmap::IndexSet;
use rand::Rng;
use crate::pos::{EdgeError, Pos2D};

#[derive(Eq, Clone)]
pub struct Edge {
    pub p1: Pos2D<usize>,
    pub p2: Pos2D<usize>
}

impl Edge {
    pub fn new(p1: Pos2D<usize>, p2: Pos2D<usize>) -> Self {
        Self { p1, p2 }
    }
    
    pub fn new_if_neighbour(p1: Pos2D<usize>, p2: Pos2D<usize>, neigh_r: u8) -> Result<Self, EdgeError> {
        let cx = p1.x.abs_diff(p2.x);
        let cy = p1.y.abs_diff(p2.y);
        let sum = cx + cy;
        if sum == 0 {
            return Err(EdgeError::SamePosition);
        }
        if sum > (neigh_r * 2) as usize {
            return Err(EdgeError::NotNeighbours);
        }
        Ok(Self { p1, p2 })
    }

    fn hash_u64(&self) -> u64 {
        let mut u1 = self.p1.pack_u32();
        let mut u2 = self.p2.pack_u32();
        if u1 > u2 {
            mem::swap(&mut u1, &mut u2);
        }
        ((u1 as u64) << 32) | (u2 as u64)
    }
}

impl PartialEq for Edge {
    fn eq(&self, other: &Self) -> bool {
        self.hash_u64() == other.hash_u64()
    }
}

// TODO: test the perfect hash algorithm that steven sent me
//       also test a perfect hash where we use row-major ordering to index a vector of all possible edges
//       the vector has size 9 * width * height (since each position has 8 neighbours), and the hash function
//       is: row_major(pos1, height) + row_major(pos2, 3)
//       this would mean replacing edge_set with a Vec
impl Hash for Edge {
    fn hash<H: Hasher>(&self, state: &mut H) {
        self.hash_u64().hash(state);
    }
}

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