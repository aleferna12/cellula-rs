//! Contains logic associated with [Edge].

use crate::positional::pos::Pos;
use std::hash::{Hash, Hasher};
use std::mem;

/// A symmetrical edge between two positions (such that `Edge(p1, p2) == Edge(p2, p1)`).
#[derive(Eq, Clone)]
#[derive(Debug)]
pub struct Edge {
    /// Position at one end of the edge.
    pub p1: Pos<usize>,
    /// Position at the other end of the edge.
    pub p2: Pos<usize>
}

impl Edge {
    /// Makes a new edge between `pos1` and `pos2` without checking if they are neighbours.
    pub fn new(p1: Pos<usize>, p2: Pos<usize>) -> Self {
        Self { p1, p2 }
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