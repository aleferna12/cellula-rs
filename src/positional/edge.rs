use std::hash::{Hash, Hasher};
use std::mem;
use crate::positional::pos::Pos2D;

#[derive(Debug)]
pub enum EdgeError {
    SamePosition,
    NotNeighbours
}

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