use std::cmp::min;
use std::ops::Index;
use indexmap::IndexSet;
use rand::random_range;
use rustc_hash::FxBuildHasher;
use crate::pos::{Edge, Pos2D};

pub struct Lattice<T> {
    pub width: usize,
    pub height: usize,
    array: Box<[T]>,
    edge_set: IndexSet<Edge, FxBuildHasher>
}

impl<T: Default + Clone> Lattice<T> {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            array: vec![T::default(); width * height].into_boxed_slice(),
            edge_set: IndexSet::with_hasher(FxBuildHasher::default())
        }
    }

    // TODO: this should take an rng
    pub fn random_pos(&self) -> Pos2D<usize> {
        Pos2D::new(random_range(0..self.width - 1), random_range(0..self.height - 1))
    }
    
    pub fn n_edges(&self) -> usize { self.edge_set.len() }
    
    pub fn insert_edge(&mut self, edge: Edge) -> bool {
        self.edge_set.insert(edge)
    }

    // TODO: this should take an rng
    pub fn remove_random_edge(&mut self) -> Edge {
        let index = random_range(0..self.edge_set.len() - 1);
        self.edge_set.swap_remove_index(index).unwrap()
    }

    // TODO: this should take an rng
    pub fn random_neighbour(&self, p: &Pos2D<usize>) -> Pos2D<usize> {
        let oldp = (p.x as i32, p.y as i32);
        let mut newp = oldp;
        while oldp == newp {
            newp.0 = oldp.0 + random_range(-min(1, oldp.0)..min(2, self.width as i32 - oldp.0));
            newp.1 = oldp.1 + random_range(-min(1, oldp.1)..min(2, self.height as i32 - oldp.1));
        }
        Pos2D::new(newp.0 as usize, newp.1 as usize)
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

mod tests {
    use super::*;
    
    #[test]
    fn test_random_neighbour() {
        let lat = Lattice::<u32>::new(20, 20);
        for _ in 0..1000 {
            let p1 = lat.random_pos();
            let p2 = lat.random_neighbour(&p1);
            assert!(Edge::new(p1, p2).is_ok())
        }
    }
}