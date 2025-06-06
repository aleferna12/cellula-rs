use std::ops::Index;
use indexmap::IndexSet;
use rand::random_range;
use crate::pos::{Edge, Pos2D};

pub struct Lattice<T> {
    pub width: usize,
    pub height: usize,
    pub array: Box<[T]>,
    edge_set: IndexSet<Edge>
}

impl<T: Default + Clone> Lattice<T> {
    pub fn new(width: usize, height: usize) -> Self {
        Self {
            width,
            height,
            array: vec![T::default(); width * height].into_boxed_slice(),
            edge_set: IndexSet::new()
        }
    }
    
    fn insert_edge(&mut self, edge: Edge) -> bool {
        self.edge_set.insert(edge)
    }

    // TODO: this should take an rng
    fn remove_random_edge(&mut self) -> Option<Edge> {
        let index = random_range(0..self.edge_set.len() - 1);
        self.edge_set.swap_remove_index(index)
    }
}

impl<T> Index<Pos2D<usize>> for Lattice<T> {
    type Output = T;

    fn index(&self, pos: Pos2D<usize>) -> &Self::Output {
        &self.array[pos.x * self.height + pos.y]
    }
}

impl<T> Index<(usize, usize)> for Lattice<T> {
    type Output = T;

    fn index(&self, pos: (usize, usize)) -> &Self::Output {
        &self[Pos2D::<usize>::from(pos)]
    }
}