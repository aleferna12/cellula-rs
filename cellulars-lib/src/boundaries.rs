use std::error::Error;
use crate::positional::boundary::{PosValidator, ToLatticeBoundary};

pub struct boundaries<B: ToLatticeBoundary> {
    pub boundary: B,
    pub lattice_boundary: B::LatticeBoundary,
}

impl<B: ToLatticeBoundary> boundaries<B> {
    pub fn new(bound: B) -> Result<Self, B::Error>
    where
        B: ToLatticeBoundary<Coord = f32>,
        B::Error: 'static + Error {
        Ok(Self {
            lattice_boundary: bound.to_lattice_boundary()?,
            boundary: bound,
        })
    }
}

impl<B: ToLatticeBoundary<Coord = f32>> PosValidator for boundaries<B> {
    type Boundary = B;

    fn boundary(&self) -> &Self::Boundary {
        &self.boundary
    }

    fn lattice_boundary(&self) -> &<Self::Boundary as ToLatticeBoundary>::LatticeBoundary {
        &self.lattice_boundary
    }
}