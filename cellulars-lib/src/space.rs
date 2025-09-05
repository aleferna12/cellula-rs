use crate::positional::boundary::ToLatticeBoundary;
use std::error::Error;

pub trait Spatial {
    type Boundary: ToLatticeBoundary;
    fn boundary(&self) -> &Self::Boundary;
    fn lattice_boundary(&self) -> &<<Self as Spatial>::Boundary as ToLatticeBoundary>::LatticeBoundary;
}

pub struct Space<B: ToLatticeBoundary> {
    bound: B,
    lat_bound: B::LatticeBoundary,
}

impl<B: ToLatticeBoundary> Space<B> {
    pub fn new(bound: B) -> Result<Self, B::Error>
    where
        B: ToLatticeBoundary<Coord = f32>,
        B::Error: 'static + Error {
        Ok(Self {
            lat_bound: bound.to_lattice_boundary()?,
            bound,
        })
    }
}

impl<B: ToLatticeBoundary> Spatial for Space<B> {
    type Boundary = B;

    fn boundary(&self) -> &Self::Boundary {
        &self.bound
    }

    fn lattice_boundary(&self) -> &<<Self as Spatial>::Boundary as ToLatticeBoundary>::LatticeBoundary {
        &self.lat_bound
    }
}