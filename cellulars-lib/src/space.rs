use crate::constants::Spin;
use crate::lattice::Lattice;
use crate::positional::boundary::AsLatticeBoundary;
use crate::positional::rect::Rect;
use crate::spatial::Spatial;
use std::error::Error;

pub struct Space<B: AsLatticeBoundary> {
    pub bound: B,
    pub lat_bound: B::LatticeBoundary,
    pub cell_lattice: Lattice<Spin>,
}

impl<B: AsLatticeBoundary> Space<B> {
    pub fn new(bound: B) -> Result<Self, Box<dyn Error>>
    where
        B: AsLatticeBoundary<Coord = f32>,
        B::Error: 'static + Error {
        let rect: Rect<usize> = bound.rect().clone().try_into()?;
        Ok(Self {
            lat_bound: bound.as_lattice_boundary()?,
            cell_lattice: Lattice::<Spin>::new(rect.clone()),
            bound,
        })
    }
}

impl<B: AsLatticeBoundary> Spatial for Space<B> {
    type Boundary = B;

    fn cell_lattice(&self) -> &Lattice<Spin> {
        &self.cell_lattice
    }

    fn cell_lattice_mut(&mut self) -> &mut Lattice<Spin> {
        &mut self.cell_lattice
    }

    fn boundary(&self) -> &Self::Boundary {
        &self.bound
    }

    fn lattice_boundary(&self) -> &<<Self as Spatial>::Boundary as AsLatticeBoundary>::LatticeBoundary {
        &self.lat_bound
    }
}