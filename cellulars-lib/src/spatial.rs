use crate::constants::Spin;
use crate::lattice::Lattice;
use crate::positional::boundary::AsLatticeBoundary;

pub trait Spatial {
    type Boundary: AsLatticeBoundary;
    fn cell_lattice(&self) -> &Lattice<Spin>;
    fn cell_lattice_mut(&mut self) -> &mut Lattice<Spin>;
    fn boundary(&self) -> &Self::Boundary;
    fn lattice_boundary(&self) -> &<<Self as Spatial>::Boundary as AsLatticeBoundary>::LatticeBoundary;
}