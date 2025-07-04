use crate::boundary::UnsafePeriodicBoundary;

/// Boundary type of the lattice.
///
/// `FixedBoundary` is approx. 16% faster than `UnsafePeriodicBoundary` (in total run time).
pub type LatticeBoundaryType = UnsafePeriodicBoundary<isize>;

/// Type of cell's spins (may also require changing the spin function of `LatticeEntity`).
pub type Spin = u32;