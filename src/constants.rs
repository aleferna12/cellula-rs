use crate::positional::boundary::UnsafePeriodicBoundary;

/// Boundary type of the lattice.
///
/// `FixedBoundary` is ~18% faster than `UnsafePeriodicBoundary` (in total run time).
pub type LatticeBoundaryType = UnsafePeriodicBoundary<isize>;

/// Type of cell's spins (determines maximum number of cells allowed in the simulation).
/// 
/// May also require changing the spin function of `LatticeEntity`.
pub type Spin = u32;