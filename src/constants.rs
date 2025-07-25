use crate::positional::boundary::UnsafePeriodicBoundary;
use crate::positional::neighbourhood::MooreNeighbourhood;

/// Boundary type of the environment.
///
/// `FixedBoundary` is ~18% faster than `UnsafePeriodicBoundary` (in total run time).
pub type BoundaryType = UnsafePeriodicBoundary<f32>;

/// Neighbourhood type of the environment.
pub type NeighbourhoodType = MooreNeighbourhood;

/// Type of cell's spins (determines maximum number of cells allowed in the simulation).
/// 
/// May also require changing the `discriminant()` method of `LatticeEntity<()>`.
pub type Spin = u32;