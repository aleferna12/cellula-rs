use crate::positional::boundary::UnsafePeriodicBoundary;
use crate::positional::neighbourhood::MooreNeighbourhood;

/// Boundary type of the environment.
///
/// `FixedBoundary` is ~18% faster than `UnsafePeriodicBoundary` (in total run time).
pub type BoundaryType = UnsafePeriodicBoundary<f32>;

/// Neighbourhood type of the environment.
pub type NeighbourhoodType = MooreNeighbourhood;