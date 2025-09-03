use cellulars_lib::environment::Environment;
use cellulars_lib::positional::boundary::UnsafePeriodicBoundary;
use cellulars_lib::positional::neighbourhood::MooreNeighbourhood;
use crate::cell::Cell;

/// Boundary type of the environment.
///
/// `FixedBoundary` is ~18% faster than `UnsafePeriodicBoundary` (in total run time).
pub type BoundaryType = UnsafePeriodicBoundary<f32>;

/// Neighbourhood type of the environment.
pub type NeighbourhoodType = MooreNeighbourhood;