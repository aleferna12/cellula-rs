use cellulars_lib::positional::neighbourhood::MooreNeighbourhood;
#[cfg(not(feature = "fixed_boundary"))]
use cellulars_lib::positional::boundaries::UnsafePeriodicBoundary;
#[cfg(feature = "fixed_boundary")]
use cellulars_lib::positional::boundaries::FixedBoundary;

/// Boundary type of the environment.
///
/// `FixedBoundary` is ~18% faster than `UnsafePeriodicBoundary` (in total run time).
#[cfg(not(feature = "fixed_boundary"))]
pub type BoundaryType = UnsafePeriodicBoundary<f32>;
#[cfg(feature = "fixed_boundary")]
pub type BoundaryType = FixedBoundary<f32>;

/// Neighbourhood type of the environment.
pub type NeighbourhoodType = MooreNeighbourhood;

/// Small value distinguishable from 0.
///
/// Used to compute cell division axis for example.
pub const EPSILON: f32 = 1e-6;