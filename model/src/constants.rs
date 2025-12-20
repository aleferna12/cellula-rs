//! Contains constants that are set at compile-time with feature flags.

#[cfg(not(feature = "von-neumann"))]
use cellulars_lib::positional::neighbourhood::MooreNeighbourhood;
#[cfg(feature = "von-neumann")]
use cellulars_lib::positional::neighbourhood::VonNeumannNeighbourhood;
#[cfg(not(feature = "fixed-boundary"))]
use cellulars_lib::positional::boundaries::UnsafePeriodicBoundary;
#[cfg(feature = "fixed-boundary")]
use cellulars_lib::positional::boundaries::FixedBoundary;

/// Boundary type of the environment.
///
/// [FixedBoundary] is ~18% faster than [UnsafePeriodicBoundary]
/// (in total run time).
#[cfg(not(feature = "fixed-boundary"))]
pub type BoundaryType = UnsafePeriodicBoundary<f32>;
#[cfg(feature = "fixed-boundary")]
pub type BoundaryType = FixedBoundary<f32>;

/// Neighbourhood type of the environment.
#[cfg(not(feature = "von-neumann"))]
pub type NeighbourhoodType = MooreNeighbourhood;
#[cfg(feature = "von-neumann")]
pub type NeighbourhoodType = VonNeumannNeighbourhood;

/// Small value distinguishable from 0.
///
/// Used to compute cell division axis for example.
pub const EPSILON: f32 = 1e-6;