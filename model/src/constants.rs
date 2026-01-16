//! Contains constants that are set at compile-time with feature flags.

use cellulars_lib::constants::FloatType;
#[cfg(feature = "fixed-boundary")]
use cellulars_lib::positional::boundaries::FixedBoundary;
#[cfg(not(feature = "fixed-boundary"))]
use cellulars_lib::positional::boundaries::UnsafePeriodicBoundary;
use cellulars_lib::positional::neighbourhood::MooreNeighbourhood;
#[cfg(feature = "von-neumann")]
use cellulars_lib::positional::neighbourhood::VonNeumannNeighbourhood;

/// Boundary type of the environment.
///
/// [`FixedBoundary`](cellulars_lib::positional::boundaries::FixedBoundary) is ~18% faster than [`UnsafePeriodicBoundary`]
/// (in total run time).
#[cfg(not(feature = "fixed-boundary"))]
pub type BoundaryType = UnsafePeriodicBoundary<FloatType>;
#[cfg(feature = "fixed-boundary")]
pub type BoundaryType = FixedBoundary<FloatType>;

/// Neighbourhood type of the environment.
#[cfg(not(feature = "von-neumann"))]
pub type NeighbourhoodType = MooreNeighbourhood;
#[cfg(feature = "von-neumann")]
pub type NeighbourhoodType = VonNeumannNeighbourhood;

/// Neighbourhood type to filter noise from the kinect.
pub type KinectNeighbourhoodType = MooreNeighbourhood;

/// Small value distinguishable from 0.
///
/// Used to compute cell division axis for example.
pub const EPSILON: FloatType = 1e-6;