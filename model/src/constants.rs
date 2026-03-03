//! Contains constants that are set at compile-time with feature flags.

use cellulars::constants::FloatType;
#[cfg(feature = "fixed-boundary")]
use cellulars::positional::boundaries::FixedBoundary;
#[cfg(not(feature = "fixed-boundary"))]
use cellulars::positional::boundaries::UnsafePeriodicBoundary;
#[cfg(not(feature = "von-neumann"))]
use cellulars::positional::neighborhood::MooreNeighborhood;
#[cfg(feature = "von-neumann")]
use cellulars::positional::neighborhood::VonNeumannNeighborhood;

/// Boundary type of the environment.
///
/// [`FixedBoundary`](cellulars::positional::boundaries::FixedBoundary) is ~18% faster than [`UnsafePeriodicBoundary`]
/// (in total run time).
#[cfg(not(feature = "fixed-boundary"))]
pub type BoundaryType = UnsafePeriodicBoundary<FloatType>;
#[cfg(feature = "fixed-boundary")]
pub type BoundaryType = FixedBoundary<FloatType>;

/// Neighborhood type of the environment.
#[cfg(not(feature = "von-neumann"))]
pub type NeighborhoodType = MooreNeighborhood;
#[cfg(feature = "von-neumann")]
pub type NeighborhoodType = VonNeumannNeighborhood;

/// Small value distinguishable from 0.
///
/// Used to compute cell division axis for example.
pub const EPSILON: FloatType = 1e-6;