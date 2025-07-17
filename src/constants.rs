use crate::positional::boundary::{FixedBoundary, UnsafePeriodicBoundary};

/// Boundary type of the environment.
///
/// `FixedBoundary` is ~18% faster than `UnsafePeriodicBoundary` (in total run time).
pub type BoundaryType = UnsafePeriodicBoundary<f32>;

/// Type of cell's spins (determines maximum number of cells allowed in the simulation).
/// 
/// May also require changing the spin function of `LatticeEntity`.
pub type Spin = u32;

// Type definition based on BoundaryType
pub type LatticeBoundaryType = <BoundaryType as LatticeBoundaryAssociate>::BoundaryType;

pub trait LatticeBoundaryAssociate {
    type BoundaryType;
}

impl LatticeBoundaryAssociate for UnsafePeriodicBoundary<f32> {
    type BoundaryType = UnsafePeriodicBoundary<isize>;
}

impl LatticeBoundaryAssociate for FixedBoundary<f32> {
    type BoundaryType = FixedBoundary<isize>;
}