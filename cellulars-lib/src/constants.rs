//! Contains constants.

/// Type of cell's indexes (also determines maximum number of cells allowed in the simulation).
pub type CellIndex = u32;

/// Type of the floats used throughout the simulation.
/// 
/// Using [f32] instead if [f64] provides some performance gain at the cost of mathematical precision.
#[cfg(feature = "high-precision")]
pub type FloatType = f64;
#[cfg(not(feature = "high-precision"))]
pub type FloatType = f32;