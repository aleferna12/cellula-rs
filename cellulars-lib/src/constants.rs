/// Type of cell's spins (determines maximum number of cells allowed in the simulation).
///
/// May also require changing the `discriminant()` method of `LatticeEntity<()>`.
pub type CellIndex = u32;