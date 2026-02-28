//! Contains logic associated with the [`Read`] trait.

/// Defines how to read data from a file.
pub trait Read<D, E> {
    /// Reads this type from `file`.
    fn read(self) -> Result<D, E>;
}