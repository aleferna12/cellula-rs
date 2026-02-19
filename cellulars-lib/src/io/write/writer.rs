//! Contains [`Writer`].

use std::path::Path;

// TODO: implement other file formats
/// Default writer.
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Writer {}

/// Defines how to write data to a file.
pub trait Write<D, E> {
    /// Writes `data` to `path`.
    fn write(&mut self, data: &D, path: impl AsRef<Path>) -> Result<(), E>;
}
