use std::path::Path;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Writer {}

pub trait Write<D, E> {
    fn write(&mut self, data: &D, path: impl AsRef<Path>) -> Result<(), E>;
}
