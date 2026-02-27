/// Defines how to write data to a file.
pub trait Write<D, E> {
    /// Writes `data`.
    fn write(self, data: &D) -> Result<(), E>;
}