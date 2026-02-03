use std::path::PathBuf;

#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Writer {
    pub outdir: PathBuf
}

pub trait Write<D, E> {
    fn write(&mut self, data: &D, time_step: u32) -> Result<PathBuf, E>;
}
