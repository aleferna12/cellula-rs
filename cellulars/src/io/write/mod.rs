//! Bundle of submodules used to write data from a simulation.

#[cfg(feature = "image-io")]
pub mod image;
#[cfg(feature = "data-io")]
pub mod parquet_writer;
pub mod write_trait;