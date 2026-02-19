//! Bundle of submodules used to write data from a simulation.

pub mod writer;
#[cfg(feature = "image-io")]
pub mod image;
#[cfg(feature = "data-io")]
pub mod data_write;