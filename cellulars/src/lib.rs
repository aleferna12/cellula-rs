//! This crate contains library functions for the [`Cellulars`](https://github.com/aleferna12/cellulars) project.

/*
TODO!:
 - Implement EnvNeighbors, which tracks neighbors
 - Make a BasePotts, BaseIoManager, BaseModel etc such that we can move the crate fully to cellulars-lib
 - Add a command to parse a grayscale PNG into a simulation
 */

pub mod lattice;
pub mod positional;
pub mod static_adhesion;
pub mod cell_container;
pub mod symmetric_table;
pub mod constants;
pub mod spin;
pub mod traits;
pub mod base;
#[cfg(any(feature = "data-io", feature = "image-io"))]
pub mod io;
pub mod prelude;
pub mod empty_cell;