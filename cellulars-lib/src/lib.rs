//! This crate contains library functions for the [Cellulars](https://github.com/aleferna12/cellulars) project.

/*
TODO!:
 - Implement EnvNeighbours, which tracks neighbours
 - Make a BasePotts, BaseIoManager, BaseModel etc such that we can move the crate fully to cellulars-lib
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
pub mod prelude;