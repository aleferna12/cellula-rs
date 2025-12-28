//! This crate contains library functions for the [Cellulars](https://github.com/aleferna12/cellulars) project.

// TODO!:
//  - Implement EnvNeighbours, which tracks neighbours
//  - Implement basic traits like Debug and Clone for all types
pub mod lattice;
pub mod positional;
pub mod environment;
pub mod static_adhesion;
pub mod cell_container;
pub mod symmetric_table;
pub mod base_cell;
pub mod constants;
pub mod spin;
pub mod pond;
pub mod traits;