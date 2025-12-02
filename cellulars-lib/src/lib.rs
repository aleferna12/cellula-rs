//! This crate contains library functions for the [Cellulars](https://github.com/aleferna12/cellulars) project.

// TODO!:
//     - Implement EnvNeighbours, which tracks neighbours
pub mod lattice;
pub mod positional;
pub mod environment;
pub mod adhesion;
pub mod cell_container;
pub mod symmetric_table;
pub mod basic_cell;
pub mod constants;
pub mod spin;
pub mod potts;
pub mod step;
pub mod habitable;