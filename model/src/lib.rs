//! This crate contains code that can be directly modified to extend the base model implementation provided.
//!
//! The base [`model::Model`] bundles a comprehensive set of IO-related features and showcases how the [`cellulars_lib`]
//! library can be extended by implementing cell chemotaxis and cell division.

pub mod model;
pub mod io;
pub mod my_pond;
pub mod potts;
pub mod constants;
pub mod my_cell;
pub mod my_environment;
pub mod kinect;