//! This crate contains code that can be directly modified to extend the base model implementation provided.
//!
//! The base [`model::Model`] bundles a comprehensive set of IO-related features and showcases how the [`cellulars`]
//! library can be extended by implementing cell chemotaxis and cell division.

pub mod model;
pub mod io;
pub mod constants;
pub mod my_cell;
pub mod my_environment;
pub mod pond;
pub mod chemotaxis_bias;