//! Bundle of submodules used for image IO.

pub mod lerper;
pub mod plot;
pub mod webp_writer;
#[cfg(feature = "movie-io")]
pub mod movie_window;