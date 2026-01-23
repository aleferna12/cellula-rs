//! Bundle of IO-related submodules.

pub mod io_manager;
pub mod plot;
pub mod parameters;
#[cfg(feature = "movie")]
pub mod movie_maker;
pub mod kinect_listener;