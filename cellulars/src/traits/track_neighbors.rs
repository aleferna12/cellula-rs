use crate::prelude::Spin;

/// Objects that can track cell neighborhoods (usually through [`NeighborTracker`]).
pub trait TrackNeighbors {
    fn neighbor_contacts(&self, spin1: Spin, spin2: Spin) -> Option<u32>;
}