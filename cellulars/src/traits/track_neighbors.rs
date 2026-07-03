use crate::prelude::Spin;

pub trait TrackNeighbors {
    fn neighbor_contacts(&self, spin1: Spin, spin2: Spin) -> Option<u32>;
}