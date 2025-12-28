//! Prelude containing commonly used items.

pub use crate::base::base_cell::BaseCell;
pub use crate::base::base_environment::BaseEnvironment;
pub use crate::base::base_pond::BasePond;
pub use crate::cell_container::{CellContainer, RelCell};
pub use crate::constants::CellIndex;
pub use crate::lattice::Lattice;
pub use crate::positional::boundaries::{FixedBoundary, UnsafePeriodicBoundary};
pub use crate::positional::edge::Edge;
pub use crate::positional::neighbourhood::{MooreNeighbourhood, Neighbourhood, VonNeumannNeighbourhood};
pub use crate::positional::pos::Pos;
pub use crate::positional::rect::Rect;
pub use crate::spin::Spin;
pub use crate::static_adhesion::StaticAdhesion;
pub use crate::symmetric_table::SymmetricTable;
pub use crate::traits::adhesion_system::AdhesionSystem;
pub use crate::traits::cellular::{Alive, Cellular};
pub use crate::traits::habitable::Habitable;
pub use crate::traits::potts_algorithm::PottsAlgorithm;
pub use crate::traits::step::Step;
