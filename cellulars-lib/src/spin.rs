use crate::constants::CellIndex;
use std::fmt::{Debug, Display, Formatter};

/// This enum represents anything that can be on the cell lattice.
#[derive(Debug, Copy, Clone, PartialEq, Eq, Hash, Default)]
pub enum Spin {
    #[default]
    Medium,
    Solid,
    Some(CellIndex),
}

impl Display for Spin {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        let s = match self {
            Spin::Medium => "m",
            Spin::Solid => "s",
            Spin::Some(ci) => &ci.to_string(),
        };
        write!(f, "{s}")
    }
}