use crate::constants::CellIndex;
use std::fmt::{Debug, Display, Formatter};
use std::num::ParseIntError;

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

pub fn spin_to_str(spin: Spin) -> String {
    match spin {
        Spin::Solid => String::from("s"),
        Spin::Medium => String::from("m"),
        Spin::Some(cell_index) => cell_index.to_string(),
    }
}

pub fn str_to_spin(s: &str) -> Result<Spin, ParseIntError> {
    Ok(match s {
        "s" => Spin::Solid,
        "m" => Spin::Medium,
        _ => {
            let cell_index = s.parse::<CellIndex>()?;
            Spin::Some(cell_index)
        },
    })
}