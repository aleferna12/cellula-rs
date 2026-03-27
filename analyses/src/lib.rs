use pyo3::prelude::*;
mod contact;
mod neighbor;

/// Rust backend module for analyses.
#[pymodule]
mod rust {
    #[pymodule_export]
    use super::contact::contact;
    #[pymodule_export]
    use super::neighbor::neighbor;
}
