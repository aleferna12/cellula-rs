use pyo3::prelude::*;
mod contact;

/// Rust backend module for analyses.
#[pymodule]
mod rust {
    #[pymodule_export]
    use super::contact::contact;
}
