use pyo3::prelude::*;
mod contact;

/// Rust backend module for analyses.
#[pymodule]
mod rust {
    use pyo3::prelude::*;

    #[pymodule_export]
    use super::contact::contact;
}
