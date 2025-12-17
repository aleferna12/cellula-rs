use pyo3::prelude::*;
use pyo3_polars::PyDataFrame;

#[pymodule]
pub mod contact {
    #[pymodule_export]
    use super::geom_contacts;
}

#[pyfunction]
pub fn geom_contacts() {
    todo!()
}

pub struct ActVecs {
    pub cell: Vec<u32>,
    pub medium: Vec<u32>,
    pub solid: Vec<u32>
}