//! Contains logic used to write data using [`Writer`].

use crate::cell_container::{CellContainer, RelCell};
use crate::io::write::writer::{Write, Writer};
use crate::lattice::Lattice;
use crate::prelude::{Cellular, Pos, Spin};
use arrow::array::{Array, ArrayRef, RecordBatch, StringArray};
use arrow::datatypes::FieldRef;
use arrow::error::ArrowError;
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, ZstdLevel};
use parquet::errors::ParquetError;
use parquet::file::metadata::ParquetMetaData;
use parquet::file::properties::WriterProperties;
use serde::{Deserialize, Serialize};
use serde_arrow::schema::{SchemaLike, TracingOptions};
use serde_arrow::to_record_batch;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

impl Write<Lattice<Spin>, LatticeWriteError> for Writer {
    fn write(&mut self, data: &Lattice<Spin>, path: impl AsRef<Path>) -> Result<(), LatticeWriteError> {
        write_lattice(data, path).map(|_| ())
    }
}

impl<'de, T> Write<CellContainer<T>, CellsWriteError> for Writer
where
    T: Cellular,
    RelCell<T>: Serialize + Deserialize<'de> {
    fn write(&mut self, data: &CellContainer<T>, path: impl AsRef<Path>) -> Result<(), CellsWriteError> {
        let cells: Box<_> = data.iter_non_empty().collect();
        let fields = Vec::<FieldRef>::from_type::<RelCell<T>>(TracingOptions::default())?;
        let batch = to_record_batch(&fields, &cells)?;
        match write_parquet(&path, &batch) {
            Ok(_) => Ok(()),
            Err(e) => Err(e.into())
        }
    }
}

/// Defines a protocol to mapping a container type into an [`arrow`] array type.
pub trait MapIntoArray<A> {
    /// Executes the mapping.
    fn map_into(self) -> A;
}

impl MapIntoArray<StringArray> for Vec<String> {
    fn map_into(self) -> StringArray {
        self.into()
    }
}

impl MapIntoArray<StringArray> for Vec<Option<String>> {
    fn map_into(self) -> StringArray {
        self.into()
    }
}

impl MapIntoArray<StringArray> for Vec<Spin> {
    fn map_into(self) -> StringArray {
        self
            .iter()
            .map(|val| match val {
                Spin::Medium => "m".into(),
                Spin::Solid => "s".into(),
                Spin::Some(ci) => ci.to_string()
            })
            .collect::<Vec<_>>()
            .into()
    }
}

/// Error thrown when writing cells to a file fails.
#[derive(thiserror::Error, Debug)]
pub enum CellsWriteError {
    /// Failed to write Parquet file.
    #[error(transparent)]
    Parquet(#[from] ParquetError),

    /// Failed to serialize the cell.
    #[error(transparent)]
    SerdeArrow(#[from] serde_arrow::Error),
}

/// Error thrown when writing a lattice fails.
#[derive(thiserror::Error, Debug)]
pub enum LatticeWriteError {
    /// Failed an operation using underlying arrow data structures.
    #[error(transparent)]
    Arrow(# [from] ArrowError),

    /// Failed to write Parquet file.
    #[error(transparent)]
    Parquet(#[from] ParquetError),
}

fn write_lattice<T, A>(data: &Lattice<T>, path: impl AsRef<Path>) -> Result<ParquetMetaData, LatticeWriteError>
where
    Vec<T>: MapIntoArray<A>,
    T: Clone,
    A: Array + 'static {
    let batch = RecordBatch::try_from_iter(
        (0..data.width()).map(move |j| {
            let vec: Vec<T> = (0..data.height()).map(|i| {
                data[Pos::new(i, j)].clone()
            }).collect();
            let arr = vec.map_into();
            (j.to_string(), Arc::new(arr) as ArrayRef)
        })
    )?;
    write_parquet(path, &batch).map_err(|e| e.into())
}

fn write_parquet(path: impl AsRef<Path>, batch: &RecordBatch) -> Result<ParquetMetaData, ParquetError> {
    let path = path.as_ref();
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(ZstdLevel::default()))
        .build();
    let file = File::create(path)?;
    let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props))?;
    writer.write(batch)?;
    writer.close()
}

mod impls {
    use super::*;
    use arrow::array::*;

    macro_rules! impl_write_lat_primitive {
        ( $( ($t1:ty, $t2:ty) ),* $(,)? ) => {
            $(
                impl MapIntoArray<$t2> for Vec<$t1> {
                    fn map_into(self) -> $t2 {
                        self.into()
                    }
                }

                impl MapIntoArray<$t2> for Vec<Option<$t1>> {
                    fn map_into(self) -> $t2 {
                        self.into()
                    }
                }

                impl Write<Lattice<$t1>, LatticeWriteError> for Writer {
                    fn write(&mut self, data: &Lattice<$t1>, path: impl AsRef<Path>) -> Result<(), LatticeWriteError> {
                        write_lattice::<$t1, $t2>(data, path).map(|_| ())
                    }
                }
            )*
        };
    }

    impl_write_lat_primitive![
        (i8,  Int8Array),
        (i32, Int32Array),
        (i64, Int64Array),
        (u8,  UInt8Array),
        (u32, UInt32Array),
        (u64, UInt64Array),
        (f32, Float32Array),
    ];
}