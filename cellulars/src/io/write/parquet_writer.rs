//! Contains logic used to write data using [`ParquetWriter`].

use crate::cell_container::{CellContainer, RelCell};
use crate::empty_cell::Empty;
use crate::io::write::write_trait::Write;
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
use serde_arrow::marrow::datatypes::Field;
use serde_arrow::schema::{SchemaLike, TracingOptions};
use serde_arrow::to_record_batch;
use std::io;
use std::sync::Arc;

/// Parquet writer.
#[derive(Clone, Debug)]
pub struct ParquetWriter<W> {
    /// Object responsible for reading data.
    pub writer: W,

    /// These are passed to [`TracingOptions::overwrite()`].
    ///
    /// First element in the tuple is the path and the second is the field.
    pub overwrites: Vec<(String, Field)>
}

impl<W: io::Write + Send> Write<Lattice<Spin>, LatticeWriteError> for ParquetWriter<W> {
    fn write(self, data: &Lattice<Spin>) -> Result<(), LatticeWriteError> {
        write_lattice(data, self.writer).map(|_| ())
    }
}

impl<'de, T, W: io::Write + Send> Write<CellContainer<T>, CellsWriteError> for ParquetWriter<W>
where
    T: Cellular + Empty,
    RelCell<T>: Serialize + Deserialize<'de> {
    fn write(self, data: &CellContainer<T>) -> Result<(), CellsWriteError> {
        let cells: Box<_> = data.iter_non_empty().collect();
        let mut options = TracingOptions::default()
            .allow_null_fields(true)
            .allow_to_string(true)
            .coerce_numbers(true)
            .enums_without_data_as_strings(true);
        for (path, field) in self.overwrites {
            options = options.overwrite(path, field)?;
        }
        let fields = Vec::<FieldRef>::from_type::<RelCell<T>>(options)?;
        let batch = to_record_batch(&fields, &cells)?;
        match write_parquet(&batch, self.writer) {
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

fn write_lattice<T, A>(data: &Lattice<T>, file: impl io::Write + Send) -> Result<ParquetMetaData, LatticeWriteError>
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
    write_parquet(&batch, file).map_err(|e| e.into())
}

fn write_parquet(batch: &RecordBatch, file: impl io::Write + Send) -> Result<ParquetMetaData, ParquetError> {
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(ZstdLevel::default()))
        .build();
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

                impl<W: io::Write + Send> Write<Lattice<$t1>, LatticeWriteError> for ParquetWriter<W> {
                    fn write(self, data: &Lattice<$t1>) -> Result<(), LatticeWriteError> {
                        write_lattice::<$t1, $t2>(data, self.writer).map(|_| ())
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