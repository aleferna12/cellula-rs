//! Contains logic associated with the default [`ParquetReader`].
//!
use crate::empty_cell::Empty;
use crate::io::read::read_trait::Read;
use crate::lattice::Lattice;
use crate::prelude::{CellContainer, Pos, RelCell};
use crate::spin::Spin;
use arrow::array::{Array, RecordBatch, StringArray};
use arrow::compute::concat_batches;
use arrow::error::ArrowError;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::errors::ParquetError;
use parquet::file::reader::ChunkReader;
use serde::de::DeserializeOwned;
use serde_arrow::from_record_batch;
use std::convert::Infallible;
use std::error::Error;
use std::num::ParseIntError;

/// Used to read data stored in the [`Parquet`](parquet) format.
#[derive(Clone, Debug)]
pub struct ParquetReader<R> {
    /// Object responsible for reading data.
    pub reader: R
}

impl<R: ChunkReader + 'static> Read<Lattice<Spin>, LatticeReadError> for ParquetReader<R> {
    fn read(self) -> Result<Lattice<Spin>, LatticeReadError> {
        read_lattice::<Spin, StringArray, _>(self.reader)
    }
}

impl<C, R: ChunkReader + 'static> Read<CellContainer<C>, CellsReadError> for ParquetReader<R>
where
    C: DeserializeOwned + Empty {
    fn read(self) -> Result<CellContainer<C>, CellsReadError> {
        let batches = read_parquet(self.reader)?;
        let Some(first) = batches.first() else {
            return Ok(CellContainer::new());
        };
        let schema = first.schema();
        let batch = concat_batches(&schema, &batches)?;
        let vec: Vec<RelCell<C>> = from_record_batch(&batch)?;
        let size = vec.iter().map(|rel_cell| rel_cell.index).max().expect("empty vec") + 1;
        let mut cells = CellContainer::with_capacity(size as usize);
        for _ in 0..size {
            cells.push(C::empty_default());
        }
        for rel_cell in vec {
            cells.replace(rel_cell);
        }
        Ok(cells)
    }
}

/// Defines a protocol to map the implementor into an iterator over `T`s.
pub trait MapFromArray<T, E> {
    /// Executes the mapping.
    fn map_from(&self) -> impl Iterator<Item = Result<T, E>>;
}

/// Error thrown when reading a file with cells fails.
#[derive(thiserror::Error, Debug)]
pub enum CellsReadError {
    /// Failed to read Parquet file.
    #[error(transparent)]
    Parquet(#[from] ParquetError),

    /// Operation using underlying arrow structures failed.
    #[error(transparent)]
    Arrow(#[from] ArrowError),

    /// Failed to deserialize cell.
    #[error(transparent)]
    SerdeArrow(#[from] serde_arrow::Error),
}

/// Error thrown when reading a lattice from a file fails.
#[derive(thiserror::Error, Debug)]
pub enum LatticeReadError {
    /// Found values of incorrect type in the file.
    #[error("encountered a value with invalid type")]
    InvalidType,

    /// File was empty.
    #[error("the file contained no records")]
    EmptyFile,

    /// Lattice representation was not a square.
    #[error("lattice width ({width}) and height ({height}) do not match")]
    NotSquare {
        /// Width of the lattice representation.
        width: usize,
        /// Height of the lattice representation.
        height: usize,
    },

    /// Found invalid value in the lattice.
    #[error("encountered an invalid value: {0}")]
    InvalidValue(#[source] Box<dyn Error + Send + Sync>),

    /// Found null value in the lattice.
    #[error(transparent)]
    Null(#[from] NullError),

    /// Failed to write Parquet file.
    #[error(transparent)]
    Parquet(#[from] ParquetError),
}

/// Error thrown when a [`Spin`] could not be parsed.
#[derive(thiserror::Error, Debug)]
pub enum SpinParseError {
    /// Failed to parse integer type.
    #[error(transparent)]
    Parse(#[from] ParseIntError),

    /// Found null value.
    #[error(transparent)]
    Null(#[from] NullError)
}

#[derive(thiserror::Error, Debug)]
#[error("encountered null value")]
/// Error thrown when encountering a null value where that should not be possible.
pub struct NullError;

impl MapFromArray<String, NullError> for StringArray {
    fn map_from(&self) -> impl Iterator<Item = Result<String, NullError>> {
        self.iter().map(|maybe_s| match maybe_s {
            Some(s) => Ok(s.to_string()),
            None => Err(NullError),
        })
    }
}

impl MapFromArray<Option<String>, Infallible> for StringArray {
    fn map_from(&self) -> impl Iterator<Item = Result<Option<String>, Infallible>> {
        self.iter().map(|maybe_s| Ok(maybe_s.map(|s| s.to_string())))
    }
}

impl MapFromArray<Spin, SpinParseError> for StringArray {
    fn map_from(&self) -> impl Iterator<Item = Result<Spin, SpinParseError>> {
        self.iter().map(|maybe_s| match maybe_s {
            Some(s) => Ok(match s {
                "m" => Spin::Medium,
                "s" => Spin::Solid,
                c => {
                    let cell_index = c.parse()?;
                    Spin::Some(cell_index)
                }
            }),
            None => Err(NullError.into()),
        })
    }
}

fn read_lattice<T, A, E>(file: impl ChunkReader + 'static) -> Result<Lattice<T>, LatticeReadError>
where
    T: Clone + Default,
    A: 'static + Array + MapFromArray<T, E>,
    E: 'static + Error + Sync + Send {
    let batches = read_parquet(file)?;
    if batches.is_empty() {
        return Err(LatticeReadError::EmptyFile);
    }

    let (width, height) = batches.iter().fold((0, 0), |count, batch| (
        count.0 + batch.num_rows(),
        count.1 + batch.num_columns()
    ));
    if width != height {
        return Err(LatticeReadError::NotSquare { width, height });
    }

    let mut row_offset = 0;
    let mut lat = Lattice::<T>::new(width, height);
    for batch in batches {
        for (j, col) in batch.columns().iter().enumerate() {
            let Some(col_array) = col.as_any().downcast_ref::<A>() else {
                return Err(LatticeReadError::InvalidType);
            };
            for (i, spin) in col_array.map_from().enumerate() {
                // We gotta do this again because is_nullable can lie...
                match spin {
                    Ok(valid_spin) => lat[Pos::new(row_offset + i, j)] = valid_spin,
                    Err(e) => return Err(LatticeReadError::InvalidValue(Box::new(e)))
                }
            }
        }
        row_offset += batch.num_rows();
    }
    Ok(lat)
}

fn read_parquet(file: impl ChunkReader + 'static) -> Result<Vec<RecordBatch>, ParquetError> {
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
    let reader = builder.build()?;

    let mut batches = Vec::new();
    for batch in reader {
        batches.push(batch?);
    }
    Ok(batches)
}

mod impls {
    use super::*;
    use arrow::array::*;

    macro_rules! impl_read_lat_primitive {
        ( $( ($t1:ty, $t2:ty) ),* $(,)? ) => {
            $(
                impl MapFromArray<Option<$t1>, Infallible> for $t2 {
                    fn map_from(&self) -> impl Iterator<Item = Result<Option<$t1>, Infallible>> {
                        self.iter().map(|x| Ok(x))
                    }
                }

                impl MapFromArray<$t1, NullError> for $t2 {
                    fn map_from(&self) -> impl Iterator<Item = Result<$t1, NullError>> {
                        self.iter().map(|x| x.ok_or(NullError))
                    }
                }

                impl<R: ChunkReader + 'static> Read<Lattice<$t1>, LatticeReadError> for ParquetReader<R> {
                    fn read(self) -> Result<Lattice<$t1>, LatticeReadError> {
                        read_lattice::<$t1, $t2, _>(self.reader)
                    }
                }
            )*
        };
    }

    impl_read_lat_primitive![
        (i8,  Int8Array),
        (i32, Int32Array),
        (i64, Int64Array),
        (u8,  UInt8Array),
        (u32, UInt32Array),
        (u64, UInt64Array),
        (f32, Float32Array),
    ];
}