use crate::cell_container;
use crate::lattice::Lattice;
use crate::prelude::{CellContainer, Pos, RelCell};
use crate::spin::Spin;
use crate::traits::cellular::EmptyCell;
use arrow::array::{Array, RecordBatch, StringArray};
use arrow::compute::concat_batches;
use arrow::error::ArrowError;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::errors::ParquetError;
use serde::de::DeserializeOwned;
use serde_arrow::from_record_batch;
use std::error::Error;
use std::fs::File;
use std::num::ParseIntError;
use std::path::Path;

pub struct Reader {}

impl Reader {}

impl Read<Lattice<Spin>, LatticeReadError> for Reader {
    fn read(&mut self, path: impl AsRef<Path>) -> Result<Lattice<Spin>, LatticeReadError> {
        read_lattice::<Spin, StringArray, _>(path)
    }
}

impl<C> Read<CellContainer<C>, CellsReadError> for Reader
where
    C: DeserializeOwned,
    EmptyCell<C>: Default + Clone {
    fn read(&mut self, path: impl AsRef<Path>) -> Result<CellContainer<C>, CellsReadError> {
        let batches = read_parquet(path)?;
        let Some(first) = batches.first() else {
            return Ok(CellContainer::new());
        };
        let schema = first.schema();
        let batch = concat_batches(&schema, &batches)?;
        let vec: Vec<RelCell<C>> = from_record_batch(&batch)?;
        let size = vec.iter().map(|rel_cell| rel_cell.index).max().expect("empty vec") + 1;
        let mut cells = cell_container![EmptyCell::<C>::default(); size as usize];
        for rel_cell in vec {
            cells.replace(rel_cell);
        }
        Ok(cells)
    }
}

#[derive(thiserror::Error, Debug)]
pub enum CellsReadError {
    #[error(transparent)]
    Parquet(#[from] ParquetError),

    #[error(transparent)]
    Arrow(#[from] ArrowError),

    #[error(transparent)]
    SerdeArrow(#[from] serde_arrow::Error),
}

fn read_lattice<T, A, E>(path: impl AsRef<Path>) -> Result<Lattice<T>, LatticeReadError>
where
    T: Clone + Default,
    A: 'static + Array + MapColumn<T, E>,
    E: 'static + Error + Sync + Send {
    let batches = read_parquet(path)?;
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
            for (i, spin) in col_array.map().enumerate() {
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

#[derive(thiserror::Error, Debug)]
pub enum LatticeReadError {
    #[error("encountered a value with invalid type")]
    InvalidType,
    #[error("the file contained no records")]
    EmptyFile,
    #[error("lattice width ({width}) and height ({height}) do not match")]
    NotSquare {
        width: usize,
        height: usize,
    },
    #[error("encountered an invalid value: {0}")]
    InvalidValue(#[source] Box<dyn Error + Send + Sync>),
    #[error(transparent)]
    Null(#[from] NullError),
    #[error(transparent)]
    Parquet(#[from] ParquetError),
}

fn read_parquet(path: impl AsRef<Path>) -> Result<Vec<RecordBatch>, ParquetError> {
    let file = File::open(path)?;
    let builder = ParquetRecordBatchReaderBuilder::try_new(file)?;
    let reader = builder.build()?;

    let mut batches = Vec::new();
    for batch in reader {
        batches.push(batch?);
    }
    Ok(batches)
}

pub trait Read<D, E> {
    fn read(&mut self, path: impl AsRef<Path>) -> Result<D, E>;
}

pub trait MapColumn<T, E> {
    fn map(&self) -> impl Iterator<Item = Result<T, E>>;
}

impl MapColumn<Option<String>, ()> for StringArray {
    fn map(&self) -> impl Iterator<Item = Result<Option<String>, ()>> {
        self.iter().map(|maybe_s| Ok(maybe_s.map(|s| s.to_string())))
    }
}

impl MapColumn<String, NullError> for StringArray {
    fn map(&self) -> impl Iterator<Item = Result<String, NullError>> {
        self.iter().map(|maybe_s| match maybe_s {
            Some(s) => Ok(s.to_string()),
            None => Err(NullError),
        })
    }
}

impl MapColumn<Spin, SpinParseError> for StringArray {
    fn map(&self) -> impl Iterator<Item = Result<Spin, SpinParseError>> {
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

#[derive(thiserror::Error, Debug)]
pub enum SpinParseError {
    #[error(transparent)]
    Parse(#[from] ParseIntError),
    #[error(transparent)]
    Null(#[from] NullError)
}

#[derive(thiserror::Error, Debug)]
#[error("encountered null value")]
pub struct NullError;

mod impls {
    use super::*;
    use arrow::array::*;

    macro_rules! impl_read_lat_primitive {
        ( $( ($t1:ty, $t2:ty) ),* $(,)? ) => {
            $(
                impl MapColumn<Option<$t1>, ()> for $t2 {
                    fn map(&self) -> impl Iterator<Item = Result<Option<$t1>, ()>> {
                        self.iter().map(|x| Ok(x))
                    }
                }

                impl MapColumn<$t1, NullError> for $t2 {
                    fn map(&self) -> impl Iterator<Item = Result<$t1, NullError>> {
                        self.iter().map(|x| x.ok_or(NullError))
                    }
                }

            impl Read<Lattice<$t1>, LatticeReadError> for Reader {
                fn read(&mut self, path: impl AsRef<Path>) -> Result<Lattice<$t1>, LatticeReadError> {
                    read_lattice::<$t1, $t2, _>(path)
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