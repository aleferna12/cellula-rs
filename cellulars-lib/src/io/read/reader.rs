use crate::lattice::Lattice;
use crate::prelude::Pos;
use crate::spin::Spin;
use arrow::array::{Array, RecordBatch, StringArray};
use arrow::error::ArrowError;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::errors::ParquetError;
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

fn read_lattice<T, A, E>(path: impl AsRef<Path>) -> Result<Lattice<T>, LatticeReadError>
where
    T: Clone + Default,
    A: 'static + Array + ValidateColumn<T, E>,
    E: 'static + Error {
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
            for (i, spin) in col_array.iter().enumerate() {
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
    #[error("encountered an invalid value: {0}")]
    InvalidValue(#[source] Box<dyn Error>),
    #[error("encountered a value with invalid type")]
    InvalidType,
    #[error("the file contained no records")]
    EmptyFile,
    #[error("lattice width ({width}) and height ({height}) do not match")]
    NotSquare {
        width: usize,
        height: usize,
    },
    #[error(transparent)]
    Null(#[from] NullError),
    #[error(transparent)]
    Parquet(#[from] ParquetError),
    #[error(transparent)]
    Arrow(#[from] ArrowError)
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

trait ValidateColumn<T, E> {
    fn iter(&self) -> impl Iterator<Item = Result<T, E>>;
}

impl ValidateColumn<Option<String>, ()> for StringArray {
    fn iter(&self) -> impl Iterator<Item = Result<Option<String>, ()>> {
        self.iter().map(|maybe_s| Ok(maybe_s.map(|s| s.to_string())))
    }
}

impl ValidateColumn<String, NullError> for StringArray {
    fn iter(&self) -> impl Iterator<Item = Result<String, NullError>> {
        self.iter().map(|maybe_s| match maybe_s {
            Some(s) => Ok(s.to_string()),
            None => Err(NullError),
        })
    }
}

impl ValidateColumn<Spin, SpinParseError> for StringArray {
    fn iter(&self) -> impl Iterator<Item = Result<Spin, SpinParseError>> {
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

// TODO! do we really need all these error types? why not just return LatticeReadErrors?
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

    macro_rules! impl_column_iter_primitive {
        ( $( ($t1:ty, $t2:ty) ),* $(,)? ) => {
            $(
                impl ValidateColumn<Option<$t1>, ()> for $t2 {
                    fn iter(&self) -> impl Iterator<Item = Result<Option<$t1>, ()>> {
                        self.iter().map(|x| Ok(x))
                    }
                }

                impl ValidateColumn<$t1, NullError> for $t2 {
                    fn iter(&self) -> impl Iterator<Item = Result<$t1, NullError>> {
                        self.iter().map(|x| x.ok_or(NullError))
                    }
                }
            )*
        };
    }

    impl_column_iter_primitive![
        (i8,  Int8Array),
        (i32, Int32Array),
        (i64, Int64Array),
        (u8,  UInt8Array),
        (u32, UInt32Array),
        (u64, UInt64Array),
        (f32, Float32Array),
    ];
}