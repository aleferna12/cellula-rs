use crate::lattice::Lattice;
use crate::prelude::Pos;
use crate::spin::Spin;
use arrow::array::{Array, RecordBatch, StringArray, UInt32Array};
use arrow::error::ArrowError;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::errors::ParquetError;
use std::fs::File;
use std::path::Path;

pub struct Reader {}

impl Reader {
    fn parse_spin(spin_str: &str) -> Result<Spin, String> {
        Ok(match spin_str {
            "m" => Spin::Medium,
            "s" => Spin::Solid,
            c => {
                let Ok(cell_index) = c.parse() else {
                    return Err(c.to_string());
                };
                Spin::Some(cell_index)
            }
        })
    }
}

impl Read<Lattice<Spin>, LatticeReadError<String>> for Reader {
    fn read(&mut self, path: impl AsRef<Path>) -> Result<Lattice<Spin>, LatticeReadError<String>> {
        read_lattice::<&str, Spin, StringArray, _, String>(path, Self::parse_spin)
    }
}

trait ColumnIter<'a, T> {
    type Iter: Iterator<Item = Option<T>>;

    fn iter(&'a self) -> Self::Iter;
}

impl<'a> ColumnIter<'a, u32> for UInt32Array {
    type Iter = std::iter::Map<
        std::slice::Iter<'a, u32>,
        fn(&u32) -> Option<u32>
    >;

    fn iter(&'a self) -> Self::Iter {
        self.iter()
    }
}

impl<'a> ColumnIter<'a, String> for StringArray {
    type Iter = std::iter::Map<std::slice::Iter<'a, &'a str>, fn(&&str) -> Option<String>>;

    fn iter(&'a self) -> Self::Iter {
        self.iter().map(|x| x.map(ToOwned::to_owned))
    }
}

fn read_lattice<T, U, A, F, E>(path: impl AsRef<Path>, map: F) -> Result<Lattice<U>, LatticeReadError<E>>
where
    A: 'static + Array,
    for<'a> &'a A: IntoIterator<Item = Option<T>>,
    F: Fn(T) -> Result<U, E>,
    U: Clone + Default {
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
    let mut lat = Lattice::<U>::new(width, height);
    for batch in batches {
        for (j, col) in batch.columns().iter().enumerate() {
            let Some(col_array) = col.as_any().downcast_ref::<A>() else {
                return Err(LatticeReadError::InvalidType);
            };
            if col_array.is_nullable() {
                return Err(LatticeReadError::NullValue);
            }
            for (i, maybe_spin_str) in col_array.into_iter().enumerate() {
                // We gotta do this again because is_nullable can lie...
                let Some(spin_val) = maybe_spin_str else {
                    return Err(LatticeReadError::NullValue);
                };
                let spin = map(spin_val).map_err(|e| LatticeReadError::InvalidValue(e))?;
                lat[Pos::new(row_offset + i, j)] = spin;
            }
        }
        row_offset += batch.num_rows();
    }
    Ok(lat)
}

#[derive(thiserror::Error, Debug)]
pub enum LatticeReadError<T> {
    #[error("encountered a null value in the file")]
    NullValue,
    #[error("encountered a value with invalid type")]
    InvalidType,
    #[error("encountered an invalid value: ({0})")]
    InvalidValue(T),
    #[error("the file contained no records")]
    EmptyFile,
    #[error("lattice width ({width}) and height ({height}) do not match")]
    NotSquare {
        width: usize,
        height: usize,
    },
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