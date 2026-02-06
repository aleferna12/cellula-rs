use crate::lattice::Lattice;
use crate::prelude::{CellIndex, Pos};
use crate::spin::Spin;
use arrow::array::{Array, RecordBatch, StringArray};
use arrow::error::ArrowError;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::errors::ParquetError;
use std::error::Error;
use std::fs::File;
use std::path::Path;
use std::str::FromStr;

pub struct Reader {}

impl Reader {
    fn parse_spin(spin_str: String) -> Result<Spin, <CellIndex as FromStr>::Err> {
        Ok(match spin_str.as_str() {
            "m" => Spin::Medium,
            "s" => Spin::Solid,
            c => {
                let cell_index = c.parse()?;
                Spin::Some(cell_index)
            }
        })
    }
}

impl Read<Lattice<Spin>, LatticeReadError> for Reader {
    fn read(&mut self, path: impl AsRef<Path>) -> Result<Lattice<Spin>, LatticeReadError> {
        read_lattice::<StringArray, _, _, _, _>(path, Self::parse_spin)
    }
}

fn read_lattice<A, T, U, F, E>(path: impl AsRef<Path>, map: F) -> Result<Lattice<U>, LatticeReadError>
where
    A: 'static + Array + ColumnIter<T>,
    U: Clone + Default,
    F: Fn(T) -> Result<U, E>,
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
    let mut lat = Lattice::<U>::new(width, height);
    for batch in batches {
        for (j, col) in batch.columns().iter().enumerate() {
            let Some(col_array) = col.as_any().downcast_ref::<A>() else {
                return Err(LatticeReadError::InvalidType);
            };
            if col_array.is_nullable() {
                return Err(LatticeReadError::NullValue);
            }
            for (i, maybe_spin_str) in col_array.iter().enumerate() {
                // We gotta do this again because is_nullable can lie...
                let Some(spin_val) = maybe_spin_str else {
                    return Err(LatticeReadError::NullValue);
                };
                let spin = map(spin_val).map_err(|e| LatticeReadError::InvalidValue(Box::new(e)))?;
                lat[Pos::new(row_offset + i, j)] = spin;
            }
        }
        row_offset += batch.num_rows();
    }
    Ok(lat)
}

#[derive(thiserror::Error, Debug)]
pub enum LatticeReadError {
    #[error("encountered a null value in the file")]
    NullValue,
    #[error("encountered a value with invalid type")]
    InvalidType,
    #[error("encountered an invalid value: {0}")]
    InvalidValue(#[source] Box<dyn Error>),
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

trait ColumnIter<T> {
    fn iter(&self) -> impl Iterator<Item = Option<T>>;
}

impl ColumnIter<String> for StringArray {
    fn iter(&self) -> impl Iterator<Item = Option<String>> {
        self.iter().map(|x| x.map(ToOwned::to_owned))
    }
}

mod impls {
    use super::*;
    use arrow::array::*;

    macro_rules! impl_column_iter_primitive {
        ( $( ($t1:ty, $t2:ty) ),* $(,)? ) => {
            $(
                impl ColumnIter<$t1> for $t2 {
                    fn iter(&self) -> impl Iterator<Item = Option<$t1>> {
                        self.iter()
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