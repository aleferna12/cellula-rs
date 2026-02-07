use crate::cell_container::{CellContainer, RelCell};
use crate::io::write::writer::{Write, Writer};
use crate::lattice::Lattice;
use crate::prelude::{Cellular, Pos, Spin};
use arrow::array::{Array, ArrayRef, RecordBatch, StringArray, UInt32Array};
use arrow::datatypes::FieldRef;
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

impl Write<Lattice<Spin>, ParquetError> for Writer {
    fn write(&mut self, data: &Lattice<Spin>, path: impl AsRef<Path>) -> Result<(), ParquetError> {
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

pub trait MapIntoColumn<A, E> {
    fn map(self) -> Result<A, E>;
}

impl MapIntoColumn<StringArray, ()> for Vec<String> {
    fn map(self) -> Result<StringArray, ()> {
        Ok(StringArray::from(self))
    }
}

impl MapIntoColumn<StringArray, ()> for Vec<Option<String>> {
    fn map(self) -> Result<StringArray, ()> {
        Ok(StringArray::from(self))
    }
}

#[derive(thiserror::Error, Debug)]
pub enum CellsWriteError {
    #[error(transparent)]
    Parquet(#[from] ParquetError),

    #[error(transparent)]
    SerdeArrow(#[from] serde_arrow::Error),
}

fn f() {
    match spin {
        // Saving as strings with one char is more compact than u32 (or longer strings)
        Spin::Medium => "m".to_string(),
        Spin::Solid => "s".to_string(),
        Spin::Some(cell_index) => cell_index.to_string()
    }
}

fn write_lattice<T, A, E>(path: impl AsRef<Path>, data: &Lattice<T>) -> Result<ParquetMetaData, E>
where
    Vec<T>: MapIntoColumn<A, E>,
    A: Array {
    let batch = RecordBatch::try_from_iter(
        // TODO!: This needs to short circuit
        (0..data.width()).find_map(move |j| {
            let vec: Vec<T> = (0..data.height()).map(|i| {
                data[Pos::new(i, j)]
            }).collect();
            let arr = vec.map()?;
            (j.to_string(), Arc::new(arr) as ArrayRef)
        })
    )?;
    write_parquet(path, &batch).map(|_| ())
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