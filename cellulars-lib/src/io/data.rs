use crate::cell_container::{CellContainer, RelCell};
use crate::io::file::file_path;
use crate::io::writer::{Write, Writer};
use crate::lattice::Lattice;
use crate::prelude::{Cellular, Pos, Spin};
use arrow::array::{ArrayRef, RecordBatch, StringArray};
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
use std::path::{Path, PathBuf};
use std::sync::Arc;

#[cfg(feature = "data-io")]
impl Write<Lattice<Spin>, ParquetError> for Writer {
    fn write(&mut self, data: &Lattice<Spin>, time_step: u32) -> Result<PathBuf, ParquetError> {
        let file_path = file_path(
            self.outdir.as_path(),
            "lattices",
            "parquet",
            time_step
        ).expect("failed to pad time step when saving cell lattice");  // This should never fail

        let batch = RecordBatch::try_from_iter(
            (0..data.width()).map(move |j| {
                let vec: Vec<_> = (0..data.height()).map(|i| {
                    let spin = data[Pos::new(i, j)];
                    match spin {
                        // Saving as strings with one char is more compact than u32 (or longer strings)
                        Spin::Medium => "m".to_string(),
                        Spin::Solid => "s".to_string(),
                        Spin::Some(cell_index) => cell_index.to_string()
                    }
                }).collect();
                let arr = StringArray::from(vec);
                (j.to_string(), Arc::new(arr) as ArrayRef)
            })
        )?;
        write_record_batch(&file_path, &batch).map(|_| file_path)
    }
}

#[cfg(feature = "data-io")]
impl<'de, T> Write<CellContainer<T>, CellsWriteError> for Writer
where
    T: Cellular,
    RelCell<T>: Serialize + Deserialize<'de> {
    fn write(&mut self, data: &CellContainer<T>, time_step: u32) -> Result<PathBuf, CellsWriteError> {
        let file_path = file_path(
            self.outdir.as_path(),
            "cells",
            "parquet",
            time_step
        ).expect("failed to pad time step when saving cells");  // This should never fail

        let cells: Box<_> = data.iter_non_empty().collect();
        let fields = Vec::<FieldRef>::from_type::<RelCell<T>>(TracingOptions::default())?;
        let batch = to_record_batch(&fields, &cells)?;
        match write_record_batch(&file_path, &batch) {
            Ok(_) => Ok(file_path),
            Err(e) => Err(e.into())
        }
    }
}

#[cfg(feature = "data-io")]
#[derive(thiserror::Error, Debug)]
pub enum CellsWriteError {
    #[error(transparent)]
    Parquet(#[from] ParquetError),

    #[error(transparent)]
    SerdeArrow(#[from] serde_arrow::Error),
}

pub fn write_record_batch(path: impl AsRef<Path>, batch: &RecordBatch) -> Result<ParquetMetaData, ParquetError> {
    let file_path = path.as_ref();
    let props = WriterProperties::builder()
        .set_compression(Compression::ZSTD(ZstdLevel::default()))
        .build();
    let file = File::create(file_path)?;
    let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props))?;
    writer.write(batch)?;
    writer.close()
}