use crate::cell_container::{CellContainer, RelCell};
use crate::io::write::writer::{Write, Writer};
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
use std::path::Path;
use std::sync::Arc;

// TODO!: make impls for primitive lattice types
impl Write<Lattice<Spin>, ParquetError> for Writer {
    fn write(&mut self, data: &Lattice<Spin>, path: impl AsRef<Path>) -> Result<(), ParquetError> {
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
        write_parquet(path, &batch).map(|_| ())
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

#[derive(thiserror::Error, Debug)]
pub enum CellsWriteError {
    #[error(transparent)]
    Parquet(#[from] ParquetError),

    #[error(transparent)]
    SerdeArrow(#[from] serde_arrow::Error),
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