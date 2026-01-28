use crate::io::file::{pad_file_name, U32_STR_LEN};
use crate::prelude::{Lattice, Pos, Spin};
use arrow_array::{ArrayRef, RecordBatch, StringArray};
use image::{ImageError, RgbaImage};
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, ZstdLevel};
use parquet::errors::ParquetError;
use parquet::file::properties::WriterProperties;
use std::fs::File;
use std::path::PathBuf;
use std::sync::Arc;

pub trait WriteData<D, E> {
    fn write(&mut self, data: &mut D, time_step: u32) -> Result<PathBuf, E>;
}

pub struct DataWriter {
    pub outdir: PathBuf
}

impl DataWriter {
    fn file_path(&self, subfolder: &str, ext: &str, time_step: u32) -> Option<PathBuf> {
        let padded = pad_file_name(
            &format!("{time_step}.{ext}"),
            U32_STR_LEN
        )?;
        Some(self.outdir.join(subfolder).join(padded))
    }
}

impl WriteData<RgbaImage, ImageError> for DataWriter {
    fn write(&mut self, data: &mut RgbaImage, time_step: u32) -> Result<PathBuf, ImageError> {
        let file_path = self.file_path(
            "images",
            "webp",
            time_step
        ).expect("failed to pad time step when saving image");  // This should never fail
        data.save(&file_path).map(|_| file_path)
    }
}

impl WriteData<Lattice<Spin>, ParquetError> for DataWriter {
    fn write(&mut self, data: &mut Lattice<Spin>, time_step: u32) -> Result<PathBuf, ParquetError> {
        let file_path = self.file_path(
            "lattices",
            "parquet",
            time_step
        ).expect("failed to pad time step when saving cell lattice");  // This should never fail

        let batch = RecordBatch::try_from_iter_with_nullable(
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
                (j.to_string(), Arc::new(arr) as ArrayRef, false)
            })
        )?;

        let props = WriterProperties::builder()
            .set_compression(Compression::ZSTD(ZstdLevel::default()))
            .build();
        let file = File::create(&file_path)?;
        let mut writer = ArrowWriter::try_new(file, batch.schema(), Some(props))?;
        writer.write(&batch)?;
        writer.close().map(|_| file_path)
    }
}