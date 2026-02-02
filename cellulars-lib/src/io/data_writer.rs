use crate::constants::FloatType;
use crate::io::file::{pad_file_name, U32_STR_LEN};
use crate::prelude::{Cell, CellContainer, Cellular, HasCenter, Lattice, Pos, Spin};
use arrow_array::{ArrayRef, Float64Array, RecordBatch, StringArray, UInt32Array};
#[cfg(feature = "image_io")]
use image::{ImageError, RgbaImage};
use parquet::arrow::ArrowWriter;
use parquet::basic::{Compression, ZstdLevel};
use parquet::errors::ParquetError;
use parquet::file::properties::WriterProperties;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use parquet::file::metadata::ParquetMetaData;

pub trait WriteData<D, E> {
    fn write(&mut self, data: &D, time_step: u32) -> Result<PathBuf, E>;
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

#[cfg(feature = "image_io")]
impl WriteData<RgbaImage, ImageError> for DataWriter {
    fn write(&mut self, data: &RgbaImage, time_step: u32) -> Result<PathBuf, ImageError> {
        let file_path = self.file_path(
            "images",
            "webp",
            time_step
        ).expect("failed to pad time step when saving image");  // This should never fail
        data.save(&file_path).map(|_| file_path)
    }
}

impl WriteData<Lattice<Spin>, ParquetError> for DataWriter {
    fn write(&mut self, data: &Lattice<Spin>, time_step: u32) -> Result<PathBuf, ParquetError> {
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
        write_record_batch(&file_path, &batch).map(|_| file_path)
    }
}

// TODO!: this REALLY should be a proc derive macro on `Cell` that creates definitions for any IntoIter<&Cell>
//  The challenge is to implement it for all fields of complex types too like `Com` (maybe use StructArray?)
//  I think the solution would be to also write a declarative macro that takes a struct or primitive and returns an Array type
//      arrow_array!(f64)  -> Float64Array
//      arrow_array!(Com)  -> StructArray[x: Float64Array, y: Float64Array, area: U32Array]
//      arrow_array!(Cell) -> StructArray[target_area: U32Array, com: StructArray]  // This can be converted into a RecordBatch via From
impl WriteData<CellContainer<Cell>, ParquetError> for DataWriter {
    fn write(&mut self, data: &CellContainer<Cell>, time_step: u32) -> Result<PathBuf, ParquetError> {
        let file_path = self.file_path(
            "cells",
            "parquet",
            time_step
        ).expect("failed to pad time step when saving cells");  // This should never fail

        let cells: Box<_> = data.iter_non_empty().collect();
        let batch = RecordBatch::try_from_iter([
            ("target_area", cells.iter().map(|rel_cell| rel_cell.cell.target_area).into_array()),
            ("area", cells.iter().map(|rel_cell| rel_cell.cell.area()).into_array()),
            ("center_x", cells.iter().map(|rel_cell| rel_cell.cell.center().x).into_array()),
            ("center_y", cells.iter().map(|rel_cell| rel_cell.cell.center().y).into_array()),
        ])?;
        write_record_batch(&file_path, &batch).map(|_| file_path)
    }
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

trait IntoArray<T> {
    fn into_array(self) -> ArrayRef;
}

impl<T: Iterator<Item = u32>> IntoArray<u32> for T {
    fn into_array(self) -> ArrayRef {
        Arc::new(UInt32Array::from(self.collect::<Vec<T::Item>>())) as ArrayRef
    }
}

impl<T: Iterator<Item = FloatType>> IntoArray<FloatType> for T {
    fn into_array(self) -> ArrayRef {
        #[cfg(feature = "high-precision")]
        return Arc::new(Float64Array::from(self.collect::<Vec<T::Item>>())) as ArrayRef;
        #[cfg(not(feature = "high-precision"))]
        return Arc::new(Float32Array::from(self.collect::<Vec<T::Item>>())) as ArrayRef;
    }
}