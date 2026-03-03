use crate::cell_container::CellContainer;
use cellulars::cell_container;
use cellulars::empty_cell::EmptyCell;
use cellulars::io::write::parquet_writer::ParquetWriter;
use cellulars::io::write::r#trait::Write;
use cellulars::prelude::{Cellular, Pos, Rect, UnsafePeriodicBoundary};
use model::my_cell::{CellType, MyCell};
use std::fs::File;
use cellulars::io::read::parquet_reader::ParquetReader;
use cellulars::io::read::r#trait::Read;

#[test]
fn test_cells_io() {
    let mut cell = MyCell::new_empty(0, 0, CellType::Migrating).into_cell();
    cell.shift_position(
        Pos::new(0, 0),
        true,
        &UnsafePeriodicBoundary::new(Rect::new(
            Pos::new(0., 0.),
            Pos::new(10., 10.)
        ))).unwrap();
    let cells = cell_container![EmptyCell::new_unchecked(cell)];
    let path = "tests/out/cells.parquet";
    ParquetWriter {
        writer: File::create(path).unwrap(),
        overwrites: vec![]
    }.write(&cells).unwrap();
    let cells: CellContainer<MyCell> = ParquetReader { reader: File::open(path).unwrap() }.read().unwrap();
    assert_eq!(cells[0].index, 0)
}