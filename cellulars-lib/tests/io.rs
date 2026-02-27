#[cfg(feature = "image-io")]
mod image_test {
    use std::fs::File;
    use image::RgbaImage;
    use cellulars_lib::io::write::image::webp_writer::WebpWriter;
    use cellulars_lib::io::write::write::Write;

    #[test]
    fn test_image_writer() {
        let image = RgbaImage::from_pixel(10, 10, [255, 0, 0, 255].into());
        WebpWriter::new(File::create("tests/out/image.webp").unwrap()).write(&image).unwrap();
    }
}

#[cfg(feature = "data-io")]
mod data_test {
    use crate::data_test::cell_container::CellContainer;
    use cellulars_lib::cell_container;
    use cellulars_lib::io::read::parquet_reader::ParquetReader;
    use cellulars_lib::io::read::read::Read;
    use cellulars_lib::prelude::{Cell, Lattice, RelCell, Spin};
    use std::fs::File;
    use cellulars_lib::io::write::parquet_writer::ParquetWriter;
    use cellulars_lib::io::write::write::Write;

    #[test]
    fn test_lattice_io() {
        let mut lattice = Lattice::<Spin>::new(10, 10);
        lattice[(0, 0).into()] = Spin::Some(0);
        let path = "tests/out/lattice.parquet";
        ParquetWriter::new(File::create(path).unwrap()).write(&lattice).unwrap();
        let lattice2 = ParquetReader::new(File::open(path).unwrap()).read().unwrap();
        assert_eq!(lattice, lattice2);
    }
    
    #[test]
    fn test_cells_io() {
        let mut cells = cell_container![Cell::new_empty(10); 5];
        let cell = Cell::new_ready(10, 10, (10., 10.).into());
        cells.replace(RelCell {
            index: 0,
            cell: cell.clone(),
        });
        cells.replace(RelCell {
            index: 4,
            cell,
        });
        let path = "tests/out/cells.parquet";
        ParquetWriter::new(File::create(path).unwrap()).write(&cells).unwrap();
        let cells2 = ParquetReader::new(File::open(path).unwrap()).read().unwrap();
        assert_eq!(cells, cells2);
    }
}