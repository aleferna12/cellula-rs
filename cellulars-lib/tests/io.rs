#![cfg(any(feature = "data-io", feature = "image-io"))]

#[cfg(feature = "image-io")]
mod image_test {
    use image::RgbaImage;
    use cellulars_lib::io::write::writer::{Write, Writer};

    #[test]
    fn test_image_writer() {
        let image = RgbaImage::from_pixel(10, 10, [255, 0, 0, 255].into());
        Writer {}.write(&image, "tests/out/image.webp").unwrap();
    }
}

#[cfg(feature = "data-io")]
mod data_test {
    use crate::data_test::cell_container::CellContainer;
    use cellulars_lib::cell_container;
    use cellulars_lib::io::read::reader::{Read, Reader};
    use cellulars_lib::io::write::writer::{Write, Writer};
    use cellulars_lib::prelude::{Cell, Lattice, RelCell, Spin};

    #[test]
    fn test_lattice_io() {
        let mut lattice = Lattice::<Spin>::new(10, 10);
        lattice[(0, 0).into()] = Spin::Some(0);
        let path = "tests/out/lattice.parquet";
        Writer {}.write(&lattice, path).unwrap();
        let lattice2 = Reader {}.read(path).unwrap();
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
        Writer {}.write(&cells, path).unwrap();
        let cells2 = Reader {}.read(path).unwrap();
        assert_eq!(cells, cells2);
    }
}