#![cfg(any(feature = "data-io", feature = "image-io"))]

use cellulars_lib::io::writer::Writer;

fn dw() -> Writer {
    Writer { outdir: "tests".into() }
}

#[cfg(feature = "image-io")]
mod image_test {
    use super::*;
    use cellulars_lib::io::writer::Write;
    use image::RgbaImage;

    #[test]
    fn test_image_writer() {
        let image = RgbaImage::from_pixel(10, 10, [255, 0, 0, 255].into());
        dw().write(&image, 0).unwrap();
    }
}

#[cfg(feature = "data-io")]
mod data_test {
    use super::*;
    use crate::data_test::cell_container::CellContainer;
    use cellulars_lib::cell_container;
    use cellulars_lib::io::writer::Write;
    use cellulars_lib::prelude::{Cell, Lattice, RelCell, Spin};

    #[test]
    fn test_lattice_writer() {
        let mut lattice = Lattice::<Spin>::new(1000, 1000);
        lattice[(0, 0).into()] = Spin::Some(0);
        dw().write(&lattice, 0).unwrap();
    }
    
    #[test]
    fn test_cells_writer() {
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
        dw().write(&cells, 0).unwrap();
    }
}