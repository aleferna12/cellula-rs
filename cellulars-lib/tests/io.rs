use image::RgbaImage;
use cellulars_lib::io::data_writer::{DataWriter, WriteData};
use cellulars_lib::lattice::Lattice;
use cellulars_lib::spin::Spin;

fn dw() -> DataWriter {
    DataWriter { outdir: "tests".into() }
}

#[test]
fn test_image_writer() {
    let mut image = RgbaImage::from_pixel(10, 10, [255, 0, 0, 255].into());
    dw().write(&mut image, 0).unwrap();
}

#[test]
fn test_lattice_writer() {
    let mut lattice = Lattice::<Spin>::new(1000, 1000);
    lattice[(0, 0).into()] = Spin::Some(0);
    dw().write(&mut lattice, 0).unwrap();
}