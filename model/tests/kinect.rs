use image::{Luma, Rgb};
use imageproc::definitions::Image;
use model::kinect::*;

#[test]
fn test_kinect() {
    unsafe {
        let h = kinect_create(true, true);
        assert!(!h.is_null());

        assert!(kinect_listen_frame(h, 10_000));

        let depth_data = kinect_depth(h);
        assert!(!depth_data.is_null());
        let depth_slice = std::slice::from_raw_parts(depth_data, 512 * 424);
        let depth_image = Image::from_fn(
            512,
            424,
            |i, j| {
                let mut val = depth_slice[(j * 512 + i) as usize];
                if val != 0. && !val.is_normal() {
                    val = 0.;
                }
                // Dividing by 20 as a hack to get more values into the appropriate range
                Luma([(val / 20.) as u8])
            },
        );

        let ir_data = kinect_ir(h);
        assert!(!ir_data.is_null());
        let ir_slice = std::slice::from_raw_parts(ir_data, 512 * 424);
        let ir_image = Image::from_fn(
            512,
            424,
            |i, j| {
                let mut val = ir_slice[(j * 512 + i) as usize];
                if val != 0. && !val.is_normal() {
                    val = 0.;
                }
                // Dividing by 100 as a hack to get more values into the appropriate range
                Luma([(val / 100.) as u8])
            },
        );

        let color_data = kinect_color(h);
        assert!(!color_data.is_null());
        let color_slice = std::slice::from_raw_parts(color_data as *const [u8; 4], 1920 * 1080);
        let color_image = Image::from_fn(
            1920,
            1080,
            |i, j| {
                let val = color_slice[(j * 1920 + i) as usize];
                Rgb([val[2], val[1], val[0]])
            },
        );

        kinect_release_frame(h);
        kinect_destroy(h);

        depth_image.save("tests/test_kinect_depth.png").unwrap();
        ir_image.save("tests/test_kinect_ir.png").unwrap();
        color_image.save("tests/test_kinect_rgb.png").unwrap();
    }
}