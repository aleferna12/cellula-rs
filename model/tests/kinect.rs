use image::Luma;
use imageproc::definitions::Image;
use model::kinect::{kinect_create, kinect_destroy, kinect_next_depth, kinect_release_frame};

#[test]
fn test_kinect() {
    unsafe {
        let h = kinect_create();
        assert!(!h.is_null());

        let data = kinect_next_depth(h, 10_000);
        let data_arr = std::ptr::read(data as *const [f32; 512 * 424]);
        let image = Image::from_fn(
            512,
            424,
            |i, j| {
                let mut val = data_arr[(j * 512 + i) as usize];
                if val != 0. && !val.is_normal() {
                    val = 0.;
                }
                // Dividing by 20 as a hack to get more values into the appropriate range
                Luma([(val / 20.) as u8])
            },
        );
        image.save("tests/test_kinect.png").unwrap();
        kinect_release_frame(h);
        kinect_destroy(h);
    }
}