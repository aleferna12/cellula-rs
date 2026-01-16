use cellulars_lib::prelude::{Habitable, Pos};
use cellulars_lib::spin::Spin;
use crate::environment::Environment;
use crate::kinect::{kinect_create, kinect_destroy, kinect_next_depth, kinect_release_frame, KinectHandle};

pub struct KinectListener {
    handle: *mut KinectHandle,
    pub min_depth: f32,
    pub max_depth: f32,
    pub frame_period: u32,
}

impl KinectListener {
    pub fn new(min_depth: f32, max_depth: f32, frame_period: u32) -> Option<Self> {
        let handle = unsafe { kinect_create() };
        if handle.is_null() {
            return None;
        }
        Some(Self {
            handle,
            min_depth,
            max_depth,
            frame_period,
        })
    }

    pub fn draw_silhouette(&mut self, env: &mut Environment) -> anyhow::Result<u32> {
        if self.handle.is_null() { anyhow::bail!("kinect handle was lost") };
        if env.base_env.width() != 512 || env.base_env.height() != 424 {
            anyhow::bail!("kinect can only be used if pond's width is 512 and height is 424");
        }

        let data_arr = unsafe {
            let data = kinect_next_depth(self.handle, 10_000);
            std::ptr::read(data as *const [f32; 512 * 424])
        };

        for spin in env.base_env.cell_lattice.iter_values_mut() {
            if matches!(spin, Spin::Solid) {
                *spin = Spin::Medium;
            }
        }

        let mut count = 0;
        for j in 0..424 {
            for i in 0..512 {
                let index = Pos::new(j, i).col_major(512);
                let depth = data_arr[index];
                if depth < self.min_depth || depth > self.max_depth {
                    continue;
                }

                env.grant_position(Pos::new(i, 423 - j), Spin::Solid);
                count += 1;
            }
        }
        unsafe { kinect_release_frame(self.handle) };
        Ok(count)
    }
}

impl Drop for KinectListener {
    fn drop(&mut self) {
        unsafe { kinect_destroy(self.handle) };
    }
}