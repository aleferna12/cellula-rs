use crate::constants::KinectNeighbourhoodType;
use crate::kinect::{KinectHandle, kinect_color_depth, kinect_create, kinect_destroy, kinect_listen_frame, kinect_release_frame};
use crate::my_environment::MyEnvironment;
use cellulars_lib::positional::boundaries::Boundary;
use cellulars_lib::prelude::{Habitable, Neighbourhood, Pos};
use cellulars_lib::spin::Spin;
use image::{Rgba, RgbaImage};
use image::imageops::flip_vertical;

pub struct KinectListener {
    handle: *mut KinectHandle,
    pub min_depth: f32,
    pub max_depth: f32,
    pub frame_period: u32,
    pub neighbourhood: KinectNeighbourhoodType
}

impl KinectListener {
    pub fn new(
        min_depth: f32,
        max_depth: f32,
        frame_period: u32,
        neighbourhood: KinectNeighbourhoodType
    ) -> Option<Self> {
        let handle = unsafe { kinect_create(true, true) };
        if handle.is_null() {
            return None;
        }
        Some(Self {
            handle,
            min_depth,
            max_depth,
            frame_period,
            neighbourhood
        })
    }

    pub fn draw_silhouette(&mut self, env: &mut MyEnvironment) -> anyhow::Result<RgbaImage> {
        if env.env.width() != 512 || env.env.height() != 424 {
            anyhow::bail!("kinect can only be used if pond's width is 512 and height is 424");
        }

        // Experimented with this being async but the thread spawn cost is not worth it
        // unless we rewrite the C part to also be async, which would prob be painful due to FFI
        let (depth_arr, color_arr) = unsafe { Self::fetch_depth_color(self.handle)? };
        let img = flip_vertical(&RgbaImage::from_raw(
            512,
            424,
            color_arr.to_vec()
        ).ok_or(anyhow::anyhow!("failed to create color image"))?);
        let mut frame = RgbaImage::new(512, 424);

        for j in 0..424 {
            for i in 0..512 {
                let pos = Pos::new(i as usize, j as usize);
                if !self.should_display(pos, depth_arr, &env.env.bounds.lattice_boundary) {
                    continue;
                }
                env.grant_position(pos, Spin::Solid);
                let pixel = img.get_pixel(i, j);
                frame.put_pixel(i, j, Rgba([pixel.0[2], pixel.0[1], pixel.0[0], 255]));
            }
        }
        unsafe { kinect_release_frame(self.handle) };
        Ok(frame)
    }

    fn should_display(
        &self,
        pos: Pos<usize>,
        data_arr: &[f32],
        lattice_boundary: &impl Boundary<Coord = isize>
    ) -> bool {
        if !self.in_bounds(data_arr[Self::flat_depth_index(pos)]) {
            return false;
        }
        for neigh in self.neighbourhood.neighbours(pos.cast_as()) {
            let Some(valid_neigh) = lattice_boundary.valid_pos(neigh) else {
                continue;
            };
            if !self.in_bounds(data_arr[Self::flat_depth_index(valid_neigh.cast_as())]) {
                return false;
            }
        }
        true
    }

    fn in_bounds(&self, depth: f32) -> bool {
        depth > self.min_depth && depth < self.max_depth
    }

    fn flat_depth_index(pos: Pos<usize>) -> usize {
        Pos::new(423 - pos.y, pos.x).col_major(512)
    }

    // Unsafe because the slice is not tied to any particular lifetime 'u
    unsafe fn fetch_depth_color<'u>(handle: *mut KinectHandle) -> anyhow::Result<(&'u [f32], &'u [u8])> {
        if handle.is_null() { anyhow::bail!("kinect handle was lost") };

        unsafe {
            if !kinect_listen_frame(handle, 10_000) {
                anyhow::bail!("failed to get frame");
            }
            let color_depth = kinect_color_depth(handle);
            if color_depth.depth.is_null() {
                anyhow::bail!("failed to fetch next depth");
            }
            let depth_slice = std::slice::from_raw_parts(color_depth.depth, 512 * 424);

            if color_depth.color.is_null() {
                anyhow::bail!("failed to fetch next color");
            }
            let color_slice = std::slice::from_raw_parts(color_depth.color, 512 * 424 * 4);
            Ok((depth_slice, color_slice))
        }
    }
}

impl Drop for KinectListener {
    fn drop(&mut self) {
        unsafe { kinect_destroy(self.handle) };
    }
}