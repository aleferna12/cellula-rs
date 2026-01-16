use cellulars_lib::positional::boundaries::Boundary;
use crate::environment::Environment;
use crate::kinect::{kinect_create, kinect_destroy, kinect_next_depth, kinect_release_frame, KinectHandle};
use cellulars_lib::prelude::{Habitable, Neighbourhood, Pos};
use cellulars_lib::spin::Spin;
use crate::constants::KinectNeighbourhoodType;

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
        let handle = unsafe { kinect_create() };
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

    pub fn draw_silhouette(&mut self, env: &mut Environment) -> anyhow::Result<u32> {
        if env.base_env.width() != 512 || env.base_env.height() != 424 {
            anyhow::bail!("kinect can only be used if pond's width is 512 and height is 424");
        }

        for spin in env.base_env.cell_lattice.iter_values_mut() {
            if matches!(spin, Spin::Solid) {
                *spin = Spin::Medium;
            }
        }

        // Experimented with this being async but the thread spawn cost is not worth it
        // unless we rewrite the C part to also be async, which would prob be painful due to FFI
        let data_arr = unsafe { Self::fetch_depth(self.handle)? };
        let mut count = 0;
        for j in 0..424 {
            for i in 0..512 {
                let pos = Pos::new(i, j);
                if !self.should_display(pos, data_arr, &env.base_env.bounds.lattice_boundary) {
                    continue;
                }
                env.grant_position(pos, Spin::Solid);
                count += 1;
            }
        }
        unsafe { kinect_release_frame(self.handle) };
        Ok(count)
    }

    fn should_display(
        &self,
        pos: Pos<usize>,
        data_arr: &[f32],
        lattice_boundary: &impl Boundary<Coord = isize>
    ) -> bool {
        if !self.in_bounds(data_arr[Self::flat_frame_index(pos)]) {
            return false;
        }
        for neigh in self.neighbourhood.neighbours(pos.cast_as()) {
            let Some(valid_neigh) = lattice_boundary.valid_pos(neigh) else {
                continue;
            };
            if !self.in_bounds(data_arr[Self::flat_frame_index(valid_neigh.cast_as())]) {
                return false;
            }
        }
        true
    }

    fn in_bounds(&self, depth: f32) -> bool {
        depth > self.min_depth && depth < self.max_depth
    }

    fn flat_frame_index(pos: Pos<usize>) -> usize {
        Pos::new(423 - pos.y, pos.x).col_major(512)
    }

    // Unsafe because the slice is not tied to any particular lifetime 'u
    unsafe fn fetch_depth<'u>(handle: *mut KinectHandle) -> anyhow::Result<&'u [f32]> {
        if handle.is_null() { anyhow::bail!("kinect handle was lost") };

        unsafe {
            let data = kinect_next_depth(handle, 10_000);
            if data.is_null() {
                anyhow::bail!("failed to fetch next depth");
            }
            Ok(std::slice::from_raw_parts(data, 512 * 424))
        }
    }
}

impl Drop for KinectListener {
    fn drop(&mut self) {
        unsafe { kinect_destroy(self.handle) };
    }
}