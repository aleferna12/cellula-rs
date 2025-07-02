use std::f32::consts::PI;
use crate::pos::Pos2D;

pub type Sigma = i32;

#[derive(Debug)]
pub struct Cell {
    pub sigma: Sigma,
    pub area: u32,
    pub target_area: u32,
    pub growth_timer: u32,
    pub center: Pos2D<f32>
}

impl Cell {
    pub fn new(sigma: Sigma, area: u32, target_area: u32, center: Pos2D<f32>) -> Self {
        Self {
            sigma,
            area,
            target_area,
            center,
            growth_timer: 0,
        }
    }
    
    pub fn add_position(&mut self, pos: Pos2D<usize>) {
        // TODO! fix
        //   should we even keep track of this? It costs 16% of performance
        //   maybe just calculate on site...

        let two_pi = 2.0 * PI;

        // Convert CoM to angles
        let theta_old_x = two_pi * self.center.x / 128f32;
        let theta_old_y = two_pi * self.center.y / 128f32;

        let sum_cos_x = self.area as f32 * theta_old_x.cos() + (2.0 * PI * (pos.x as f32) / 128f32).cos();
        let sum_sin_x = self.area as f32 * theta_old_x.sin() + (2.0 * PI * (pos.x as f32) / 128f32).sin();

        let sum_cos_y = self.area as f32 * theta_old_y.cos() + (2.0 * PI * (pos.y as f32) / 128f32).cos();
        let sum_sin_y = self.area as f32 * theta_old_y.sin() + (2.0 * PI * (pos.y as f32) / 128f32).sin();

        let new_theta_x = sum_sin_x.atan2(sum_cos_x);
        let new_theta_y = sum_sin_y.atan2(sum_cos_y);

        let new_cx = 128f32 * new_theta_x.rem_euclid(two_pi) / two_pi;
        let new_cy = 128f32 * new_theta_y.rem_euclid(two_pi) / two_pi;
        
        self.center = Pos2D::new(new_cx, new_cy);
        self.area += 1;
    }

    pub fn remove_position(&mut self, pos: Pos2D<usize>) {
        let two_pi = 2.0 * PI;

        // Convert CoM to angles
        let theta_old_x = two_pi * self.center.x / 128f32;
        let theta_old_y = two_pi * self.center.y / 128f32;

        let sum_cos_x = self.area as f32 * theta_old_x.cos() - (2.0 * PI * (pos.x as f32) / 128f32).cos();
        let sum_sin_x = self.area as f32 * theta_old_x.sin() - (2.0 * PI * (pos.x as f32) / 128f32).sin();

        let sum_cos_y = self.area as f32 * theta_old_y.cos() - (2.0 * PI * (pos.y as f32) / 128f32).cos();
        let sum_sin_y = self.area as f32 * theta_old_y.sin() - (2.0 * PI * (pos.y as f32) / 128f32).sin();

        let new_theta_x = sum_sin_x.atan2(sum_cos_x);
        let new_theta_y = sum_sin_y.atan2(sum_cos_y);

        let new_cx = 128f32 * new_theta_x.rem_euclid(two_pi) / two_pi;
        let new_cy = 128f32 * new_theta_y.rem_euclid(two_pi) / two_pi;

        self.center = Pos2D::new(new_cx, new_cy);
        self.area -= 1;
    }
}