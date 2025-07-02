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
        self.area += 1;
        self.center = Pos2D::new(
            self.center.x + (pos.x as f32 - self.center.x) / self.area as f32,
            self.center.y + (pos.y as f32 - self.center.y) / self.area as f32,
        );
    }

    pub fn remove_position(&mut self, pos: Pos2D<usize>) {
        self.area -= 1;
        self.center = Pos2D::new(
            self.center.x - (pos.x as f32 - self.center.x) / self.area as f32,
            self.center.y - (pos.y as f32 - self.center.y) / self.area as f32,
        );
    }
}