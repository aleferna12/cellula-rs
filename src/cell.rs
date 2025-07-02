pub type Sigma = i32;

#[derive(Debug)]
pub struct Cell {
    pub sigma: Sigma,
    pub area: u32,
    pub target_area: u32,
    pub growth_timer: u32
}

impl Cell {
    pub fn new(sigma: Sigma, area: u32, target_area: u32) -> Self {
        Self {
            sigma,
            area,
            target_area,
            growth_timer: 0
        }
    }
}