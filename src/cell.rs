#[derive(Debug)]
pub struct Cell {
    pub area: u32,
    pub target_area: u32
}
impl Cell {
    pub fn new(area: u32, target_area: u32) -> Cell {
        Self { area, target_area }
    }
}