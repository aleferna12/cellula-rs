pub trait Step {
    fn step(&mut self);

    fn run_for(&mut self, steps: u32) {
        for _ in 0..steps {
            self.step();
        }
    }
}