//! Contains logic associated with [Step].

/// This trait describes objects that can run a simulation by successively executing [Step::step()].
pub trait Step {
    /// Execute one step of the simulation.
    fn step(&mut self);

    /// Runs the simulation for a specific number of `steps`.
    fn run_for(&mut self, steps: u32) {
        for _ in 0..steps {
            self.step();
        }
    }
}