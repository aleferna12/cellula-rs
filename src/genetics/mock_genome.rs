use crate::genetics::genome::Genome;
use rand::Rng;

/// This is a fake genome that just cycles through a boolean `state`.
#[derive(Clone, Debug)]
pub struct MockGenome {
    period_updates: u32,
    counter: u32,
    pub state: bool
}

impl MockGenome {
    /// Makes a new `MockGenome` with a specified period.
    ///
    /// `period_updates` is the period for which each cell type will last for.
    /// The unit is the number of `update_expression()` calls, not MCS.
    pub fn new(period_updates: u32) -> Self {
        Self {
            period_updates,
            counter: 0,
            state: false,
        }
    }
}

impl Genome for MockGenome {
    fn attempt_mutate(&mut self, _rng: &mut impl Rng) -> bool {
        false
    }

    fn update_expression(&mut self) {
        self.counter += 1;
        if self.counter > self.period_updates {
            self.state = !self.state;
            self.counter = 0;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cycles_state_after_period() {
        let mut g = MockGenome::new(2);
        assert_eq!(g.state, false);

        g.update_expression(); // counter=1
        assert_eq!(g.state, false);

        g.update_expression(); // counter=2
        assert_eq!(g.state, false);

        g.update_expression(); // counter resets, state flips
        assert_eq!(g.state, true);
        assert_eq!(g.counter, 0);
    }
}