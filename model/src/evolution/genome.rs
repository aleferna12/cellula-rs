use rand::Rng;

pub trait Genome {
    fn attempt_mutate(&mut self, rng: &mut impl Rng) -> u32;
}