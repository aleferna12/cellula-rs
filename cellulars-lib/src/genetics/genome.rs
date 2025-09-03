use rand::Rng;

pub trait Genome {
    fn attempt_mutate(&mut self, rng: &mut impl Rng) -> bool;
    fn update_expression(&mut self);
}