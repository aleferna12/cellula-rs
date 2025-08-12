use rand::Rng;
use crate::cell::Fit;

pub trait Selector {
    // This can instead return an impl Iterator<Item = F> if needed for generalisation
    fn select<'f, F: Fit>(&mut self, selectable: &'f [F]) -> Vec<&'f F>;
}

pub struct WeightedSelection<R> {
    select_n: u32,
    rng: R
}

impl<R: Rng> WeightedSelection<R> {
    fn select_random_fit<'f, F: Fit>(&mut self, selectable: &'f [F], tot_fit: f32) -> &'f F {
        let mut rand_fit = self.rng.random::<f32>() * tot_fit;
        for s in selectable {
            let this_fit = s.fitness();
            if rand_fit < this_fit {
                return s
            } else {
                rand_fit -= this_fit;
            }
        }
        &selectable[self.rng.random_range(0..selectable.len())]
    }
}

impl<R: Rng> Selector for WeightedSelection<R> {
    fn select<'f, F: Fit>(&mut self, selectable: &'f [F]) -> Vec<&'f F> {
        let tot_fit = selectable
            .iter()
            .map(|s| { s.fitness() })
            .sum();
        let mut selected = vec![];
        for _ in 0..self.select_n {
            selected.push(self.select_random_fit(selectable, tot_fit));
        }
        selected
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_weighted_selection() {
        
    }
}