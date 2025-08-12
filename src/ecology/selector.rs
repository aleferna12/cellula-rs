use crate::cell::Fit;
use rand::Rng;

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
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256StarStar;
    use std::collections::HashSet;

    #[derive(Debug, Hash, PartialEq, Eq, Copy, Clone)]
    struct FitTest(usize);

    impl Fit for FitTest {
        fn fitness(&self) -> f32 {
            self.0 as f32
        }
    }
    
    #[test]
    fn test_weighted_selection() {
        let fits: [FitTest; 100] = core::array::from_fn(|i| FitTest(i));
        let fit_sum = fits.iter().map(|f| f.fitness()).sum::<f32>();
        assert_eq!(fit_sum, 4950.);

        for seed in 0..100 {
            let mut sel = WeightedSelection {
                select_n: 100,
                rng: Xoshiro256StarStar::seed_from_u64(seed)
            };
            let new_fits = sel.select(&fits);
            let new_fit_sum = new_fits.iter().map(|f| f.fitness()).sum::<f32>();
            assert!(new_fit_sum > fit_sum);

            let unique: HashSet<FitTest> = new_fits.into_iter().copied().collect();
            assert!(unique.len() < sel.select_n as usize);
        }
    }
}