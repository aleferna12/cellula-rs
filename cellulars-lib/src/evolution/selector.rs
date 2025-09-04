use rand::Rng;

pub trait Selector {
    fn select<F: Fit>(&mut self, selectable: &[F]) -> Vec<usize>;
}

/// This trait is an assurance that the output of `select` will have the same number of elements as the input
/// and that all organisms that didn't die will be represented in the same place.
pub trait PreservesOrder: Selector {}

pub struct WeightedSelection<R> {
    select_n: u32,
    rng: R
}

impl<R: Rng> WeightedSelection<R> {
    fn select_random_fit<F: Fit>(&mut self, selectable: &[F], tot_fit: f32) -> usize {
        let mut rand_fit = self.rng.random::<f32>() * tot_fit;
        for (i, s) in selectable.iter().enumerate() {
            let this_fit = s.fitness();
            if rand_fit < this_fit {
                return i
            } else {
                rand_fit -= this_fit;
            }
        }
        self.rng.random_range(0..selectable.len())
    }
}

impl<R: Rng> Selector for WeightedSelection<R> {
    fn select<F: Fit>(&mut self, selectable: &[F]) -> Vec<usize> {
        let tot_fit = selectable
            .iter()
            .map(|s| { s.fitness() })
            .sum();
        (0..self.select_n)
            .map(|_| self.select_random_fit(selectable, tot_fit))
            .collect()
    }
}

pub struct WeightedOrderedSelection<R> {
    pub rng: R
}

impl<R: Rng> Selector for WeightedOrderedSelection<R> {
    fn select<F: Fit>(&mut self, selectable: &[F]) -> Vec<usize> {
        let mut selector = WeightedSelection {
            select_n: selectable.len() as u32,
            rng: &mut self.rng
        };
        let selected: Vec<_> = selector.select(selectable);
        // TODO: it might be faster to iterate selected while removing elements, but then for each
        //  `s` in selectable we need to iterate selected once while comparing elements which might be slow
        //  Benchmark
        let mut selection_count = vec![0u32; selected.len()];
        for &s in &selected {
            selection_count[s] += 1
        }
        let mut dead = vec![];
        let mut parents = vec![0; selected.len()];
        let mut offspring = vec![0; selected.len()];
        for (i, count) in selection_count.into_iter().enumerate() {
            if count == 0 {
                dead.push(i);
                continue;
            }
            offspring[i] = count - 1;
            parents[i] = i;
        }

        for i in dead {
            let parent = offspring.iter().position(|&count| count > 0).unwrap();
            parents[i] = parent;
            offspring[parent] -= 1;
        }

        parents
    }
}

impl<R: Rng> PreservesOrder for WeightedOrderedSelection<R> {}

pub trait Fit {
    fn fitness(&self) -> f32;
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
            let new_fit_sum = new_fits.iter().map(|&f| fits[f].fitness()).sum::<f32>();
            assert!(new_fit_sum > fit_sum);

            let unique: HashSet<_> = new_fits.iter().collect();
            assert!(unique.len() < new_fits.len());
        }
    }

    #[test]
    fn test_weighted_ordered_selection() {
        let fits: [FitTest; 100] = core::array::from_fn(|i| FitTest(i));
        let fit_sum = fits.iter().map(|f| f.fitness()).sum::<f32>();
        assert_eq!(fit_sum, 4950.);

        for seed in 0..100 {
            let mut sel = WeightedOrderedSelection {
                rng: Xoshiro256StarStar::seed_from_u64(seed)
            };
            let mut sel_ = WeightedSelection {
                select_n: 100,
                rng: Xoshiro256StarStar::seed_from_u64(seed)
            };

            let new_fits = sel.select(&fits);
            let new_fit_sum = new_fits.iter().map(|&f| fits[f].fitness()).sum::<f32>();
            let new_fits_ = sel_.select(&fits);
            let new_fit_sum_ = new_fits_.iter().map(|&f| fits[f].fitness()).sum::<f32>();
            // Same organisms were selected
            assert_eq!(new_fit_sum, new_fit_sum_);

            for &f in &new_fits {
                // All survivors preserved their position
                assert_eq!(f, new_fits[f]);
            }
        }
    }
}