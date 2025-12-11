#![allow(dead_code)]  // Want to keep the WightedOrder stuff for now

use rand::Rng;

pub trait Selector {
    fn select<'a, F: Fit>(&mut self, selectable: &'a [F]) -> Vec<&'a F>;
}

pub struct WeightedSelection<'r, R> {
    pub select_n: u32,
    pub rng: &'r mut R
}

impl<R: Rng> WeightedSelection<'_, R> {
    fn select_random_fit<'a, F: Fit>(&mut self, selectable: &'a [F], tot_fit: f32) -> &'a F {
        let mut rand_fit = self.rng.random::<f32>() * tot_fit;
        for s in selectable.iter() {
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

impl<R: Rng> Selector for WeightedSelection<'_, R> {
    fn select<'a, F: Fit>(&mut self, selectable: &'a [F]) -> Vec<&'a F> {
        let tot_fit = selectable
            .iter()
            .map(|s| { s.fitness() })
            .sum();
        (0..self.select_n)
            .map(|_| self.select_random_fit(selectable, tot_fit))
            .collect()
    }
}

pub struct WeightedOrderedSelection<'r, R> {
    pub rng: &'r mut R
}

impl<R: Rng> Selector for WeightedOrderedSelection<'_, R> {
    fn select<'a, F: Fit>(&mut self, selectable: &'a [F]) -> Vec<&'a F> {
        let mut selector = WeightedSelection {
            select_n: selectable.len() as u32,
            rng: &mut self.rng
        };
        let indexed = selectable.iter()
            .enumerate()
            .map(|(i, f)| IndexedFit {
                index: i,
                fit: f
            })
            .collect::<Vec<_>>();
        let selected: Vec<_> = selector.select(&indexed);

        // TODO: it might be faster to iterate selected while removing elements, but then for each
        //  `s` in selectable we need to iterate selected once while comparing elements which might be slow
        //  Benchmark
        let mut selection_count = vec![0u32; selected.len()];
        for &s in &selected {
            selection_count[s.index] += 1
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

        parents.iter().map(|&i| &selectable[i]).collect()
    }
}

pub trait Fit {
    fn fitness(&self) -> f32;
}

pub struct IndexedFit<'f, F> {
    pub fit: &'f F,
    pub index: usize,
}

impl<F: Fit> Fit for IndexedFit<'_, F> {
    fn fitness(&self) -> f32 {
        self.fit.fitness()
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
        let fits: [FitTest; 100] = core::array::from_fn(FitTest);
        let fit_sum = fits.iter().map(|f| f.fitness()).sum::<f32>();
        assert_eq!(fit_sum, 4950.);

        for seed in 0..100 {
            let mut sel = WeightedSelection {
                select_n: 100,
                rng: &mut Xoshiro256StarStar::seed_from_u64(seed)
            };
            let selected = sel.select(&fits);
            let selected_sum = selected.iter().map(|&f| f.fitness()).sum::<f32>();
            assert!(selected_sum > fit_sum);

            let unique: HashSet<_> = selected.iter().collect();
            assert!(unique.len() < selected.len());
        }
    }

    #[test]
    fn test_weighted_ordered_selection() {
        let fits: [FitTest; 100] = core::array::from_fn(FitTest);
        let fit_sum = fits.iter().map(|f| f.fitness()).sum::<f32>();
        assert_eq!(fit_sum, 4950.);

        for seed in 0..100 {
            let mut order_sel = WeightedOrderedSelection {
                rng: &mut Xoshiro256StarStar::seed_from_u64(seed)
            };
            let mut sel = WeightedSelection {
                select_n: 100,
                rng: &mut Xoshiro256StarStar::seed_from_u64(seed)
            };

            let ordered_selected = order_sel.select(&fits);
            let ordered_selected_sum = ordered_selected.iter().map(|&f| f.fitness()).sum::<f32>();
            let selected = sel.select(&fits);
            let selected_sum = selected.iter().map(|&f| f.fitness()).sum::<f32>();
            // Same organisms were selected
            assert_eq!(ordered_selected_sum, selected_sum);

            for f in ordered_selected.iter() {
                let FitTest(index) = f;
                // All survivors preserved their position
                assert_eq!(f, &ordered_selected[*index]);
            }

            // Some organisms reproduced
            let set: HashSet<&FitTest> = HashSet::from_iter(ordered_selected.iter().copied());
            assert!(ordered_selected.len() > set.len())
        }
    }
}