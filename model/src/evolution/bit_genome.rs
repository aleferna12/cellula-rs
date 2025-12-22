use crate::evolution::genome::Genome;
use rand::Rng;

#[derive(Debug, Clone)]
pub struct BitGenome {
    ligands: u64,
    receptors: u64,
    pub mut_rate: f32,
    pub length: u8
}

impl BitGenome {
    /// Creates a new [BitGenome] using `length` significant bits of `ligands` and `receptors` if 0 < `length` < 64.
    pub fn new(ligands: u64, receptors: u64, mut_rate: f32, length: u8) -> Option<Self> {
        if length == 0 || length > 64 {
            return None;
        }
        Some(Self {
            ligands: Self::truncate(ligands, length),
            receptors: Self::truncate(receptors, length),
            mut_rate,
            length
        })
    }

    /// Creates a new [BitGenome] with random ligands and receptors if 0 < `length` < 64.
    pub fn new_random(mut_rate: f32, length: u8, rng: &mut impl Rng) -> Option<Self> {
        Self::new(
            rng.random::<u64>(),
            rng.random::<u64>(),
            mut_rate,
            length
        )
    }

    pub fn ligands(&self) -> u64 {
        self.ligands
    }

    pub fn receptors(&self) -> u64 {
        self.receptors
    }

    fn truncate(protein: u64, length: u8) -> u64 {
        protein & (u64::MAX >> (64 - length))
    }

    fn flip_bit(protein: u64, bit_index: u8) -> u64 {
        protein ^ (1 << bit_index)
    }
}

impl Genome for BitGenome {
    fn attempt_mutate(&mut self, rng: &mut impl Rng) -> u32 {
        let mut mut_count = 0;
        for i in 0..self.length {
            if rng.random::<f32>() < self.mut_rate {
                self.ligands = Self::flip_bit(self.ligands, i);
                mut_count += 1;
            }
            if rng.random::<f32>() < self.mut_rate {
                self.receptors = Self::flip_bit(self.receptors, i);
                mut_count += 1;
            }
        }
        mut_count
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_bit_genome() {
        assert!(BitGenome::new(0, 0, 0., 0).is_none());
        assert!(BitGenome::new(0, 0, 0., 65).is_none());

        let bit_genome = BitGenome::new(
            255,
            257,
            0.,
            8
        ).unwrap();
        assert_eq!(bit_genome.ligands(), 255);
        assert_eq!(bit_genome.receptors(), 1);
    }

    #[test]
    fn test_flip_bit() {
        let ligand = 0b1101u64;
        assert_eq!(0b1100, BitGenome::flip_bit(ligand, 0));
        assert_eq!(0b1111, BitGenome::flip_bit(ligand, 1));
    }

    #[test]
    fn test_mut() {
        let mut rng = rand::rng();
        let mut bit_genome = BitGenome::new_random(0.1, 4, &mut rng).unwrap();
        let mutated = (0..100000).reduce(|i, _| i + bit_genome.attempt_mutate(&mut rng)).unwrap();
        // 10% acceptable error
        assert!(mutated > 9000 * 8 && mutated < 11000 * 8);
    }
}