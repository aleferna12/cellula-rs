use crate::evolution::genome::Genome;
use bitvec::BitArr;
use bitvec::slice::BitSlice;
use rand::Rng;
use crate::constants::MAX_GENOME_LENGTH;

type ArrayType = BitArr!(for MAX_GENOME_LENGTH as usize, in u64);

#[derive(Debug, Clone)]
pub struct BitGenome {
    ligands: ArrayType,
    receptors: ArrayType,
    pub mut_rate: f32,
    pub length: u16
}

impl BitGenome {
    pub fn new_random(mut_rate: f32, length: u16, rng: &mut impl Rng) -> Self {
        Self {
            ligands: ArrayType::new(rng.random()),
            receptors: ArrayType::new(rng.random()),
            mut_rate,
            length
        }
    }

    pub fn new_empty(mut_rate: f32, length: u16) -> Self {
        Self {
            ligands: ArrayType::ZERO,
            receptors: ArrayType::ZERO,
            mut_rate,
            length
        }
    }
    
    pub fn from_iterators(
        ligands: impl IntoIterator<Item = bool>,
        receptors: impl IntoIterator<Item = bool>,
        mut_rate: f32,
        length: u16
    ) -> Option<Self> {
        let mut ligands = ligands.into_iter();
        let mut receptors = receptors.into_iter();
        let mut ligands_array = ArrayType::ZERO;
        let mut receptor_array = ArrayType::ZERO;

        for i in 0..length as usize{
            ligands_array.set(i, ligands.next()?);
            receptor_array.set(i, receptors.next()?);
        }

        Some(Self {
            ligands: ligands_array,
            receptors: receptor_array,
            mut_rate,
            length
        })
    }

    pub fn ligands(&self) -> &BitSlice<u64> {
        &self.ligands.as_bitslice()[..self.length as usize]
    }

    pub fn receptors(&self) -> &BitSlice<u64> {
        &self.receptors.as_bitslice()[..self.length as usize]
    }

    pub fn ligands_string(&self) -> String {
        self.ligands()
            .iter()
            .map(|b| if *b { "1" } else { "0" })
            .collect::<String>()
    }

    pub fn receptors_string(&self) -> String {
        self.receptors()
            .iter()
            .map(|b| if *b { "1" } else { "0" })
            .collect::<String>()
    }
}

impl Genome for BitGenome {
    fn attempt_mutate(&mut self, rng: &mut impl Rng) -> u32 {
        let mut mut_count = 0;
        for i in 0..self.length as usize {
            if rng.random::<f32>() < self.mut_rate {
                let prev = self.ligands[i];
                self.ligands.set(i, !prev);
                mut_count += 1;
            }
            if rng.random::<f32>() < self.mut_rate {
                let prev = self.receptors[i];
                self.receptors.set(i, !prev);
                mut_count += 1;
            }
        }
        mut_count
    }
}