//! Contains different copy biases applied during a [`Potts`](crate::potts::Potts) step.

// TODO!: PerimeterBias, ActBias

use crate::constants::FloatType;
use crate::lattice::Lattice;
use crate::prelude::Pos;

/// Defines a bias in the energy functional.
pub trait CopyBias<C> {
    /// Computes the bias in the Hamiltonian energy functional from the two positions and a given `context`.
    fn bias(&self, pos_source: Pos<usize>, pos_target: Pos<usize>, context: &C) -> FloatType;
}

/// Computes no biases besides the size and adhesion terms.
///
/// [`NoBias::bias()`] returns 0.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct NoBias;

impl<C> CopyBias<C> for NoBias {
    fn bias(&self, _pos_source: Pos<usize>, _pos_target: Pos<usize>, _context: &C) -> FloatType {
        0.
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
/// Bias copies towards the source of a chemical stored in a [`Lattice<FloatType>`].
pub struct ChemotaxisBias {
    /// Strength of the chemotaxis constraint on the energy functional.
    pub lambda: FloatType
}

impl CopyBias<Lattice<FloatType>> for ChemotaxisBias {
    fn bias(&self, pos_source: Pos<usize>, pos_target: Pos<usize>, chem_lattice: &Lattice<FloatType>) -> FloatType {
        -self.lambda * (chem_lattice[pos_target] - chem_lattice[pos_source])
    }
}