//! Contains different copy biases applied during a [`Potts`](crate::potts::Potts) step.

// TODO!: PerimeterBias, ActBias

use crate::constants::FloatType;
use crate::lattice::Lattice;
use crate::prelude::{Boundary, Pos};

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

#[derive(Debug, Clone, Copy, PartialEq)]
/// Bias copies towards the source of a chemical stored in a [`Lattice<FloatType>`].
pub struct DirectionBias<B> {
    /// Strength of the chemotaxis constraint on the energy functional.
    pub lambda: FloatType,
    pub boundary: B
}

impl<B: Boundary<Coord = FloatType>> DirectionBias<B> {
    pub fn angle_from_positions(&self, cell: Pos<FloatType>, target: Pos<FloatType>) -> FloatType {
        let (dx, dy) = self.boundary.displacement(cell, target);
        dy.atan2(dx)
    }
}

impl<B: Boundary<Coord = FloatType>> CopyBias<FloatType> for DirectionBias<B> {
    fn bias(&self, pos_source: Pos<usize>, pos_target: Pos<usize>, angle: &FloatType) -> FloatType {
        let angle_pos = self.angle_from_positions(pos_source.cast_as(), pos_target.cast_as());
        -self.lambda * (angle - angle_pos).cos()
    }
}