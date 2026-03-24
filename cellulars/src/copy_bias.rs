//! Contains different copy biases applied during a [`Potts`](crate::potts::Potts) step.

/*
 TODO!: PerimeterBias, ActBias
    For PerBias:
        Best way I can think to implement this would be through a new trait TrackPerimeter which can be implemented
        for cells. Methods are perimeter() and possibly delta_perimeter() -> Option<i32>.
        This last method is needed if we care to implement the optimization where the cells neighbors dont need to
        be iterated ove both in copy_bias() and in transfer_position(). Maybe there is a way to prevent this double
        iterations for other operations besides calculating the perimeter? Like storing a Vec<Pos<usize>> in cells
        or something of sorts (although a Vec is most def too slow).
 */

use crate::constants::FloatType;
use crate::lattice::Lattice;
use crate::prelude::{Boundary, Pos, Spin};

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
    pub lambda: FloatType,
    pub dir_params: DirectionalOptions
}

impl CopyBias<ChemContext<'_>> for ChemotaxisBias {
    fn bias(&self, pos_source: Pos<usize>, pos_target: Pos<usize>, context: &ChemContext) -> FloatType {
        directional_bias(
            context.cell_lattice[pos_source],
            context.cell_lattice[pos_target],
            -self.lambda * (context.chem_lattice[pos_target] - context.chem_lattice[pos_source]),
            &self.dir_params
        )
    }
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ChemContext<'a> {
    pub cell_lattice: &'a Lattice<Spin>,
    pub chem_lattice: &'a Lattice<FloatType>,
}

#[derive(Debug, Clone, Copy, PartialEq)]
/// Bias copies towards the source of a chemical stored in a [`Lattice<FloatType>`].
pub struct DirectionBias<B> {
    /// Strength of the chemotaxis constraint on the energy functional.
    pub lambda: FloatType,
    /// Boundary conditions used to evaluate directionality.
    pub boundary: B,
    pub dir_params: DirectionalOptions
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DirectionContext<'a> {
    pub cell_lattice: &'a Lattice<Spin>,
    pub angle: FloatType,
}

impl<B: Boundary<Coord = FloatType>> DirectionBias<B> {
    pub fn angle_from_positions(&self, cell: Pos<FloatType>, target: Pos<FloatType>) -> FloatType {
        let (dx, dy) = self.boundary.displacement(cell, target);
        dy.atan2(dx)
    }
}

impl<B: Boundary<Coord = FloatType>> CopyBias<DirectionContext<'_>> for DirectionBias<B> {
    fn bias(&self, pos_source: Pos<usize>, pos_target: Pos<usize>, context: &DirectionContext<'_>) -> FloatType {
        let angle_pos = self.angle_from_positions(pos_source.cast_as(), pos_target.cast_as());
        directional_bias(
            context.cell_lattice[pos_source],
            context.cell_lattice[pos_target],
            -self.lambda * (context.angle - angle_pos).cos(),
            &self.dir_params
        )
    }
}

/// Parameters determining how to handle directional biases. 
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DirectionalOptions {
    /// Whether cell protrusions count towards the energy functional.
    pub protrusions: bool,
    /// Whether cell retractions count towards the energy functional.
    pub retractions: bool,
    /// Whether cell protrusions/retractions should be discarded if the target position is not the medium.
    pub contact_inhibition: bool,
}

fn directional_bias(
    spin_source: Spin,
    spin_target: Spin,
    energy_diff: FloatType,
    dir_params: &DirectionalOptions
) -> FloatType {
    if dir_params.contact_inhibition
        && !matches!(spin_source, Spin::Medium)
        && !matches!(spin_target, Spin::Medium) {
        return 0.;
    }
    let mut energy = 0.;
    if dir_params.protrusions && matches!(spin_source, Spin::Some(_)) {
        energy += energy_diff;
    }
    if dir_params.retractions && matches!(spin_target, Spin::Some(_)) {
        energy += energy_diff;
    }
    energy
}