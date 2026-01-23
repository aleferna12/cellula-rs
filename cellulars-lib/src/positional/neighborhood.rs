//! Contains logic associated with neighborhoods for a discrete lattice.

use crate::positional::pos::Pos;

const MAX_NEIGH_R: u8 = 16;
const MOORE_SIZE: usize = 4 * MAX_NEIGH_R as usize * (MAX_NEIGH_R as usize + 1);
const MOORE_NEIGHS: [(i16, i16); MOORE_SIZE] = {
    let mut ret = [(0i16, 0i16); MOORE_SIZE];
    let mut r = 1;
    let mut flat_index = 0usize;
    while r <= MAX_NEIGH_R as i16 {
        let mut i = -r;
        while i <= r {
            let mut j = -r;
            while j <= r {
                let max_abs = if i.abs() > j.abs() { i.abs() } else { j.abs() };
                if max_abs == r {
                    ret[flat_index] = (i, j);
                    flat_index += 1;
                }
                j += 1;
            }
            i += 1;
        }
        r += 1;
    }
    ret
};

const VON_NEUMANN_SIZE: usize = 2 * (MAX_NEIGH_R as usize) * (MAX_NEIGH_R as usize + 1);
const VON_NEUMANN_NEIGHS: [(i16, i16); VON_NEUMANN_SIZE] = {
    let mut ret = [(0i16, 0i16); VON_NEUMANN_SIZE];
    let mut flat_index = 0usize;

    let mut r = 1;
    while r <= MAX_NEIGH_R as i16 {
        let mut dx = -r;
        while dx <= r {
            let dy_abs = r - dx.abs();

            if dy_abs == 0 {
                ret[flat_index] = (dx, 0);
                flat_index += 1;
            } else {
                ret[flat_index] = (dx, -dy_abs);
                flat_index += 1;
                ret[flat_index] = (dx, dy_abs);
                flat_index += 1;
            }

            dx += 1;
        }
        r += 1;
    }

    ret
};

// NEVER REMOVE THIS INLINE
#[inline(always)]
fn fetch_neighs(
    pos: Pos<isize>,
    neigh_array: &[(i16, i16)],
    n_neighs: u16
) -> impl Iterator<Item = Pos<isize>> {
    neigh_array[..n_neighs.into()]
        .iter()
        .map(move |(i, j)| {
            Pos::new(
                pos.x + *i as isize,
                pos.y + *j as isize,
            )
        })
}

/// Describes a neighborhood of a square, discrete lattice.
pub trait Neighborhood {
    /// Returns the radius of the neighborhood.
    fn radius(&self) -> u8;

    /// Returns the number of positions in this neighborhood.
    fn n_neighs(&self) -> u16;

    /// Returns the positions in the neighborhood of `pos`.
    fn neighbors(&self, pos: Pos<isize>) -> impl Iterator<Item = Pos<isize>>;
}

/// Moore neighborhood with variable radius.
#[derive(Clone, Debug, PartialEq)]
pub struct MooreNeighborhood {
    radius: u8
}

impl MooreNeighborhood {
    /// Makes a new neighborhood with an associated `radius`.
    pub fn new(radius: u8) -> Self {
        Self { radius }
    }
}

impl Neighborhood for MooreNeighborhood {
    fn radius(&self) -> u8 {
        self.radius
    }

    #[inline]
    fn n_neighs(&self) -> u16 {
        4 * self.radius as u16 * (self.radius as u16 + 1)
    }

    #[inline]
    fn neighbors(&self, pos: Pos<isize>) -> impl Iterator<Item=Pos<isize>> {
        fetch_neighs(pos, &MOORE_NEIGHS, self.n_neighs())
    }
}

/// VonNeumann neighborhood with variable radius.
#[derive(Clone, Debug, PartialEq)]
pub struct VonNeumannNeighborhood {
    radius: u8,
}

impl VonNeumannNeighborhood {
    /// Makes a new neighborhood with an associated `radius`.
    pub fn new(radius: u8) -> Self {
        Self { radius }
    }
}

impl Neighborhood for VonNeumannNeighborhood {
    fn radius(&self) -> u8 {
        self.radius
    }

    #[inline]
    fn n_neighs(&self) -> u16 {
        2 * self.radius as u16 * (self.radius as u16 + 1)
    }

    #[inline]
    fn neighbors(&self, pos: Pos<isize>) -> impl Iterator<Item = Pos<isize>> {
        fetch_neighs(pos, &VON_NEUMANN_NEIGHS, self.n_neighs())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::assert;
    use std::collections::HashSet;
    #[test]
    fn test_moore() {
        let first_8 = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)];
        assert_eq!(first_8, MOORE_NEIGHS[..8]);
    }

    #[test]
    fn test_neighborhood_symmetry_moore() {
        let nh = MooreNeighborhood::new(3);
        let pos = Pos::new(0, 0);
        let neighs: HashSet<_> = nh.neighbors(pos).collect();
        for p in &neighs {
            assert!(neighs.contains(&Pos::new(-p.x, -p.y)), "Asymmetric Moore offset: {p:?}");
        }
    }

    #[test]
    fn test_neighborhood_symmetry_von_neumann() {
        let nh = VonNeumannNeighborhood::new(5);
        let pos = Pos::new(0, 0);
        let neighs: HashSet<_> = nh.neighbors(pos).collect();
        for p in &neighs {
            assert!(neighs.contains(&Pos::new(-p.x, -p.y)), "Asymmetric Von Neumann offset: {p:?}");
        }
    }
    
    #[test]
    fn test_radius_zero_returns_empty() {
        let moore = MooreNeighborhood::new(0);
        let von = VonNeumannNeighborhood::new(0);
        let pos = Pos::new(123, 456);
        assert_eq!(moore.neighbors(pos).count(), 0);
        assert_eq!(von.neighbors(pos).count(), 0);
    }

    #[test]
    fn test_neighs_do_not_include_center() {
        let pos = Pos::new(100, 100);
        let moore = MooreNeighborhood::new(1);
        let von = VonNeumannNeighborhood::new(1);

        assert!(!moore.neighbors(pos).any(|p| p == pos), "Moore included center");
        assert!(!von.neighbors(pos).any(|p| p == pos), "Von Neumann included center");
    }
}