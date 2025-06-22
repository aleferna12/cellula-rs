use crate::pos::Pos2D;

const MAX_NEIGH_R: u8 = 16;
const NEIGHBOURHOOD_SIZE: usize = 4 * MAX_NEIGH_R as usize * (MAX_NEIGH_R as usize + 1);
pub(crate) const MOORE_NEIGHS: [(i16, i16); NEIGHBOURHOOD_SIZE] = {
    let mut ret = [(0i16, 0i16); NEIGHBOURHOOD_SIZE];
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

pub trait Neighbourhood {
    fn radius(&self) -> u8;
    
    fn n_neighs(&self) -> u16;
    
    fn neighbours(&self, pos: Pos2D<isize>) -> impl Iterator<Item = Pos2D<isize>>;
}

pub struct MooreNeighbourhood {
    radius: u8
}
impl MooreNeighbourhood {
    pub fn new(radius: u8) -> Self {
        Self { radius }
    }
}

impl Neighbourhood for MooreNeighbourhood {
    fn radius(&self) -> u8 {
        self.radius
    }

    fn n_neighs(&self) -> u16 {
        4 * self.radius as u16 * (self.radius as u16 + 1)
    }

    fn neighbours(&self, pos: Pos2D<isize>) -> impl Iterator<Item=Pos2D<isize>> {
        let vec_size = 4 * self.radius as u16 * (self.radius as u16 + 1);
        MOORE_NEIGHS[..vec_size as usize]
            .iter()
            .map(move |(i, j)| {
                Pos2D::new(
                    pos.x + *i as isize,
                    pos.y + *j as isize,
                )
            })
    }
}

#[cfg(test)]
mod tests {
    use crate::edge::Edge;
    use super::*;
    use crate::pos::Pos2D;

    #[test]
    fn test_neighbours_are_edges() {
        let p1 = Pos2D::from((100, 100));
        for r in 1..9 {
            let neigh = MooreNeighbourhood::new(r);
            for p2 in neigh.neighbours(p1) {
                assert!(Edge::new_if_neighbour(p1.into(), p2.into(), r).is_ok());
            }
        }
    }

    #[test]
    fn test_moore() {
        let first_8 = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)];
        assert_eq!(first_8, MOORE_NEIGHS[..8]);
    }
}