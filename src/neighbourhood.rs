use crate::pos::Pos2D;

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
        crate::pos::MOORE_NEIGHS[..vec_size as usize]
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
}