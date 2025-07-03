use crate::boundary::Boundary;
use crate::cell::{Cell, CellCenter};
use crate::edge::{Edge, EdgeBook};
use crate::environment::LatticeEntity::*;
use crate::lattice::Lattice;
use crate::model::{LatticeBoundaryType, Spin};
use crate::neighbourhood::{MooreNeighbourhood, Neighbourhood};
use crate::pos::{Pos2D, Rect};
use std::borrow::Borrow;
use std::collections::hash_map::Entry::Vacant;
use std::collections::hash_map::IntoKeys;
use std::collections::{HashMap, VecDeque};

pub struct Environment {
    pub cell_lattice: Lattice<Spin, LatticeBoundaryType>,
    pub(crate) cell_vec: Vec<Cell>,
    pub edge_book: EdgeBook,
    pub neighbourhood: MooreNeighbourhood,
    pub cell_target_area: u32,
    pub cell_growth_period: u32,
    pub cell_div_area: u32
}

impl Environment {
    pub fn new(
        width: usize, 
        height: usize, 
        neigh_r: u8,
        cell_target_area: u32,
        cell_div_area: u32,
        cell_growth_period: u32,
    ) -> Self {
        let rect = Rect::new(
            (0, 0).into(),
            (width as isize, height as isize).into()
        );
        
         Self {
            cell_lattice: Lattice::new(LatticeBoundaryType::new(rect)),
            cell_vec: vec![],
            edge_book: EdgeBook::new(),
            neighbourhood: MooreNeighbourhood::new(neigh_r),
             cell_target_area,
             cell_div_area,
             cell_growth_period
        }
    }

    pub fn width(&self) -> usize {
        self.cell_lattice.width()
    }

    pub fn height(&self) -> usize {
        self.cell_lattice.height()
    }

    pub fn get_entity(&self, spin: Spin) -> LatticeEntity<&Cell> {
        if spin == Medium::<&Cell>.spin() {
            return Medium;
        }
        if spin == Solid::<&Cell>.spin() {
            return Solid;
        }
        SomeCell(&self.cell_vec[spin as usize - LatticeEntity::first_cell_spin() as usize])
    }

    pub fn get_entity_mut(&mut self, spin: Spin) -> LatticeEntity<&mut Cell> {
        if spin == Medium::<&Cell>.spin() {
            return Medium;
        }
        if spin == Solid::<&Cell>.spin() {
            return Solid;
        }
        SomeCell(&mut self.cell_vec[spin as usize - LatticeEntity::first_cell_spin() as usize])
    }
    
    pub fn n_cells(&self) -> usize {
        self.cell_vec.len()
    }

    pub fn spawn_rect_cell(&mut self, rect: Rect<usize>) -> Option<&Cell> {
        let spin = self.n_cells() as Spin + LatticeEntity::first_cell_spin();
        let center = self.cell_lattice.bound.valid_pos(Pos2D::new(
            rect.min.x as isize,
            rect.min.y as isize
        ));
        let mut cell = Cell::new(
            spin,
            0,
            self.cell_target_area,
            CellCenter::new(Pos2D::new(center?.x as f32, center?.y as f32), self.width(), self.height())
        );
        
        for pos in rect.iter_positions() {
            let trans_pos = self.cell_lattice.bound.valid_pos(pos.into());
            if trans_pos.is_none() {
                continue;
            }
            let valid_pos: Pos2D<usize> = trans_pos.unwrap().into();
            if self.cell_lattice[valid_pos] != Medium::<&Cell>.spin() {
                continue
            }
            self.cell_lattice[valid_pos] = spin;
            self.update_edges(valid_pos);
            cell.shift_position::<LatticeBoundaryType>(pos, self.width(), self.height(), true);
        }
        if cell.area == 0 { 
            return None;
        }
        self.cell_vec.push(cell);
        Some(self.get_entity(spin).unwrap_cell())
    }
    
    pub fn spawn_solid(&mut self, positions: impl Iterator<Item = Pos2D<usize>>) -> usize {
        let mut area = 0;
        for pos in positions {
            if self.cell_lattice[pos] != Medium::<&Cell>.spin() {
                continue
            }
            self.cell_lattice[pos] = Solid::<&Cell>.spin();
            area += 1;
        }
        area
    }
    
    pub fn make_border(&mut self) {
        let mut border_positions = Vec::<Pos2D<usize>>::new();
        for x in 0..self.width() {
            border_positions.push((x, 0).into());
        }
        for y in 1..self.height() {
            border_positions.push((self.width() - 1, y).into());
        }
        if self.width() > 1 {
            for y in (1..self.height() - 1).rev() {
                border_positions.push((0, y).into());
            }
        }
        if self.height() > 1 {
            for x in (0..self.width() - 1).rev() {
                border_positions.push((x, self.height() - 1).into());
            }
        }
        
        self.spawn_solid(border_positions.into_iter());
    }
    
    pub fn update_edges(&mut self, pos: Pos2D<usize>) -> (u16, u16) {
        let mut removed = 0;
        let mut added = 0;
        let spin = self.cell_lattice[pos];
        let valid_neighs = self.cell_lattice
            .bound
            .valid_positions(self.neighbourhood.neighbours(pos.into()))
            .map(|neigh| neigh.into());
        for neigh in valid_neighs {
            let edge = Edge::new(pos, neigh);
            let spin_neigh = self.cell_lattice[neigh];
            if spin == spin_neigh {
                self.edge_book.remove(&edge);
                // Also representing the inverse edge
                removed += 2;
                continue;
            }
            if spin < LatticeEntity::first_cell_spin() && spin_neigh < LatticeEntity::first_cell_spin() {
                continue;
            }
            if self.edge_book.insert(edge) {
                // Also representing the inverse edge
                added += 2;
            }
        }
        (removed, added)
    }
    
    pub fn update_cells(&mut self) {
        for cell in &mut self.cell_vec {
            if cell.area < self.cell_div_area {
                if cell.growth_timer >=  self.cell_growth_period {
                    cell.target_area += 1;
                    cell.growth_timer = 0;
                } else { 
                    cell.growth_timer += 1;
                }
            } else { // Divide
                cell.target_area = self.cell_target_area;
            }
        }
    }
    
    // TODO!: make a version where we iterate all positions in a large box and benchmark against this
    pub fn contiguous_cell_positions(&self, cell: &Cell) -> IntoKeys<Pos2D<usize>, ()> {
        let mut found = HashMap::<Pos2D<usize>, ()>::default();
        let mut deque = VecDeque::from([Pos2D::new(
            cell.center.pos.x as isize,
            cell.center.pos.y as isize
        )]);

        while !deque.is_empty() {
            let pos = deque.pop_front().unwrap();
            let lat_pos = Pos2D::<usize>::from(pos);
            if cell.spin != self.cell_lattice[lat_pos] {
                continue;
            }
            
            if let Vacant(entry) = found.entry(lat_pos) {
                let neighs = self
                    .cell_lattice
                    .bound
                    .valid_positions(self.neighbourhood.neighbours(pos));
                for neigh in neighs {
                    deque.push_back(neigh);
                }
                entry.insert(());
            }
        }
        found.into_keys()
    }
    
    /// Got tired of refactoring test and benchmark code
    pub fn empty_test(width:usize, height: usize) -> Self {
        Environment::new(
            width,
            height,
            1,
            64,
            1,
            64
        )
    }
}

/// This enum represents anything that can be on the cell lattice.
#[derive(Debug, Copy, Clone)]
pub enum LatticeEntity<C> {
    Solid,
    Medium,
    SomeCell(C),
}

impl<C> LatticeEntity<C> {
    pub fn map<D, F: FnOnce(C) -> D>(self, f: F) -> LatticeEntity<D> {
        match self {
            SomeCell(c) => SomeCell(f(c)),
            Medium => Medium,
            Solid => Solid,
        }
    }
}

impl<C: Borrow<Cell>> LatticeEntity<C> {
    /// Maps the `LatticeEntity` to a unique spin value.
    pub fn spin(&self) -> Spin {
        match self {
            SomeCell(cell) => cell.borrow().spin,
            Medium => 0,
            Solid => 1
        }
    }
}

impl<C: std::fmt::Debug> LatticeEntity<C> {
    pub fn unwrap_cell(self) -> C {
        match self {
            SomeCell(cell) => cell,
            _ => panic!("called `LatticeEntity::unwrap_cell()` on a `{:?}` value", self)
        }
    }
}

impl LatticeEntity<()> {
    /// Returns the first spin that corresponds to a cell.
    /// 
    /// This is required to be larger than `Medium::spin()` and `Solid::spin()`.
    pub fn first_cell_spin() -> Spin {
        2
    }
}

#[cfg(test)]
pub mod tests {
    use super::*;
    
    #[test]
    fn test_spawn_rect_cell() {
        let mut env = Environment::empty_test(100, 100);
        env.spawn_rect_cell(
            Rect::new(
                Pos2D::new(10, 10),
                Pos2D::new(20, 20)
            )
        );
        assert_eq!(env.edge_book.len(), 8 * 4 * 3 + 4 * 5);
        let entity1 = env.get_entity(LatticeEntity::first_cell_spin());
        assert!(matches!(entity1, SomeCell(_)));
        
        env.spawn_rect_cell(
            Rect::new(
                Pos2D::new(15, 15),
                Pos2D::new(25, 25)
            )
        );
        
        let entity2 = env.get_entity(LatticeEntity::first_cell_spin() + 1);
        assert!(matches!(entity2, SomeCell(_)));

        let cell2 = entity2.unwrap_cell();
        assert_eq!(cell2.area, 75);
        assert_eq!(env.contiguous_cell_positions(cell2).count(), 75);
    }

    #[test]
    fn test_lattice_entity_spin() {
        assert!(LatticeEntity::first_cell_spin() > Medium::<&Cell>.spin());
        assert!(LatticeEntity::first_cell_spin() > Solid::<&Cell>.spin());
    }
}
