use crate::boundary::Boundary;
use crate::cell::{Cell, CellCenter};
use crate::cell_container::CellContainer;
use crate::constants::{LatticeBoundaryType, Spin};
use crate::edge::{Edge, EdgeBook};
use crate::environment::LatticeEntity::*;
use crate::lattice::CellLattice;
use crate::neighbourhood::{MooreNeighbourhood, Neighbourhood};
use crate::pos::{AngularProjection, Pos2D, Rect};
use std::borrow::Borrow;
use crate::environment::DivisionError::{NewCellTooBig, NewCellTooSmall};

pub struct Environment {
    pub cell_lattice: CellLattice<LatticeBoundaryType>,
    pub cells: CellContainer,
    pub edge_book: EdgeBook,
    pub neighbourhood: MooreNeighbourhood
}

impl Environment {
    pub fn new(
        width: usize, 
        height: usize, 
        neigh_r: u8
    ) -> Self {
        let rect = Rect::new(
            (0, 0).into(),
            (width as isize, height as isize).into()
        );
        
         Self {
             cell_lattice: CellLattice::new(LatticeBoundaryType::new(rect)),
             cells: CellContainer::new(),
             edge_book: EdgeBook::new(),
             neighbourhood: MooreNeighbourhood::new(neigh_r)
        }
    }

    pub fn width(&self) -> usize {
        self.cell_lattice.width()
    }

    pub fn height(&self) -> usize {
        self.cell_lattice.height()
    }

    pub fn spawn_rect_cell(&mut self, rect: Rect<usize>, cell_target_area: u32) -> Option<&Cell> {
        let spin = self.cells.n_cells() as Spin + LatticeEntity::first_cell_spin();
        let center = self.cell_lattice.bound.valid_pos(Pos2D::new(
            rect.min.x as isize,
            rect.min.y as isize
        ));
        let mut cell = Cell::new(
            spin,
            0,
            cell_target_area,
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
        self.cells.push(cell);
        Some(self.cells.get_entity(spin).unwrap_cell())
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
    
    pub fn reproduce(&mut self, cell_target_area: u32, cell_div_area: u32) {
        let mut divide = vec![];
        for cell in &self.cells {
            if cell.area >= cell_div_area {
                divide.push(cell.spin);
            }
        }
        for spin in divide {
            if let Err(e) = self.divide_cell(spin, cell_target_area) {
                log::warn!("Failed to divide cell with spin {} with error `{:?}`", spin, e);
            }
        }
    }
    
    // We take spin here because this operation is not safe with &Cell (pushing to vec can cause reallocation)
    pub fn divide_cell(&mut self, spin: Spin, cell_target_area: u32) -> Result<&Cell, DivisionError> {
        let new_spin = self.cells.next_spin();
        let mom_cell = self
            .cells
            .get_entity_mut(spin)
            .expect_cell(&format!("passed non-cell with spin {} to `divide_cel()`", spin));
        let mut mom_clone = mom_cell.clone();
        
        let mut new_cell = Cell::new(
            new_spin,
            0,
            cell_target_area,
            CellCenter::origin()
        );
        let new_positions: Vec<_> = self
            .cell_lattice
            // TODO!: parameterise search radius
            .box_cell_positions(&mom_cell, 2.5)
            .into_iter()
            .filter(|pos| { 
                let proj = AngularProjection::from_pos(
                    Pos2D::new(pos.x as f32, pos.y as f32),
                    self.cell_lattice.width(),
                    self.cell_lattice.height()
                );
                // TODO!: use principal component to determine division axis
                //  current algorithm hands out all x positions to the right of the cell centre to the new cell
                proj.delta_angles(&mom_cell.center.projection).0 > 0.
            })
            .collect();
        for pos in new_positions {
            if mom_cell.area == 1 {
                return Err(NewCellTooBig);
            }
            self.cell_lattice[pos] = new_spin;
            new_cell.shift_position::<LatticeBoundaryType>(
                pos,
                self.cell_lattice.width(),
                self.cell_lattice.height(),
                true
            );
            mom_clone.shift_position::<LatticeBoundaryType>(
                pos,
                self.cell_lattice.width(),
                self.cell_lattice.height(),
                false
            );
        }
        if new_cell.area > 0 {
            mom_clone.target_area = cell_target_area;
            self.cells.replace(mom_clone);
            Ok(self.cells.push(new_cell))
        } else { 
            Err(NewCellTooSmall)
        }
    }

    /// Got tired of refactoring test and benchmark code
    pub fn empty_test(width:usize, height: usize) -> Self {
        Environment::new(
            width,
            height,
            1
        )
    }
}

#[derive(Debug)]
pub enum DivisionError {
    NewCellTooSmall,
    NewCellTooBig
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

    pub fn expect_cell(self, message: &str) -> C {
        match self {
            SomeCell(cell) => cell,
            _ => panic!("{}", message)
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
            ),
            0
        );
        assert_eq!(env.edge_book.len(), 8 * 4 * 3 + 4 * 5);
        let entity1 = env.cells.get_entity(LatticeEntity::first_cell_spin());
        assert!(matches!(entity1, SomeCell(_)));

        env.spawn_rect_cell(
            Rect::new(
                Pos2D::new(15, 15),
                Pos2D::new(25, 25)
            ),
            0
        );

        let entity2 = env.cells.get_entity(LatticeEntity::first_cell_spin() + 1);
        assert!(matches!(entity2, SomeCell(_)));

        let cell2 = entity2.unwrap_cell();
        assert_eq!(cell2.area, 75);
        assert_eq!(
            env.cell_lattice.contiguous_cell_positions(cell2, &env.neighbourhood).len(),
            75
        );
        assert_eq!(
            env.cell_lattice.box_cell_positions(cell2, 2.).len(),
            75
        );
    }

    #[test]
    fn test_lattice_entity_spin() {
        assert!(LatticeEntity::first_cell_spin() > Medium::<&Cell>.spin());
        assert!(LatticeEntity::first_cell_spin() > Solid::<&Cell>.spin());
    }
}
