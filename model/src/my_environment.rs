//! Contains logic associated with [`MyEnvironment`].

use crate::constants::{BoundaryType, NeighborhoodType, EPSILON};
use crate::my_cell::MyCell;
use cellulars::prelude::*;
use rand::RngExt;

/// An environment that contains a chemical gradient and limits cell growth to [`MyEnvironment::max_cells`].
#[derive(Debug, Clone)]
pub struct MyEnvironment {
    /// Inner [`Environment`].
    pub env: Environment<MyCell, NeighborhoodType, BoundaryType>,
    /// Lattice containing the chemical gradient.
    pub chem_lattice: Lattice<FloatType>,
    /// Scaler used to determine the radius of search for cell positions starting from its center.
    pub cell_search_scaler: FloatType,
    /// Maximum number of cells supported in the environment.
    pub max_cells: CellIndex,
    population_exploded: bool
}

impl MyEnvironment {
    /// Make a new [`MyEnvironment`] from an existing [`Environment`].
    pub fn new(
        env: Environment<MyCell, NeighborhoodType, BoundaryType>,
        max_cells: CellIndex,
        cell_search_scaler: FloatType
    ) -> Self {
        let mut env_ = Self {
            chem_lattice: Lattice::from(env.cell_lattice.rect.clone()),
            env,
            cell_search_scaler,
            max_cells,
            population_exploded: false
        };
        env_.make_chem_gradient();
        env_
    }

    /// Creates a chemical gradient spanning from the top to the bottom of the environment.
    pub fn make_chem_gradient(&mut self) {
        for i in 0..self.env.width() {
            for j in 0..self.env.height() {
                self.chem_lattice[(i, j).into()] = j as FloatType;
            }
        }
    }

    /// Returns whether the environment supports additional cells based on [`MyEnvironment::max_cells`].
    pub fn can_add_cell(&mut self) -> bool {
        if self.env.cells.n_non_empty() < self.max_cells {
            return true;
        }
        if !self.population_exploded {
            log::warn!(
                        "Population exceeded maximum threshold `max-cells={}` during cell division",
                        {self.max_cells}
                    );
            log::warn!("This warning will be suppressed from now on");
            self.population_exploded = true;
        }
        false
    }

    /// Spawns a square cell centered at a random position with area = `cell_area`.
    ///
    /// Uses [`MyEnvironment::spawn_cell_checked()`] to restrict spawns to the medium.
    pub fn spawn_cell_random(
        &mut self,
        empty_cell: EmptyCell<MyCell>,
        cell_area: u32,
        rng: &mut impl RngExt,
    ) -> &RelCell<MyCell> {
        let pos_isize = self
            .env
            .cell_lattice
            .random_pos(rng)
            .cast_as::<isize>();
        let cell_side = ((cell_area as FloatType).sqrt() / 2.).floor() as isize;
        let rect = Rect::new(
            Pos::new(pos_isize.x - cell_side, pos_isize.y - cell_side),
            Pos::new(pos_isize.x + cell_side, pos_isize.y + cell_side)
        );
        self.spawn_cell_checked(
            empty_cell,
            rect.iter_positions()
        )
    }

    /// Divides a cell along its minor axis.
    pub fn divide_cell(&mut self, mom_index: CellIndex) -> &RelCell<MyCell> {
        let rel_mom = &self.env.cells[mom_index];
        let div_axis = self.find_division_axis(rel_mom);
        let normal = (-div_axis.1, div_axis.0);
        let new_positions: Box<_> = self
            .env
            .search_cell_box(rel_mom, self.cell_search_scaler)
            .into_iter()
            .filter(|pos| {
                let (dx, dy) = self.env.bounds.boundary.displacement(
                    pos.cast_as(),
                    rel_mom.cell.center()
                );
                dx * normal.0 + dy * normal.1 < 0.
            })
            .collect();
        
        let newborn_ta = rel_mom.cell.newborn_target_area;
        let newborn = rel_mom.cell.birth();
        let new_index = self.env.cells.add(newborn).index;
        for pos in new_positions {
            self.transfer_position(
                pos,
                Spin::Some(new_index),
            );
        }
        self.env.cells[mom_index].cell.cell.target_area = newborn_ta;
        &self.env.cells[new_index]
    }

    /// Checks which cells should divide and executes cell divisions.
    pub fn reproduce(&mut self) {
        let mut divide = vec![];
        for rel_cell in self.env.cells.iter() {
            if !rel_cell.cell.is_alive() {
                continue;
            }
            // Cells don't need to express the dividing type to divide, they just need to be big enough
            if rel_cell.cell.area() >= rel_cell.cell.divide_area {
                divide.push(rel_cell.index);
            }
        }
        for cell_index in divide {
            if !self.can_add_cell() {
                return;
            }
            self.divide_cell(cell_index);
        }
    }

    /// Finds the minor axis along which to split the cell.
    pub fn find_division_axis(&self, rel_cell: &RelCell<MyCell>) -> (FloatType, FloatType) {
        // Compute covariance elements relative to centroid
        let mut sum_xx = 0.0;
        let mut sum_yy = 0.0;
        let mut sum_xy = 0.0;

        for p in &self.env.search_cell_box(rel_cell, self.cell_search_scaler) {
            let (dx, dy) = self.env.bounds.boundary.displacement(
                p.cast_as(),
                rel_cell.cell.center()
            );
            sum_xx += dx * dx;
            sum_yy += dy * dy;
            sum_xy += dx * dy;
        }

        let n = rel_cell.cell.area() as FloatType;
        let cov_xx = sum_xx / n;
        let cov_yy = sum_yy / n;
        let cov_xy = sum_xy / n;

        // Eigenvalues of covariance matrix:
        // λ = (trace ± sqrt((cov_xx - cov_yy)^2 + 4*cov_xy^2)) / 2
        let trace = cov_xx + cov_yy;
        let determinant = cov_xx * cov_yy - cov_xy * cov_xy;
        let discriminant = (trace * trace - 4.0 * determinant).sqrt();
        let lambda2 = (trace - discriminant) / 2.0; // smaller eigenvalue

        // Eigenvector for the minor axis (lambda2)
        let (vec_x, vec_y) = if cov_xy.abs() > EPSILON {
            // Solve (C - λI)v = 0
            (lambda2 - cov_yy, cov_xy)
        } else {
            // Axis-aligned case
            if cov_xx < cov_yy {
                (1.0, 0.0) // x-axis is minor
            } else {
                (0.0, 1.0) // y-axis is minor
            }
        };

        // Normalize vector
        let norm = (vec_x * vec_x + vec_y * vec_y).sqrt();
        let vec_x = vec_x / norm;
        let vec_y = vec_y / norm;

        (vec_x, vec_y)
    }

    /// Removes all cells from the environment and restore it to a clean state.
    pub fn wipe_out(&mut self) {
        self.env.wipe_out();
    }

    /// Creates a border of [`Spin::Solid`] around the environment.
    pub fn make_border(
        &mut self,
        bottom: bool,
        top: bool,
        left: bool,
        right: bool,
    ) {
        let mut border_positions = Vec::<Pos<usize>>::new();
        if bottom {
            for x in 0..self.env.width() {
                border_positions.push((x, 0).into());
            }
        }
        if top {
            for x in (0..self.env.width() - 1).rev() {
                border_positions.push((x, self.env.height() - 1).into());
            }
        }
        if left {
            for y in (1..self.env.height() - 1).rev() {
                border_positions.push((0, y).into());
            }
        }
        if right {
            for y in 1..self.env.height() {
                border_positions.push((self.env.width() - 1, y).into());
            }
        }

        self.spawn_solid(border_positions.into_iter());
    }

    /// Spawns an `empty_cell` on valid `positions` that belong to the medium,
    /// while ignoring positions owned by solids or other cells.
    pub fn spawn_cell_checked(
        &mut self,
        empty_cell: EmptyCell<MyCell>,
        positions: impl IntoIterator<Item = Pos<isize>>
    ) -> &RelCell<MyCell> {
        let med_positions = positions.into_iter().filter_map(|pos| {
            let valid_pos = self.env.bounds.lattice_boundary.valid_pos(pos)?;
            let lat_pos = valid_pos.cast_as();
            if !matches!(self.env.cell_lattice[lat_pos], Spin::Medium) {
                return None;
            }
            Some(lat_pos)
        }).collect::<Box<[_]>>();
        self.spawn_cell(empty_cell, med_positions)
    }
}

impl TransferPosition for MyEnvironment {
    fn transfer_position(
        &mut self,
        pos: Pos<usize>,
        to: Spin
    ) -> EdgesUpdate {
        let chem_at_pos = self.chem_lattice[pos];
        if let Spin::Some(index) = to {
            self.env.cells[index].cell.shift_chem(pos, chem_at_pos, true, &self.env.bounds.boundary);
        }
        if let Spin::Some(index) = self.env.cell_lattice[pos] {
            let from_rel_cell = &mut self.env.cells[index];
            from_rel_cell.cell.shift_chem(pos, chem_at_pos, false, &self.env.bounds.boundary);
            // If the copy kills the cell
            if from_rel_cell.cell.area() == 0 {
                from_rel_cell.cell.apoptosis();
            }
        }
        self.env.transfer_position(pos, to)
    }
}

impl AsEnv for MyEnvironment {
    type Cell = MyCell;

    fn env(&self) -> &Environment<Self::Cell, impl Neighborhood, impl ToLatticeBoundary> {
        &self.env
    }

    fn env_mut(&mut self) -> &mut Environment<Self::Cell, impl Neighborhood, impl ToLatticeBoundary> {
        &mut self.env
    }
}

impl Spawn for MyEnvironment {}