//! Contains logic associated with [`MyEnvironment`].

use crate::my_cell::MyCell;
use crate::constants::{BoundaryType, NeighbourhoodType, EPSILON};
use cellulars_lib::base::environment::{Environment, EdgesUpdate};
use cellulars_lib::cell_container::RelCell;
use cellulars_lib::constants::{CellIndex, FloatType};
use cellulars_lib::lattice::Lattice;
use cellulars_lib::positional::boundaries::{Boundary, ToLatticeBoundary};
use cellulars_lib::positional::neighbourhood::Neighbourhood;
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::positional::rect::Rect;
use cellulars_lib::spin::Spin;
use cellulars_lib::traits::cellular::{Alive, Cellular, EmptyCell};
use cellulars_lib::traits::habitable::Habitable;
use rand::Rng;

/// An environment that contains a chemical gradient and limits cell growth to [`MyEnvironment::max_cells`].
#[derive(Clone)]
pub struct MyEnvironment {
    /// Inner [`Environment`].
    pub env: Environment<MyCell, NeighbourhoodType, BoundaryType>,
    /// Lattice containing the chemical gradient.
    pub chem_lattice: Lattice<u32>,
    /// Scaler used to determine the radius of search for cell positions starting from its center.
    pub cell_search_scaler: FloatType,
    /// Maximum number of cells supported in the environment.
    pub max_cells: CellIndex,
    population_exploded: bool
}

impl MyEnvironment {
    /// Make a new [`MyEnvironment`] from an existing [`Environment`].
    pub fn new(
        env: Environment<MyCell, NeighbourhoodType, BoundaryType>,
        max_cells: CellIndex,
        cell_search_scaler: FloatType
    ) -> Self {
        let mut env_ = Self {
            chem_lattice: Lattice::new(env.cell_lattice.rect.clone()),
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
                self.chem_lattice[(i, j).into()] = j.try_into().expect("lattice is too big");
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

    // TODO: make spawn as a circle with center at pos
    /// Spawns a square cell centered at a random position with area = `cell_area`.
    ///
    /// Uses [`MyEnvironment::spawn_cell_checked()`] to restrict spawns to the medium.
    pub fn spawn_cell_random(
        &mut self,
        empty_cell: EmptyCell<MyCell>,
        cell_area: u32,
        rng: &mut impl Rng,
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

    /// Forces a cell to execute cell division.
    pub fn divide_cell(&mut self, mom_index: CellIndex) -> &RelCell<MyCell> {
        let rel_mom = &self.env.cells[mom_index];
        // TODO!: This searches cell positions twice (once to find div axis).
        let div_axis = self.find_division_axis(rel_mom, self.cell_search_scaler);
        let new_positions: Box<_> = self
            .env
            .search_cell_box(rel_mom, self.cell_search_scaler)
            .into_iter()
            .filter(|pos| {
                let y = div_axis.slope * pos.x as FloatType + div_axis.intercept;
                (pos.y as FloatType) < y
            })
            .collect();
        
        let newborn_ta = rel_mom.cell.newborn_target_area;
        let newborn = rel_mom.cell.birth();
        let new_index = self.env.cells.add(newborn).index;
        for pos in new_positions {
            self.grant_position(
                pos,
                Spin::Some(new_index),
            );
        }
        self.env.cells[mom_index].cell.cell.target_area = newborn_ta;
        &self.env.cells[new_index]
    }

    // With some unsafe code we can return Vec<&RelCell> from this function, but it would
    // require that self.divide_cell never invalidates any references to self.cells
    // we need thorough testing of self.divide_cells to make this change, and the performance
    // gain is minimal (although the ergonomic gains are significant)
    /// Checks which cells should divide and executes cell divisions.
    pub fn reproduce(&mut self) {
        let mut divide = vec![];
        for rel_cell in self.env.cells.iter() {
            if !rel_cell.cell.is_alive() {
                continue;
            }
            // Currently cells don't need to express the dividing type to divide, they just need to be big enough
            if rel_cell.cell.area() >= rel_cell.cell.divide_area {
                divide.push(rel_cell.index);
            }
        }
        for cell_index in divide {
            if !self.can_add_cell() {
                return;
            }

            let mom = &self.env.cells[cell_index];
            self.divide_cell(mom.index);
        }
    }

    // TODO!: add plot to make sure this is right
    /// Finds the minor axis along which to split the cell.
    pub fn find_division_axis(&self, rel_cell: &RelCell<MyCell>, search_scaler: FloatType) -> SplitLine {
        // Compute covariance elements relative to centroid
        let mut sum_xx = 0.0;
        let mut sum_yy = 0.0;
        let mut sum_xy = 0.0;

        for p in &self.env.search_cell_box(rel_cell, search_scaler) {
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

        // Line equation through centroid with this direction
        let slope = if vec_x.abs() > EPSILON {
            vec_y / vec_x
        } else {
            FloatType::INFINITY // vertical line
        };
        let intercept = rel_cell.cell.center().y - slope * rel_cell.cell.center().x;

        SplitLine { slope, intercept }
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

impl Habitable for MyEnvironment {
    type Cell = MyCell;

    fn env(&self) -> &Environment<Self::Cell, impl Neighbourhood, impl ToLatticeBoundary> {
        &self.env
    }

    fn env_mut(&mut self) -> &mut Environment<Self::Cell, impl Neighbourhood, impl ToLatticeBoundary> {
        &mut self.env
    }

    fn grant_position(
        &mut self,
        pos: Pos<usize>,
        to: Spin
    ) -> EdgesUpdate {
        let chem_at_pos = self.chem_lattice[pos];
        if let Spin::Some(index) = to {
            let to_rel_cell = &mut self.env.cells[index];
            to_rel_cell.cell.shift_position(pos, true, &self.env.bounds.boundary);
            to_rel_cell.cell.shift_chem(pos, chem_at_pos, true, &self.env.bounds.boundary);
        }
        if let Spin::Some(index) = self.env.cell_lattice[pos] {
            let from_rel_cell = &mut self.env.cells[index];
            from_rel_cell.cell.shift_position(pos, false, &self.env.bounds.boundary);
            from_rel_cell.cell.shift_chem(pos, chem_at_pos, false, &self.env.bounds.boundary);
            // If the copy kills the cell
            if from_rel_cell.cell.area() == 0 {
                from_rel_cell.cell.apoptosis();
            }
        }
        // Executes the copy
        self.env.cell_lattice[pos] = to;
        self.env.update_edges(pos)
    }
}

/// Represents a split line through the centroid: y = `slope` * x + `intercept`.
#[derive(Debug)]
pub struct SplitLine {
    /// Slope of the linear equation.
    pub slope: FloatType,
    /// Intercept of the linear equation.
    pub intercept: FloatType,
}