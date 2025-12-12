//! Contains logic associated with [MyEnvironment].

use crate::cell::Cell;
use crate::constants::{BoundaryType, EPSILON};
use cellulars_lib::basic_cell::{Alive, Cellular, RelCell};
use cellulars_lib::constants::CellIndex;
use cellulars_lib::environment::{EdgesUpdate, Environment};
use cellulars_lib::habitable::Habitable;
use cellulars_lib::lattice::Lattice;
use cellulars_lib::positional::boundaries::{Boundary, ToLatticeBoundary};
use cellulars_lib::positional::neighbourhood::{MooreNeighbourhood, Neighbourhood};
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::positional::rect::Rect;
use cellulars_lib::spin::Spin;
use rand::Rng;

/// An environment that contains a chemical gradient and limits cell growth to [MyEnvironment::max_cells].
#[derive(Clone)]
pub struct MyEnvironment {
    env: Environment<Cell, MooreNeighbourhood, BoundaryType>,
    /// Lattice containing the chemical gradient.
    pub chem_lattice: Lattice<u32>,
    /// Scaler used to determine the radius of search for cell positions starting from its center.
    pub cell_search_scaler: f32,
    /// Maximum number of cells supported in the environment.
    pub max_cells: CellIndex,
    population_exploded: bool
}

impl MyEnvironment {
    /// Make a new [MyEnvironment] from an existing [Environment].
    pub fn new(
        env: Environment<Cell, MooreNeighbourhood, BoundaryType>, 
        max_cells: CellIndex,
        cell_search_scaler: f32
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

    /// Returns a reference to the inner [cellulars_lib::environment::Environment](Environment).
    pub fn env(&self) -> &Environment<Cell, MooreNeighbourhood, BoundaryType> {
        &self.env
    }

    /// Returns a mutable reference to the inner [cellulars_lib::environment::Environment](Environment).
    pub fn env_mut(&mut self) -> &mut Environment<Cell, MooreNeighbourhood, BoundaryType> {
        &mut self.env
    }

    /// Creates a chemical gradient spanning from the top to the bottom of the environment.
    pub fn make_chem_gradient(&mut self) {
        for i in 0..self.env.width() {
            for j in 0..self.env.height() {
                self.chem_lattice[(i, j).into()] = j.try_into().expect("lattice is too big");
            }
        }
    }

    /// Returns whether the environment supports additional cells based on [MyEnvironment::max_cells].
    pub fn can_add_cell(&mut self) -> bool {
        if self.env.cells.n_valid() < self.max_cells {
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
    pub fn spawn_cell_random(
        &mut self,
        empty_cell: Cell,
        cell_area: u32,
        rng: &mut impl Rng,
    ) -> &RelCell<Cell> {
        let pos_isize = self.env.cell_lattice.random_pos(rng).to_isize();
        let cell_side = ((cell_area as f32).sqrt() / 2.) as isize;
        let rect = Rect::new(
            Pos::new(pos_isize.x - cell_side, pos_isize.y - cell_side),
            Pos::new(pos_isize.x + cell_side, pos_isize.y + cell_side)
        );
        let positions: Box<_> = rect
            .iter_positions()
            .filter_map(|pos| self.env.bounds.lattice_boundary.valid_pos(pos))
            .map(|pos| pos.to_usize())
            .collect();
        self.spawn_cell(
            empty_cell,
            positions
        )
    }

    /// Forces a cell to execute cell division.
    pub fn divide_cell(&mut self, mom_index: CellIndex) -> &RelCell<Cell> {
        let mom = self
            .env
            .cells
            .get_cell(mom_index);
        // TODO!: This searches cell positions twice (once to find div axis).
        let div_axis = self.find_division_axis(mom, self.cell_search_scaler);
        let new_positions: Box<_> = self
            .env
            .search_cell_box(mom, self.cell_search_scaler)
            .into_iter()
            .filter(|pos| {
                let y = div_axis.slope * pos.x as f32 + div_axis.intercept;
                (pos.y as f32) < y
            })
            .collect();
        
        let newborn = mom.birth();
        let newborn_ta = mom.newborn_target_area;
        let new_index = self.env.cells.add(newborn).index;
        for pos in new_positions {
            self.grant_position(
                pos,
                Spin::Some(new_index),
            );
        }
        self.env.cells.get_cell_mut(mom_index).basic_cell_mut().target_area = newborn_ta;
        self.env.cells.get_cell(new_index)
    }

    // With some unsafe code we can return Vec<&RelCell> from this function, but it would
    // require that self.divide_cell never invalidates any references to self.cells
    // we need thorough testing of self.divide_cells to make this change, and the performance
    // gain is minimal (although the ergonomic gains are significant)
    /// Checks which cells should divide and executes cell divisions.
    pub fn reproduce(&mut self) {
        let mut divide = vec![];
        for cell in self.env.cells.iter() {
            if !cell.is_alive() {
                continue;
            }
            // Currently cells don't need to express the dividing type to divide, they just need to be big enough
            if cell.area() >= cell.divide_area {
                divide.push(cell.index);
            }
        }
        for cell_index in divide {
            if !self.can_add_cell() {
                return;
            }

            let mom = self
                .env
                .cells
                .get_cell(cell_index);
            self.divide_cell(mom.index);
        }
    }

    // TODO!: add plot to make sure this is right
    /// Finds the minor axis along which to split the cell.
    pub fn find_division_axis(&self, cell: &RelCell<Cell>, search_scaler: f32) -> SplitLine {
        // Compute covariance elements relative to centroid
        let mut sum_xx = 0.0;
        let mut sum_yy = 0.0;
        let mut sum_xy = 0.0;

        for p in &self.env.search_cell_box(cell, search_scaler) {
            let (dx, dy) = self.env.bounds.boundary.displacement(p.to_f32(), cell.center());
            sum_xx += dx * dx;
            sum_yy += dy * dy;
            sum_xy += dx * dy;
        }

        let n = cell.area() as f32;
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
            f32::INFINITY // vertical line
        };
        let intercept = cell.center().y - slope * cell.center().x;

        SplitLine { slope, intercept }
    }

    /// Removes all cells from the environment and restore it to a clean state.
    pub fn wipe_out(&mut self) {
        self.env.wipe_out();
    }

    /// Creates a border of [Spin::Solid] around the environment.
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
}

impl Habitable for MyEnvironment {
    type Cell = Cell;

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
            let to_cell = self.env.cells.get_cell_mut(index);
            to_cell.shift_position(pos, true, &self.env.bounds.boundary);
            to_cell.shift_chem(pos, chem_at_pos, true, &self.env.bounds.boundary);
        }
        if let Spin::Some(index) = self.env.cell_lattice[pos] {
            let from_cell = self.env.cells.get_cell_mut(index);
            from_cell.shift_position(pos, false, &self.env.bounds.boundary);
            from_cell.shift_chem(pos, chem_at_pos, false, &self.env.bounds.boundary);
            // If the copy kills the cell
            if from_cell.area() == 0 {
                from_cell.apoptosis();
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
    pub slope: f32,
    /// Intercept of the linear equation.
    pub intercept: f32,
}