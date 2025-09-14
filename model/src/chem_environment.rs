use crate::cell::Cell;
use crate::constants::{BoundaryType, EPSILON};
use cellulars_lib::basic_cell::{Alive, Cellular, RelCell};
use cellulars_lib::cell_container::CellContainer;
use cellulars_lib::constants::Spin;
use cellulars_lib::environment::{EdgesUpdate, Environment, Habitable};
use cellulars_lib::lattice::Lattice;
use cellulars_lib::lattice_entity::LatticeEntity::SomeCell;
use cellulars_lib::positional::neighbourhood::MooreNeighbourhood;
use cellulars_lib::positional::pos::Pos;
use std::ops::{Deref, DerefMut};
use rand::Rng;
use cellulars_lib::positional::boundaries::Boundary;
use cellulars_lib::positional::rect::Rect;

#[derive(Clone)]
pub struct ChemEnvironment {
    env: Environment<Cell, MooreNeighbourhood, BoundaryType>,
    pub chem_lattice: Lattice<u32>,
    pub max_cells: Spin,
    population_exploded: bool
}

impl ChemEnvironment {
    pub fn new(env: Environment<Cell, MooreNeighbourhood, BoundaryType>, max_cells: Spin) -> Self {
        let mut env_ = Self {
            chem_lattice: env.cell_lattice.clone(),
            env,
            max_cells,
            population_exploded: false
        };
        env_.make_chem_gradient();
        env_
    }

    pub fn make_chem_gradient(&mut self) {
        for row in 0..self.height() {
            for col in 0..self.width() {
                self.chem_lattice[(col, row).into()] = row.try_into().expect("lattice is too big");
            }
        }
    }

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
        let positions = rect
            .iter_positions()
            .filter_map(|pos| self.env.bounds.lattice_boundary.valid_pos(pos))
            .map(|pos| pos.to_usize())
            .collect::<Vec<_>>();
        self.env.spawn_cell(
            empty_cell,
            positions
        )
    }

    pub fn divide_cell(&mut self, mom_spin: Spin, search_scaler: f32) -> &RelCell<Cell> {
        let mom = self
            .env
            .cells
            .get_entity(mom_spin)
            .expect_cell("retrieved non-cell during cell division");
        let new_positions: Vec<_> = self
            .env
            .search_cell_box(mom, search_scaler)
            .into_iter()
            .filter(|pos| {
                // TODO!: use principal component to determine division axis
                //  current algorithm hands out all x positions to the right of the cell centre to the new cell
                self.env.bounds.boundary.displacement(Pos::new(pos.x as f32, pos.y as f32), mom.center()).0 > 0.
            })
            .collect();

        let newborn = mom.birth();
        let new_spin = self.env.cells.add(newborn, Some(mom_spin)).spin;
        for pos in new_positions {
            self.env.grant_position(
                pos,
                new_spin,
            );
        }
        self.env.cells.get_entity(new_spin).expect_cell("retrieved non-cell during cell division")
    }

    // Should this also replace some of the cell's positions with Medium?
    pub fn kill_cell(&mut self, cell: &mut RelCell<Cell>) {
        cell.apoptosis();
    }

    // With some unsafe code we can return Vec<&RelCell> from this function, but it would
    // require that self.divide_cell never invalidates any references to self.cells
    // we need thorough testing of self.divide_cells to make this change, and the performance
    // gain is minimal (although the ergonomic gains are significant)
    pub fn reproduce(&mut self, search_scaler: f32) -> Vec<Spin> {
        let mut divide = vec![];
        for cell in self.env.cells().iter() {
            if !cell.is_alive() {
                continue;
            }
            // Currently cells don't need to express the dividing type to divide, they just need to be big enough
            if cell.area() >= cell.divide_area() {
                divide.push(cell.spin);
            }
        }
        divide.into_iter().filter_map(|spin| {
            if !self.can_add_cell() {
                return None;
            }
            let mom = self.env
                .cells
                .get_entity(spin)
                .expect_cell("retrieved non-cell during reproduction");
            Some(self.divide_cell(mom.spin, search_scaler).spin)
        }).collect()
    }

    /// Find the minor axis along which to split the cell.
    pub fn find_division_axis(&self, cell: &RelCell<Cell>, search_scaler: f32) -> SplitLine {
        let center_x = cell.center().x;
        let center_y = cell.center().y;

        // Compute covariance elements relative to centroid
        let mut sum_xx = 0.0;
        let mut sum_yy = 0.0;
        let mut sum_xy = 0.0;

        for p in &self.search_cell_box(cell, search_scaler) {
            let dx = p.x as f32 - center_x;
            let dy = p.y as f32 - center_y;
            sum_xx += dx * dx;
            sum_yy += dy * dy;
            sum_xy += dx * dy;
        }

        let n = cell.area as f32;
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
        let intercept = center_y - slope * center_x;

        SplitLine { slope, intercept }
    }
}

impl Deref for ChemEnvironment {
    type Target = Environment<Cell, MooreNeighbourhood, BoundaryType>;
    fn deref(&self) -> &Self::Target {
        &self.env
    }
}

impl DerefMut for ChemEnvironment {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.env
    }
}

impl Habitable for ChemEnvironment {
    type Cell = Cell;

    fn cells(&self) -> &CellContainer<Self::Cell> {
        self.env.cells()
    }

    fn cells_mut(&mut self) -> &mut CellContainer<Self::Cell> {
        self.env.cells_mut()
    }

    fn grant_position(
        &mut self,
        pos: Pos<usize>,
        to: Spin
    ) -> EdgesUpdate {
        // TODO! chem should always be u32
        let chem_at_pos = self.chem_lattice[pos] as f32;
        if let SomeCell(to_cell) = self.env.cells.get_entity_mut(to) {
            to_cell.shift_position(pos, true, &self.env.bounds.boundary);
            to_cell.shift_chem(pos, chem_at_pos, true, &self.env.bounds.boundary);
        }
        let from = self.cell_lattice[pos];
        if let SomeCell(from_cell) = self.env.cells.get_entity_mut(from) {
            from_cell.shift_position(pos, false, &self.env.bounds.boundary);
            from_cell.shift_chem(pos, chem_at_pos, false, &self.env.bounds.boundary);
        }
        // Executes the copy
        self.cell_lattice[pos] = to;
        self.update_edges(pos)
    }
}

/// Represents a split line through the centroid:
///   y = slope * x + intercept
#[derive(Debug)]
pub struct SplitLine {
    pub slope: f32,
    pub intercept: f32,
}