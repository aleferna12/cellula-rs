use crate::cell::Cell;
use crate::constants::{BoundaryType, EPSILON};
use crate::evolution::genome::Genome;
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
use cellulars_lib::symmetric_table::SymmetricTable;
use rand::Rng;
use std::collections::HashSet;
use std::ops::{Deref, DerefMut};

#[derive(Clone)]
pub struct ChemEnvironment {
    env: Environment<Cell, MooreNeighbourhood, BoundaryType>,
    pub chem_lattice: Lattice<u32>,
    pub act_lattice: Lattice<u32>,
    pub clones_table: SymmetricTable<bool>,
    pub cell_search_scaler: f32,
    pub max_cells: CellIndex,
    act_max: u32,
    population_exploded: bool
}

impl ChemEnvironment {
    pub fn new(env: Environment<Cell, MooreNeighbourhood, BoundaryType>, max_cells: CellIndex, act_max: u32, cell_search_scaler: f32) -> Self {
        let lat = Lattice::new(env.cell_lattice.rect.clone());
        let mut env_ = Self {
            chem_lattice: lat.clone(),
            act_lattice: lat,
            clones_table: SymmetricTable::new(max_cells as usize),
            population_exploded: false,
            env,
            cell_search_scaler,
            max_cells,
            act_max
        };
        env_.make_chem_gradient();
        env_
    }

    pub fn make_chem_gradient(&mut self) {
        for i in 0..self.width() {
            for j in 0..self.height() {
                self.chem_lattice[(i, j).into()] = j.try_into().expect("lattice is too big");
            }
        }
    }

    pub fn can_add_cell(&mut self) -> bool {
        if self.cells.n_valid() < self.max_cells {
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
        let pos_isize = self.cell_lattice.random_pos(rng).to_isize();
        let cell_side = ((cell_area as f32).sqrt() / 2.) as isize;
        let rect = Rect::new(
            Pos::new(pos_isize.x - cell_side, pos_isize.y - cell_side),
            Pos::new(pos_isize.x + cell_side, pos_isize.y + cell_side)
        );
        let positions = rect
            .iter_positions()
            .filter_map(|pos| self.bounds.lattice_boundary.valid_pos(pos))
            .map(|pos| pos.to_usize())
            .collect::<Vec<_>>();
        self.spawn_cell(
            empty_cell,
            positions
        )
    }

    pub fn divide_cell(&mut self, mom_index: CellIndex) -> &RelCell<Cell> {
        let mom = self
            .env
            .cells
            .get_cell(mom_index);
        let div_axis = self.find_division_axis(mom, self.cell_search_scaler);
        let new_positions: Vec<_> = self
            .search_cell_box(mom, self.cell_search_scaler)
            .into_iter()
            .filter(|pos| {
                let y = div_axis.slope * pos.x as f32 + div_axis.intercept;
                (pos.y as f32) < y
            })
            .collect();
        
        let mut newborn = mom.birth();
        newborn.ancestor = Some(mom_index);
        let newborn_ta = mom.newborn_target_area;
        let new_index = self.env.cells.add(newborn).index;
        for pos in new_positions {
            self.grant_position(
                pos,
                Spin::Some(new_index),
            );
        }
        self.env.cells.get_cell_mut(mom_index).set_target_area(newborn_ta);
        self.cells.get_cell(new_index)
    }

    // Should this also replace some of the cell's positions with Medium?
    pub fn kill_cell(&mut self, cell: &mut RelCell<Cell>) {
        cell.apoptosis();
    }
    
    pub fn invalidate_cell(&mut self, cell_index: CellIndex) {
        for index in 0..self.clones_table.length() {
            self.clones_table[(index, cell_index as usize)] = false;
        }
    }

    // With some unsafe code we can return Vec<&RelCell> from this function, but it would
    // require that self.divide_cell never invalidates any references to self.cells
    // we need thorough testing of self.divide_cells to make this change, and the performance
    // gain is minimal (although the ergonomic gains are significant)
    pub fn reproduce(&mut self, rng: &mut impl Rng) {
        let mut divide = vec![];
        for cell in self.cells.iter() {
            if !cell.is_alive() {
                continue;
            }
            // Currently cells don't need to express the dividing type to divide, they just need to be big enough
            if cell.area() >= cell.divide_area() {
                divide.push(cell.index);
            }
        }
        for cell_index in divide {
            if !self.can_add_cell() {
                return;
            }

            let mom = self
                .cells
                .get_cell(cell_index);
            let new_cell = self.divide_cell(mom.index);
            if new_cell.is_valid() {
                let new_index = new_cell.index;
                self.update_clones(new_index, cell_index);
                // We could also instead choose to mutate at a fix rate throughout the cell's life cycle
                self.env.cells.get_cell_mut(new_index).genome.attempt_mutate(rng);
            }
        }
    }

    pub fn update_clones(
        &mut self,
        cell_index: CellIndex,
        mom_index: CellIndex
    ) {
        let cell = self.cells.get_cell(cell_index);
        let cell_neighs = self.cell_neighbours(cell, 2.0);

        let mom_cell = self.cells.get_cell(mom_index);
        let mom_neighs = self.cell_neighbours(
            mom_cell,
            2.0
        );

        let mom_clones = HashSet::from_iter(
            (0..self.cells.n_cells())
                .filter_map(|index| {
                    if self.clones_table[(mom_cell.index as usize, index as usize)] {
                        Some(Spin::Some(index))
                    } else {
                        None
                    }
                })
        );
        for spin in mom_clones.difference(&mom_neighs) {
            if let Spin::Some(index) = spin {
                self.clones_table[(mom_index as usize, *index as usize)] = false;
            }
        }
        let clones: Vec<_> = cell_neighs.intersection(&mom_clones).collect();
        for spin in &clones {
            if let Spin::Some(index) = spin {
                self.clones_table[(cell_index as usize, *index as usize)] = true;
            }
        }
        self.clones_table[(cell_index as usize, mom_index as usize)] = true;
    }

    // TODO!: add plot to make sure this is right
    /// Find the minor axis along which to split the cell.
    pub fn find_division_axis(&self, cell: &RelCell<Cell>, search_scaler: f32) -> SplitLine {
        // Compute covariance elements relative to centroid
        let mut sum_xx = 0.0;
        let mut sum_yy = 0.0;
        let mut sum_xy = 0.0;

        for p in &self.search_cell_box(cell, search_scaler) {
            let (dx, dy) = self.bounds.boundary.displacement(p.to_f32(), cell.center);
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
        let intercept = cell.center.y - slope * cell.center.x;

        SplitLine { slope, intercept }
    }
    
    pub fn wipe_out(&mut self) {
        self.env.wipe_out();
        self.clones_table.clear();
    }

    pub fn make_border(
        &mut self,
        bottom: bool,
        top: bool,
        left: bool,
        right: bool,
    ) {
        let mut border_positions = Vec::<Pos<usize>>::new();
        if bottom {
            for x in 0..self.width() {
                border_positions.push((x, 0).into());
            }
        }
        if top {
            for x in (0..self.width() - 1).rev() {
                border_positions.push((x, self.height() - 1).into());
            }
        }
        if left {
            for y in (1..self.height() - 1).rev() {
                border_positions.push((0, y).into());
            }
        }
        if right {
            for y in 1..self.height() {
                border_positions.push((self.width() - 1, y).into());
            }
        }

        self.spawn_solid(border_positions.into_iter());
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
            self.act_lattice[pos] = self.act_max;
        } else {
            self.act_lattice[pos] = 0;
        }
        if let Spin::Some(index) = self.cell_lattice[pos] {
            let from_cell = self.env.cells.get_cell_mut(index);
            from_cell.shift_position(pos, false, &self.env.bounds.boundary);
            from_cell.shift_chem(pos, chem_at_pos, false, &self.env.bounds.boundary);
            // If the copy kills the cell
            if from_cell.area() == 0 {
                from_cell.apoptosis();
                self.invalidate_cell(index);
            }
        }
        // Executes the copy
        self.cell_lattice[pos] = to;
        self.update_edges(pos)
    }

    fn spawn_cell(
        &mut self,
        empty_cell: Self::Cell,
        positions: impl IntoIterator<Item = Pos<usize>>
    ) -> &RelCell<Self::Cell> {
        let cell_index = self.cells.add(empty_cell).index;
        let new_spin = Spin::Some(cell_index);
        for pos in positions {
            self.grant_position(pos, new_spin);
        }
        self.cells.get_cell_mut(cell_index).ancestor =  Some(cell_index);
        self.cells.get_cell(cell_index)
    }
}

/// Represents a split line through the centroid:
///   y = slope * x + intercept
#[derive(Debug)]
pub struct SplitLine {
    pub slope: f32,
    pub intercept: f32,
}