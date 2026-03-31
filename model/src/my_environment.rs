use std::f64::consts::E;
use crate::cell::Cell;
use crate::constants::{BoundaryType, EPSILON};
use bon::bon;
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
use std::ops::{Deref, DerefMut};

#[derive(Clone)]
pub struct MyEnvironment {
    env: Environment<Cell, MooreNeighbourhood, BoundaryType>,
    pub chem_lattice: Lattice<u32>,
    pub act_lattice: Lattice<u32>,
    pub cell_search_scaler: f32,
    pub max_cells: CellIndex,
    pub act_max: u32,
    pub max_chem: u32,
    population_exploded: bool,
    corners: [Pos<usize>; 4],
    current_corner: usize,
}

#[bon]
impl MyEnvironment {
    #[builder]
    pub fn new(
        env: Environment<Cell, MooreNeighbourhood, BoundaryType>,
        max_cells: CellIndex,
        act_max: u32,
        cell_search_scaler: f32
    ) -> Self {
        let lat = Lattice::new(env.cell_lattice.rect.clone());
        Self::new_from_backup()
            .env(env)
            .chem_lattice(lat.clone())
            .act_lattice(lat)
            .max_cells(max_cells)
            .act_max(act_max)
            .cell_search_scaler(cell_search_scaler)
            .call()
    }
    
    #[builder]
    pub fn new_from_backup(
        env: Environment<Cell, MooreNeighbourhood, BoundaryType>,
        chem_lattice: Lattice<u32>,
        act_lattice: Lattice<u32>,
        max_cells: CellIndex,
        act_max: u32,
        cell_search_scaler: f32
    ) -> Self {
        Self {
            max_chem: Self::distance(
                (0, 0).into(),
                (env.cell_lattice.width(), env.cell_lattice.height()).into()
            ).round() as u32 + 1,
            population_exploded: false,
            corners: [
                (0, 0).into(), (env.width() - 1, 0).into(),
                (0, env.height() - 1).into(),
                (env.width() - 1, env.height() - 1).into()
            ],
            current_corner: 0,
            env,
            cell_search_scaler,
            max_cells,
            act_max,
            chem_lattice,
            act_lattice
        }
    }

    pub fn chem_center(&self) -> Pos<usize> {
        self.corners[self.current_corner]
    }

    fn distance2(pos1: Pos<usize>, pos2: Pos<usize>) -> usize {
        (pos1.x as isize - pos2.x as isize).pow(2) as usize
            + (pos1.y as isize - pos2.y as isize).pow(2) as usize
    }

    fn distance(pos1: Pos<usize>, pos2: Pos<usize>) -> f32 {
        (Self::distance2(pos1, pos2) as f32).sqrt()
    }

    pub fn chem_signal(&self, pos: Pos<usize>, chem_center: Pos<usize>) -> u32 {
        self.max_chem - Self::distance(pos, chem_center).round() as u32
    }

    pub fn update_rel_chem(&mut self) {
        let mut min = f32::INFINITY;
        let mut max = 0.;
        for rel_cell in self.env.cells.iter() {
            if !rel_cell.is_valid() {
                continue;
            }
            let cell_chem = rel_cell.cell.chem_mass as f32 / rel_cell.cell.area as f32;
            if cell_chem < min {
                min = cell_chem;
            }
            if cell_chem > max {
                max = cell_chem;
            }
        }
        for rel_cell in self.env.cells.iter_mut() {
            let cell_chem = rel_cell.cell.chem_mass as f32 / rel_cell.cell.area as f32;
            rel_cell.rel_chem = (cell_chem - min) / (max - min);
        }
    }

    pub fn update_neighbours(&mut self) {
        for rel_cell in self.cells.iter_mut() {
            rel_cell.neighbors.clear();
        }
        for edge in self.env.edge_book.iter().cloned() {
            let spin1: Spin = self.cell_lattice[edge.p1];
            let spin2: Spin = self.cell_lattice[edge.p2];
            if let Spin::Some(cell_index) = spin1 {
                let rel_cell = self.env.cells.get_cell_mut(cell_index);
                let entry = rel_cell.neighbors.entry(spin2).or_insert(0);
                *entry += 1;
            }
            if let Spin::Some(cell_index) = spin2 {
                let rel_cell = self.env.cells.get_cell_mut(cell_index);
                let entry = rel_cell.neighbors.entry(spin1).or_insert(0);
                *entry += 1;
            }
        }
    }

    pub fn update_act(&mut self) {
        let mut act_pairs = vec![];
        for cell in self.cells.iter() {
            let mut act = 0;
            let mut kact = 0.;
            for pos in self.search_cell_box(cell, self.cell_search_scaler) {
                act += self.act_lattice[pos];
                kact += self.kact(pos);
            }
            act_pairs.push((cell.index, (act, kact)));
        }
        for (index, (act, kact)) in act_pairs {
            let cell = self.cells.get_cell_mut(index);
            cell.tot_act = act;
            cell.tot_kact = kact;
        }
    }

    pub fn kact(&self, pos: Pos<usize>) -> f64 {
        let mut local_act = (self.act_lattice[pos] as f64).ln();
        let mut owned_neighs = 1;
        for neigh in self.valid_neighbours(pos) {
            if self.cell_lattice[neigh] != self.cell_lattice[pos] {
                continue;
            }
            local_act += (self.act_lattice[neigh] as f64).ln();
            owned_neighs += 1;
        }
        E.powf(local_act / owned_neighs as f64)
    }
    
    pub fn reset_act(&mut self) {
        for cell in self.cells.iter_mut() {
            cell.tot_act = 0;
            cell.tot_kact = 0.;
        }
    }

    fn make_chem_gradient(&mut self, chem_center: Pos<usize>) {
        for pos in self.chem_lattice.iter_positions() {
            let new_chem = self.chem_signal(pos, chem_center);
            if let Spin::Some(cell_index) = self.cell_lattice[pos] {
                let prev_chem = self.chem_lattice[pos];
                let cell = self.env.cells.get_cell_mut(cell_index);
                // Remove position from center_chem/chem_mass
                cell.shift_chem(
                    pos,
                    prev_chem,
                    false,
                    &self.env.bounds.boundary
                );
                // Add position to center_chem/chem_mass
                cell.shift_chem(
                    pos,
                    new_chem,
                    true,
                    &self.env.bounds.boundary
                );
            }
            self.chem_lattice[pos] = new_chem;
        }
    }

    pub fn make_next_chem_gradient(&mut self, rng: &mut impl Rng) {
        let curr_corner = self.current_corner;
        while self.current_corner == curr_corner {
            self.current_corner = rng.random_range(0..self.corners.len());
        }
        self.make_chem_gradient(self.corners[curr_corner]);
    }

    pub fn can_add_cell(&mut self) -> bool {
        if self.cells.n_valid() < self.max_cells {
            return true;
        }
        if !self.population_exploded {
            log::warn!(
                        "Population exceeded maximum threshold `max-cells={}`, cell is not going to be added",
                        {self.max_cells}
                    );
            log::warn!("This warning will be suppressed from now on");
            self.population_exploded = true;
        }
        false
    }

    /// Spawns an `empty_cell` on valid `positions` that belong to the medium,
    /// while ignoring positions owned by solids or other cells.
    pub fn spawn_cell_checked(
        &mut self,
        empty_cell: Cell,
        positions: impl IntoIterator<Item = Pos<isize>>
    ) -> &RelCell<Cell> {
        let med_positions = positions.into_iter().filter_map(|pos| {
            let valid_pos = self.bounds.lattice_boundary.valid_pos(pos)?;
            let lat_pos = valid_pos.to_usize();
            if !matches!(self.cell_lattice[lat_pos], Spin::Medium) {
                return None;
            }
            Some(lat_pos)
        }).collect::<Box<[_]>>();
        self.spawn_cell(empty_cell, med_positions)
    }

    // TODO: make spawn as a circle with center at pos
    // TODO: make spawn as a circle with center at pos
    /// Spawns a square cell centered at a random position with area = `cell_area`.
    ///
    /// Uses [Environment::spawn_cell_checked()] to restrict spawns to medium.
    pub fn spawn_cell_random(
        &mut self,
        empty_cell: Cell,
        cell_area: u32,
        rng: &mut impl Rng,
    ) -> &RelCell<Cell> {
        let pos_isize = self
            .cell_lattice
            .random_pos(rng)
            .to_isize();
        let cell_side = ((cell_area as f32).sqrt() / 2.).floor() as isize;
        let rect = Rect::new(
            Pos::new(pos_isize.x - cell_side, pos_isize.y - cell_side),
            Pos::new(pos_isize.x + cell_side, pos_isize.y + cell_side)
        );
        self.spawn_cell_checked(
            empty_cell,
            rect.iter_positions()
        )
    }

    pub fn divide_cell(&mut self, mom_index: CellIndex) -> &RelCell<Cell> {
        let mom = self
            .env
            .cells
            .get_cell(mom_index);
        let div_axis = self.find_division_axis(mom, self.cell_search_scaler);
        let new_positions: Box<_> = self
            .search_cell_box(mom, self.cell_search_scaler)
            .into_iter()
            .filter(|pos| {
                let y = div_axis.slope * pos.x as f32 + div_axis.intercept;
                (pos.y as f32) < y
            })
            .collect();

        let new_index = self.env.cells.add(mom.birth()).index;
        for pos in new_positions {
            let prev_act = self.act_lattice[pos];
            self.update_delta_perimeter(false, mom_index, pos);
            self.update_delta_perimeter(true, new_index, pos);
            self.grant_position(
                pos,
                Spin::Some(new_index),
            );
            // We dont want to model this position transfer as a copy
            self.act_lattice[pos] = prev_act;
        }
        self.cells.get_cell(new_index)
    }

    pub fn erase_cell(&mut self, cell_index: CellIndex) {
        let cell = self.env.cells.get_cell(cell_index);
        // It's imperative to find all positions here as otherwise we are left over with a partial cell
        let mut cell_positions = self.search_cell_box(cell, self.cell_search_scaler);
        if cell.area as usize != cell_positions.len() {
            cell_positions = self.search_cell_box(cell, 2. * self.cell_search_scaler);
        }
        if cell.area as usize != cell_positions.len() {
            log::error!("Critical: cell with index {} was only partially erased after two attempts", cell.index);
        }
        for pos in cell_positions {
            self.update_delta_perimeter(false, cell_index, pos);
            self.grant_position(pos, Spin::Medium);
        }
        self.env.cells.get_cell_mut(cell_index).apoptosis();
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

    /// Updates the perimeter to be shifted on the cell corresponding to `cell_index` on the next call to [MyEnvironment::grant_position()].
    /// 
    /// This function should be called exactly once on the source and the target of the site copy, before the call to [MyEnvironment::grant_position()].   
    pub fn update_delta_perimeter(
        &mut self,
        source: bool,
        cell_index: CellIndex,
        pos: Pos<usize>
    ) {
        let shift_when_eq = if source { -1 } else { 1 };
        let cell_spin = Spin::Some(cell_index);
        self.cells.get_cell_mut(cell_index).delta_perimeter = Some(self
            .neighbour_spins(pos)
            .map(|spin| if spin == cell_spin { shift_when_eq } else { -shift_when_eq } )
            .sum()
        );
    }

    // This is slightly faster than allocating a Vec that can be reutilized
    pub fn neighbour_spins(&self, pos: Pos<usize>) -> impl Iterator<Item = Spin> {
        self.valid_neighbours(pos).map(|neigh| {
            self.cell_lattice[neigh]
        })
    }
}

impl Deref for MyEnvironment {
    type Target = Environment<Cell, MooreNeighbourhood, BoundaryType>;
    fn deref(&self) -> &Self::Target {
        &self.env
    }
}

impl DerefMut for MyEnvironment {
    fn deref_mut(&mut self) -> &mut Self::Target {
        &mut self.env
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
            to_cell.delta_perimeter = None;
            self.act_lattice[pos] = self.act_max;
        } else {
            self.act_lattice[pos] = 0;
        }
        if let Spin::Some(index) = self.cell_lattice[pos] {
            let from_cell = self.env.cells.get_cell_mut(index);
            from_cell.shift_position(pos, false, &self.env.bounds.boundary);
            from_cell.shift_chem(pos, chem_at_pos, false, &self.env.bounds.boundary);
            from_cell.delta_perimeter = None;
            // If the copy kills the cell
            if from_cell.area() == 0 {
                from_cell.apoptosis();
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
            if let Spin::Some(target_index) = self.cell_lattice[pos] {
                self.update_delta_perimeter(false, target_index, pos);
            }
            self.update_delta_perimeter(true, cell_index, pos);
            self.grant_position(pos, new_spin);
        }
        self.cells.get_cell_mut(cell_index).ancestor =  Some(cell_index);
        self.cells.get_cell(cell_index)
    }

    fn spawn_solid(&mut self, positions: impl Iterator<Item = Pos<usize>>) {
        for pos in positions {
            if let Spin::Some(target_index) = self.cell_lattice[pos] {
                self.update_delta_perimeter(false, target_index, pos);
            }
            self.grant_position(pos, Spin::Solid);
        }
    }
}

/// Represents a split line through the centroid:
///   y = slope * x + intercept
#[derive(Debug)]
pub struct SplitLine {
    pub slope: f32,
    pub intercept: f32,
}
