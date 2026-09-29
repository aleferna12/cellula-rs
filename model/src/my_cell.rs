//! Contains logic associated with [`MyCell`].

use crate::constants::EPSILON;
use bon::Builder;
use cellulars_lib::base::cell::Cell;
use cellulars_lib::constants::FloatType;
use cellulars_lib::positional::boundaries::Boundary;
use cellulars_lib::positional::com::Com;
use cellulars_lib::positional::pos::Pos;
use cellulars_lib::traits::cellular::{Alive, Cellular, EmptyCell};
use cellulars_lib::traits::track_perimeter::TrackPerimeter;
use rand::Rng;
use strum_macros::{Display, EnumString};
#[cfg(not(feature = "high-precision"))]
use std::f32::consts::TAU;
#[cfg(feature = "high-precision")]
use std::f64::consts::TAU;

/// A cell that can track a chemical concentration and migrate towards its source.
#[derive(Clone, Debug, Builder)]
pub struct MyCell {
    /// Area at which the cell divides.
    pub divide_area: u32,
    /// Target area for newborns of this cell (see [`Alive::birth()`]).
    pub newborn_target_area: u32,
    /// Perimeter the cell strives for
    /// (see [`PerimeterConstraint`](cellulars_lib::perimeter_constraint::PerimeterConstraint)).
    pub target_perimeter: u32,
    /// For how many time-steps the cell keeps pushing in the same direction before realigning itself
    /// with the direction it actually travelled (see [`MyCell::update_persistence()`]).
    ///
    /// Zero disables persistent migration.
    pub persistence_duration: u32,
    /// Current type of the cell.
    pub cell_type: CellType,
    /// Inner base cell.
    pub cell: Cell,
    /// Center of mass of the cell's perceived chemical concentration.
    chem_com: Com,
    /// Kept up to date by
    /// [`MyEnvironment::grant_position()`](crate::my_environment::MyEnvironment::grant_position()),
    /// which is the only place that knows the lattice around the cell.
    #[builder(default)]
    perimeter: u32,
    /// Unit vector along which the cell pushes its protrusions.
    #[builder(default = Pos::new(0., 0.))]
    target_vec: Pos<FloatType>,
    /// Center of the cell the last time [`MyCell::update_target_vec()`] ran.
    #[builder(default = Pos::new(0., 0.))]
    prev_center: Pos<FloatType>,
    /// Counts the time-steps elapsed since the last realignment of [`MyCell::target_vec()`].
    #[builder(default)]
    persistence_time: u32
}

impl MyCell {
    /// Initialises an empty [`MyCell`] to be filled progressively with [`MyCell::shift_position()`].
    pub fn new_empty(
        target_area: u32,
        target_perimeter: u32,
        divide_area: u32,
        persistence_duration: u32,
        cell_type: CellType
    ) -> EmptyCell<Self> {
        EmptyCell::new(Self {
            cell: Cell::new_empty(target_area).into_cell(),
            chem_com: Com { pos: Pos::new(0., 0.), mass: 0 },
            newborn_target_area: target_area,
            perimeter: 0,
            target_vec: Pos::new(0., 0.),
            prev_center: Pos::new(0., 0.),
            persistence_time: 0,
            target_perimeter,
            persistence_duration,
            divide_area,
            cell_type,
        }).expect("cell was not empty")
    }
    
    /// Returns the total concentration of the chemical perceived by the cell.
    pub fn chem_mass(&self) -> u32 {
        self.chem_com.mass
    }

    /// Returns the center of the cell weighted by the chemical concentration at each cell position.
    pub fn chem_center(&self) -> Pos<FloatType> {
        self.chem_com.pos
    }

    /// Sets the area at which the cell divides when
    /// [`Environment::reproduce()`](crate::my_environment::MyEnvironment::reproduce()) is called.
    pub fn set_divide_area(&mut self, value: u32) {
        self.divide_area = value;
    }

    /// Adds or removes the chemical concentration `chem_at` at position `pos` from the cell.
    pub fn shift_chem<B: Boundary<Coord = FloatType>>(&mut self, pos: Pos<usize>, chem_at: u32, adding: bool, boundary: &B) {
        let shifted = self.chem_com.shift(
            Com { pos: pos.cast_as(), mass: chem_at },
            adding,
            boundary
        );
        match shifted {
            Ok(new_com) => self.chem_com = new_com,
            Err(e) => log::warn!("Failed to shift chem center: {e}")
        }
    }

    /// Returns the unit vector along which the cell pushes its protrusions.
    ///
    /// Is the zero vector for cells that have not been given a direction yet
    /// (see [`MyCell::has_direction()`]).
    pub fn target_vec(&self) -> Pos<FloatType> {
        self.target_vec
    }

    /// Returns the center the cell had the last time [`MyCell::update_target_vec()`] ran.
    pub fn prev_center(&self) -> Pos<FloatType> {
        self.prev_center
    }

    /// Returns how many time-steps have elapsed since the last realignment of [`MyCell::target_vec()`].
    pub fn persistence_time(&self) -> u32 {
        self.persistence_time
    }

    /// Returns whether the cell has a direction to migrate towards.
    ///
    /// Is `false` for cells that were never given one by [`MyCell::randomize_direction()`], such as
    /// freshly spawned cells or cells restored from data files written before persistent migration existed.
    pub fn has_direction(&self) -> bool {
        self.target_vec.x.hypot(self.target_vec.y) > EPSILON
    }

    /// Points the cell in a uniformly random direction and staggers its persistence clock, so that cells
    /// initialised together do not all turn on the same time-step.
    ///
    /// Ports `Cell::startTarVec()` and the `perstime` initialisation of
    /// [Colizzi, 2020](https://doi.org/10.7554/eLife.56349).
    pub fn randomize_direction(&mut self, rng: &mut impl Rng) {
        let angle = rng.random::<FloatType>() * TAU;
        self.target_vec = Pos::new(angle.cos(), angle.sin());
        self.prev_center = self.center();
        self.persistence_time = if self.persistence_duration == 0 {
            0
        } else {
            rng.random_range(0..self.persistence_duration)
        };
    }

    /// Realigns [`MyCell::target_vec()`] with the direction the cell actually travelled since this method
    /// last ran.
    ///
    /// A cell that did not move keeps pushing in the same direction.
    pub fn update_target_vec(&mut self, boundary: &impl Boundary<Coord = FloatType>) {
        let (dx, dy) = boundary.displacement(self.prev_center, self.center());
        let hyp = dx.hypot(dy);
        if hyp > EPSILON {
            self.target_vec = Pos::new(dx / hyp, dy / hyp);
        }
        self.prev_center = self.center();
    }

    /// Advances the cell's persistence clock, realigning [`MyCell::target_vec()`] once every
    /// [`MyCell::persistence_duration`] time-steps.
    ///
    /// Since the realized direction of motion is a noisy version of the direction the cell was pushing
    /// towards, a longer `persistence_duration` makes the cell turn more slowly.
    ///
    /// Ports `Cell::updatePersTime()` of [Colizzi, 2020](https://doi.org/10.7554/eLife.56349).
    pub fn update_persistence(&mut self, boundary: &impl Boundary<Coord = FloatType>) {
        if self.persistence_duration == 0 {
            return;
        }
        self.persistence_time += 1;
        if self.persistence_time >= self.persistence_duration {
            self.update_target_vec(boundary);
            self.persistence_time = 0;
        }
    }

    /// Updates parameters of the cell (called by [`Pond::step()`](cellulars_lib::traits::step::Step::step())).
    pub fn update(&mut self) {
        if let CellType::Dividing = self.cell_type && self.target_area() < self.divide_area {
            let new_target_area = self.target_area() + 1;
            self.cell.target_area = new_target_area;
        }
    }
}

impl Cellular for MyCell {
    fn target_area(&self) -> u32 {
        self.cell.target_area()
    }

    fn area(&self) -> u32 {
        self.cell.area()
    }

    fn center(&self) -> Pos<FloatType> {
        self.cell.center()
    }

    fn is_empty(&self) -> bool {
        self.cell.is_empty()
    }

    fn shift_position(&mut self, pos: Pos<usize>, adding: bool, bound: &impl Boundary<Coord = FloatType>) {
        self.cell.shift_position(pos, adding, bound)
    }
}

impl TrackPerimeter for MyCell {
    fn perimeter(&self) -> u32 {
        self.perimeter
    }

    fn target_perimeter(&self) -> u32 {
        self.target_perimeter
    }

    fn shift_perimeter(&mut self, delta: i32) {
        self.perimeter = self.perimeter.saturating_add_signed(delta);
    }
}

impl Alive for MyCell {
    fn is_alive(&self) -> bool {
        self.cell.is_alive()
    }

    fn apoptosis(&mut self) {
        self.cell.apoptosis()
    }

    fn birth(&self) -> EmptyCell<Self> {
        let mut basic_cell = self.cell.birth().into_cell();
        basic_cell.target_area = self.newborn_target_area;
        EmptyCell::new(Self {
            chem_com: Com { pos: basic_cell.center(), mass: 0 },
            prev_center: basic_cell.center(),
            perimeter: 0,
            persistence_time: 0,
            cell: basic_cell,
            ..self.clone()
        }).expect("failed to create empty cell")
    }
}

/// A cell is either migrating or dividing.
#[derive(Clone, Debug, EnumString, Display)]
#[strum(serialize_all = "kebab-case")]
pub enum CellType {
    /// A cell that is migrating.
    Migrating,
    /// A cell that is dividing.
    Dividing
}

#[cfg(test)]
mod tests {
    use super::*;
    use cellulars_lib::positional::boundaries::UnsafePeriodicBoundary;
    use cellulars_lib::positional::rect::Rect;
    use rand::SeedableRng;
    use rand_xoshiro::Xoshiro256StarStar;

    fn make_unsafe_boundary() -> UnsafePeriodicBoundary<FloatType> {
        UnsafePeriodicBoundary::new(Rect::new((0., 0.).into(), (100., 100.).into()))
    }
    
    fn make_test_cell() -> MyCell {
        MyCell::new_empty(
            100,
            120,
            200,
            10,
            CellType::Migrating,
        ).into_cell()
    }

    #[test]
    fn test_shift_position_area_and_center() {
        let mut cell = make_test_cell();
        let bound = make_unsafe_boundary();

        cell.shift_position(Pos::new(10, 10), true, &bound);
        assert_eq!(cell.area(), 1);
        assert_eq!(cell.center(), Pos::new(10.0, 10.0));

        cell.shift_position(Pos::new(20, 20), true, &bound);
        assert_eq!(cell.area(), 2);
        assert_eq!(cell.center(), Pos::new(15.0, 15.0));

        cell.shift_position(Pos::new(10, 10), false, &bound);
        assert_eq!(cell.area(), 1);
        assert_eq!(cell.center(), Pos::new(20.0, 20.0));
    }

    #[test]
    fn test_shift_position_chem_center_and_mass() {
        let bound = make_unsafe_boundary();
        let mut cell = make_test_cell();

        // Add chem at (2, 3) with value 10
        cell.shift_chem(Pos::new(2, 3), 10, true, &bound);
        assert_eq!(cell.chem_com.mass, 10);
        assert_eq!(cell.chem_com.pos, Pos::new(2., 3.));

        // Add chem at (4, 5) with value 10
        cell.shift_chem(Pos::new(4, 5), 10, true, &bound);
        assert_eq!(cell.chem_com.mass, 20);
        assert_eq!(cell.chem_com.pos, Pos::new(3., 4.));

        // Remove chem from (2, 3)
        cell.shift_chem(Pos::new(2, 3), 10, false, &bound);
        assert_eq!(cell.chem_com.mass, 10);
        assert_eq!(cell.chem_com.pos, Pos::new(4., 5.));
    }

    #[test]
    fn test_randomize_direction_gives_a_unit_vector() {
        let mut rng = Xoshiro256StarStar::seed_from_u64(42);
        let mut cell = make_test_cell();
        assert!(!cell.has_direction());

        for _ in 0..32 {
            cell.randomize_direction(&mut rng);
            assert!(cell.has_direction());
            let hyp = cell.target_vec().x.hypot(cell.target_vec().y);
            assert!((hyp - 1.).abs() < EPSILON, "target vector is not a unit vector");
            // Cells start with a staggered clock so they don't all turn on the same time-step
            assert!(cell.persistence_time() < cell.persistence_duration);
        }
    }

    #[test]
    fn test_target_vec_follows_the_realized_displacement() {
        let bound = make_unsafe_boundary();
        let mut cell = make_test_cell();
        cell.shift_position(Pos::new(10, 10), true, &bound);
        // prev_center is only set by `randomize_direction` and `update_target_vec`
        cell.update_target_vec(&bound);
        assert_eq!(cell.prev_center(), Pos::new(10., 10.));

        // Moving the cell one position to the right points it along +x
        cell.shift_position(Pos::new(11, 10), true, &bound);
        cell.shift_position(Pos::new(10, 10), false, &bound);
        cell.update_target_vec(&bound);
        assert_eq!(cell.target_vec(), Pos::new(1., 0.));
        assert_eq!(cell.prev_center(), Pos::new(11., 10.));
    }

    #[test]
    fn test_persistence_only_turns_the_cell_every_duration_steps() {
        let bound = make_unsafe_boundary();
        let mut cell = make_test_cell();
        cell.shift_position(Pos::new(10, 10), true, &bound);
        cell.target_vec = Pos::new(0., 1.);
        cell.prev_center = cell.center();

        // Move the cell along +x while its target vector points along +y
        cell.shift_position(Pos::new(11, 10), true, &bound);
        cell.shift_position(Pos::new(10, 10), false, &bound);
        for _ in 0..cell.persistence_duration - 1 {
            cell.update_persistence(&bound);
            assert_eq!(cell.target_vec(), Pos::new(0., 1.), "cell turned before its persistence ran out");
        }
        cell.update_persistence(&bound);
        assert_eq!(cell.target_vec(), Pos::new(1., 0.));
        assert_eq!(cell.persistence_time(), 0);
    }

    #[test]
    fn test_zero_duration_disables_persistence() {
        let bound = make_unsafe_boundary();
        let mut cell = MyCell::new_empty(100, 120, 200, 0, CellType::Migrating).into_cell();
        cell.shift_position(Pos::new(10, 10), true, &bound);
        cell.target_vec = Pos::new(0., 1.);

        cell.shift_position(Pos::new(11, 10), true, &bound);
        for _ in 0..64 {
            cell.update_persistence(&bound);
        }
        assert_eq!(cell.target_vec(), Pos::new(0., 1.));
        assert_eq!(cell.persistence_time(), 0);
    }

    #[test]
    fn test_newborns_start_with_a_clean_perimeter() {
        let bound = make_unsafe_boundary();
        let mut cell = make_test_cell();
        cell.shift_position(Pos::new(10, 10), true, &bound);
        cell.shift_perimeter(8);
        assert_eq!(cell.perimeter(), 8);

        let newborn = cell.birth().into_cell();
        assert_eq!(newborn.perimeter(), 0);
        assert_eq!(newborn.target_perimeter(), cell.target_perimeter());
    }
}
