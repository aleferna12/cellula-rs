use std::f32::consts::TAU;

/// 2D position in space.
#[derive(PartialEq, Eq, Copy, Clone, Debug)]
#[derive(Hash)]
pub struct Pos<T> {
    pub x: T,
    pub y: T
}

impl<T> Pos<T> {
    pub fn new(x: T, y: T) -> Self {
        Self { x, y }
    }
}

impl<T> From<(T, T)> for Pos<T> {
    fn from(value: (T, T)) -> Self {
        Pos::<T>::new(value.0, value.1)
    }
}

impl Pos<f32> {
    /// Unwraps the `AngularProjection` into a position.
    pub(crate) fn from_projection(proj: &AngularProjection, width: usize, height: usize) -> Self {
        let (angle_x, angle_y) = proj.angles();
        Self {
            x: width as f32 * angle_x / TAU,
            y: height as f32 * angle_y / TAU
        }
    }
}

impl Pos<usize> {
    pub(crate) fn pack_u32(self) -> u32 {
        ((self.x as u32) << 16) | self.y as u32
    }

    pub fn row_major(self, height: usize) -> usize {
        self.x * height + self.y
    }
}

impl From<Pos<usize>> for Pos<isize> {
    fn from(value: Pos<usize>) -> Self {
        Pos::new(value.x as isize, value.y as isize)
    }
}

impl From<Pos<isize>> for Pos<usize> {
    fn from(value: Pos<isize>) -> Self {
        let message = "overflow when translating position from general to lattice coordinates";
        Pos::new(
            value.x.try_into().expect(message), 
            value.y.try_into().expect(message)
        )
    }
}

#[derive(Debug, Clone)]
pub struct WrappedPos {
    pub(crate) pos: Pos<f32>,
    pub(crate) projection: AngularProjection
}

impl WrappedPos {
    pub fn new(pos: Pos<f32>, width: usize, height: usize) -> Self {
        Self {
            pos,
            projection: AngularProjection::from_pos(pos, width, height)
        }
    }

    /// Represents the origin of the lattice, at 0, 0.
    pub fn origin() -> Self {
        Self {
            pos: (0., 0.).into(),
            projection: AngularProjection {
                x_sin: 0.,
                x_cos: 1.,
                y_sin: 0.,
                y_cos: 1.,
            }
        }
    }

    pub fn pos(&self) -> Pos<f32> {
        self.pos
    }
}

#[derive(Debug, Clone)]
pub(crate) struct AngularProjection {
    pub(crate) x_sin: f32,
    pub(crate) x_cos: f32,
    pub(crate) y_sin: f32,
    pub(crate) y_cos: f32
}

impl AngularProjection {
    pub(crate) fn from_pos(pos: Pos<f32>, width: usize, height: usize) -> Self {
        let cx = TAU * pos.x / width as f32;
        let cy = TAU * pos.y  / height as f32;
        Self {
            x_sin: cx.sin(),
            x_cos: cx.cos(),
            y_sin: cy.sin(),
            y_cos: cy.cos()
        }
    }
    
    /// Returns the angles associated with this projection.
    pub(crate) fn angles(&self) -> (f32, f32) {
        (self.x_sin.atan2(self.x_cos).rem_euclid(TAU), self.y_sin.atan2(self.y_cos).rem_euclid(TAU))
    }
    
    // TODO!: this can be optimised to not require atan2 
    //  (and has a significant impact in performance due to hot call in CA)
    pub(crate) fn delta_angles(&self, other: &AngularProjection) -> (f32, f32) {
        // This avoids unecessary multiple calls to rem_euclid
        let x = (self.x_sin * other.x_cos - self.x_cos * other.x_sin)
            .atan2(self.x_cos * other.x_cos + self.x_sin * other.x_sin);
        let y = (self.y_sin * other.y_cos - self.y_cos * other.y_sin)
            .atan2(self.y_cos * other.y_cos + self.y_sin * other.y_sin);
        (x, y)
    }
}