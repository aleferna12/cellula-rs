//! Contains logic associated with [`Lerper`].

use num::{Num, One, Zero};
use palette::Mix;

/// Used to interpolate colors with [`Lerper::lerp()`].
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Lerper<C> {
    /// Color returned when the interpolated value is the minimum.
    pub min_color: C,
    /// Color returned when the interpolated value is the maximum.
    pub max_color: C
}

impl<C> Lerper<C> {
    /// Linearly interpolates colors at `p`.
    pub fn lerp(&self, p: C::Scalar) -> Result<C, LerpError>
    where
        C: Mix + Clone,
        C::Scalar: Num + PartialOrd + Clone {
        if p < C::Scalar::zero() {
            return Err(LerpError::ValueTooSmall);
        }
        if p > C::Scalar::one() {
            return Err(LerpError::ValueTooLarge);
        }

        let blended = self.min_color.clone().mix(self.max_color.clone(), p);
        Ok(blended)
    }
}

/// Error thrown when linear interpolation fails.
#[derive(thiserror::Error, Debug)]
pub enum LerpError {
    /// The value falls outside the range because it's too small.
    #[error("cannot interpolate at p < 0")]
    ValueTooSmall,

    /// The value falls outside the range because it's too large.
    #[error("cannot interpolate at p > 0")]
    ValueTooLarge
}