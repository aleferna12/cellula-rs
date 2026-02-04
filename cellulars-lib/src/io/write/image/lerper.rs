use num::{Num, Zero};
use palette::Mix;

/// Used to interpolate colors with [`Lerper::lerp()`].
#[derive(Clone, Debug, PartialEq, Eq, Hash)]
pub struct Lerper<C> {
    pub min_color: C,
    pub max_color: C
}

impl<C> Lerper<C> {
    /// Linearly interpolates `value` between `min` and `max`.
    pub fn lerp(&self, value: C::Scalar, min: C::Scalar, max: C::Scalar) -> Result<C, LerpError>
    where
        C: Mix + Clone,
        C::Scalar: Num + PartialOrd + Clone {
        if max < min {
            return Err(LerpError::NegativeRange);
        }
        if value < min {
            return Err(LerpError::ValueTooSmall);
        }
        if value > max {
            return Err(LerpError::ValueTooLarge);
        }

        let p = if min == max { C::Scalar::zero() } else { (value - min.clone()) / (max - min) };
        let blended = self.min_color.clone().mix(self.max_color.clone(), p);
        Ok(blended)
    }
}

/// Error thrown when linear interpolation fails.
#[derive(thiserror::Error, Debug)]
pub enum LerpError {
    #[error("value falls outside the range because it's too small")]
    ValueTooSmall,
    #[error("value falls outside the range because it's too large")]
    ValueTooLarge,
    #[error("minimum value passed is larger than maximum")]
    NegativeRange
}