pub use nalgebra::{RealField, Scalar};

pub type Matrix3<T = f64> = nalgebra::Matrix3<T>;
pub type Vector3<T = f64> = nalgebra::Vector3<T>;
pub type Vector6<T = f64> = nalgebra::Vector6<T>;

#[derive(Debug, Clone, thiserror::Error)]
pub enum MatrixError {
    #[error("Try to inverse a singular matrix")]
    InverseSingular,
    #[error("Try to normalize a zero matrix")]
    NormalizeZero,
}

mod orbit_coefficient;
pub use orbit_coefficient::*;

mod pos_vel;
pub use pos_vel::*;

pub mod vector3;

#[cfg(test)]
mod test;
