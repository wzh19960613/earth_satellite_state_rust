pub use nalgebra::{ComplexField, Matrix3, RealField, Scalar, Vector3, Vector6};

#[derive(Debug, thiserror::Error)]
#[error("Try to inverse a singular matrix")]
pub struct InverseOfSingularMatrix;

mod impl_isomorphic;
use impl_isomorphic::*;

mod orbit_coefficient;
pub use orbit_coefficient::*;

mod pos_vel_state;
pub use pos_vel_state::*;

pub mod vector3;

#[cfg(test)]
mod test;
