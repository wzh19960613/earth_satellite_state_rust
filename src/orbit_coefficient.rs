use isomorphic::{impl_mut, impl_ref, impl_transmute};

use crate::{J2000, Matrix3, RealField, Scalar, Vector3, Vector6};

#[derive(Debug, Clone, PartialEq)]
pub struct OrbitCoefficients<T: Scalar = f64> {
    /// Semi-major axis in `m`
    pub a: T,
    /// Eccentricity
    pub e: T,
    /// Inclination in `rad`
    pub i: T,
    /// Longitude of ascending node in `rad`
    pub omega: T,
    /// Argument of perigee in `rad`
    pub w: T,
    /// True anomaly in `rad`
    pub theta: T,
}

impl_ref!(OrbitCoefficients<T> = [T; 6] where T: Scalar);
impl_ref!(OrbitCoefficients<T> = Vector6<T> where T: Scalar);
impl_mut!(OrbitCoefficients<T> = [T; 6] where T: Scalar);
impl_mut!(OrbitCoefficients<T> = Vector6<T> where T: Scalar);
impl_transmute!(OrbitCoefficients<T> = [T; 6] where T: Scalar);
impl_transmute!(OrbitCoefficients<T> = Vector6<T> where T: Scalar);

pub trait CanBeOrbitCoefficient: RealField + Copy {
    /// Gravitational constant in `m^3/s^2`
    const MU: Self;
}

impl CanBeOrbitCoefficient for f32 {
    const MU: Self = 3.9860045e14;
}

impl CanBeOrbitCoefficient for f64 {
    const MU: Self = 3.986004418e14;
}

impl<T: CanBeOrbitCoefficient> From<OrbitCoefficients<T>> for J2000<T> {
    fn from(orbit: OrbitCoefficients<T>) -> Self {
        fn cos_sin<T: CanBeOrbitCoefficient>(rad: T) -> (T, T) {
            (rad.cos(), rad.sin())
        }
        let e = orbit.e;
        let h = (T::MU * orbit.a * (T::one() - e * e)).sqrt();
        let (cos_i, sin_i) = cos_sin(orbit.i);
        let (cos_omega, sin_omega) = cos_sin(orbit.omega);
        let (cos_w, sin_w) = cos_sin(orbit.w);
        let (cos_theta, sin_theta) = cos_sin(orbit.theta);
        let r_factor = (h * h / T::MU) / (T::one() + e * cos_theta);
        let a = Matrix3::<T>::from_row_slice(&[
            cos_omega * cos_w - sin_omega * sin_w * cos_i,
            -cos_omega * sin_w - sin_omega * cos_w * cos_i,
            sin_omega * sin_i,
            sin_omega * cos_w + cos_omega * sin_w * cos_i,
            -sin_omega * sin_w + cos_omega * cos_w * cos_i,
            -cos_omega * sin_i,
            sin_w * sin_i,
            cos_w * sin_i,
            cos_i,
        ]);
        let v_factor = T::MU / h;
        let zero = T::zero();
        let pos = a * Vector3::new(r_factor * cos_theta, r_factor * sin_theta, zero);
        let vel = a * Vector3::new(v_factor * (-sin_theta), v_factor * (e + cos_theta), zero);
        J2000::new(pos, vel)
    }
}
