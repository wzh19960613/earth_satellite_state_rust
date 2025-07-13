use crate::{ComplexField, RealField, Scalar};
use crate::{InverseOfSingularMatrix, J2000, Matrix3, PosVelState, Vector3};

pub fn to_array<T>(vec: Vector3<T>) -> [T; 3] {
    let [arr] = vec.data.0;
    arr
}

pub fn lvlh_from_vvlh<T: Scalar + std::ops::Neg<Output = T>>(vvlh: Vector3<T>) -> Vector3<T> {
    let [x, y, z] = to_array(vvlh);
    Vector3::new(-z, x, -y)
}

pub fn vvlh_from_lvlh<T: Scalar + std::ops::Neg<Output = T>>(lvlh: Vector3<T>) -> Vector3<T> {
    let [x, y, z] = to_array(lvlh);
    Vector3::new(y, -z, -x)
}

pub fn try_j2000_vel_from_lvlh<T: ComplexField>(
    lvlh_vel: Vector3<T>,
    ref_state: &J2000<T>,
) -> Result<Vector3<T>, InverseOfSingularMatrix> {
    let J2000(PosVelState { pos: rou_ref, vel }) = ref_state;
    let rou_hat = rou_ref.normalize();
    let h_hat = rou_ref.cross(vel).normalize();
    let theta_hat = h_hat.cross(&rou_hat);
    let a = Matrix3::from_columns(&[rou_hat, theta_hat, h_hat]).transpose();
    Ok(a.try_inverse().ok_or(InverseOfSingularMatrix)? * lvlh_vel)
}

pub fn try_j2000_vel_from_vvlh<T: ComplexField>(
    vvlh_vel: Vector3<T>,
    ref_state: &J2000<T>,
) -> Result<Vector3<T>, InverseOfSingularMatrix> {
    try_j2000_vel_from_lvlh(lvlh_from_vvlh(vvlh_vel), ref_state)
}

pub fn lvlh_vel_from_j2000<T: RealField + Copy>(
    j2000_vel: Vector3<T>,
    ref_state: &J2000<T>,
) -> Vector3<T> {
    lvlh_vel_from_j2000_impl(j2000_vel, ref_state).lvlh_vel
}

pub(crate) struct LvlhVelFromJ2000ImplResult<T> {
    pub lvlh_vel: Vector3<T>,
    pub a: Matrix3<T>,
    pub rou_rel: Vector3<T>,
}

pub(crate) fn lvlh_vel_from_j2000_impl<T: RealField + Copy>(
    j2000_vel: Vector3<T>,
    ref_state: &J2000<T>,
) -> LvlhVelFromJ2000ImplResult<T> {
    let J2000(PosVelState { pos, vel }) = ref_state;
    let (rou_ref, v_ref) = (pos, vel);
    let rou_norm = rou_ref.norm();
    let rou_hat = rou_ref / rou_norm;
    let h = rou_ref.cross(v_ref);
    let h_hat = h.normalize();
    let theta_hat = h_hat.cross(&rou_hat);
    let a = Matrix3::from_columns(&[rou_hat, theta_hat, h_hat]).transpose();
    let [x, y, z] = to_array(h / (rou_norm * rou_norm)); // omega
    let zero = T::zero();
    let omega_w = Matrix3::from_row_slice(&[zero, -z, y, z, zero, -x, -y, x, zero]);
    let rou_rel = pos - rou_ref;
    let lvlh_vel = a * (j2000_vel - v_ref - (-a) * omega_w * rou_rel);
    LvlhVelFromJ2000ImplResult {
        lvlh_vel,
        a,
        rou_rel,
    }
}
