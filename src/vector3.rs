use crate::{J2000, Matrix3, MatrixError, Vector3};
use crate::{RealField, Scalar};

macro_rules! r { () => { Result<Vector3<T>, MatrixError> } }

pub fn norm_not_zero<T: RealField>(vec: &Vector3<T>) -> Result<T, MatrixError> {
    match vec.norm() {
        zero if zero.is_zero() => Err(MatrixError::NormalizeZero),
        norm => Ok(norm),
    }
}

pub fn normalize<T: RealField>(vec: &Vector3<T>) -> r!() {
    Ok(vec.unscale(norm_not_zero(vec)?))
}

pub fn ref_from_array<T>(arr: &[T; 3]) -> &Vector3<T> {
    unsafe { &*(arr as *const [T; 3] as *const Vector3<T>) }
}

pub fn to_array<T>(vec: Vector3<T>) -> [T; 3] {
    let [arr] = vec.data.0;
    arr
}

pub fn lvlh_from_vvlh<T: Scalar + std::ops::Neg<Output = T>>(vvlh: &Vector3<T>) -> Vector3<T> {
    Vector3::new(-vvlh.z.clone(), vvlh.x.clone(), -vvlh.y.clone())
}

pub fn vvlh_from_lvlh<T: Scalar + std::ops::Neg<Output = T>>(lvlh: &Vector3<T>) -> Vector3<T> {
    Vector3::new(lvlh.y.clone(), -lvlh.z.clone(), -lvlh.x.clone())
}

pub fn try_j2000_vel_from_lvlh<T: RealField>(lvlh_vel: &Vector3<T>, ref_state: &J2000<T>) -> r!() {
    let rou_hat = normalize(&ref_state.pos)?;
    let h_hat = normalize(&ref_state.pos.cross(&ref_state.vel))?;
    let theta_h = h_hat.cross(&rou_hat);
    let a = Matrix3::from_columns(&[rou_hat, theta_h, h_hat]).transpose();
    Ok(a.try_inverse().ok_or(MatrixError::InverseSingular)? * lvlh_vel)
}

pub fn try_j2000_vel_from_vvlh<T: RealField>(vvlh_vel: &Vector3<T>, ref_state: &J2000<T>) -> r!() {
    try_j2000_vel_from_lvlh(&lvlh_from_vvlh(vvlh_vel), ref_state)
}

pub fn try_lvlh_pos_from_j2000<T: RealField>(j2000_pos: &Vector3<T>, ref_state: &J2000<T>) -> r!() {
    let rou_hat = normalize(&ref_state.pos)?;
    let h_hat = normalize(&ref_state.pos.cross(&ref_state.vel))?;
    let theta_hat = h_hat.cross(&rou_hat);
    let a = Matrix3::from_columns(&[rou_hat, theta_hat, h_hat]).transpose();
    Ok(a * (j2000_pos - &ref_state.pos))
}

pub fn try_vvlh_pos_from_j2000<T: RealField>(j2000_pos: &Vector3<T>, ref_state: &J2000<T>) -> r!() {
    let lvlh_pos = try_lvlh_pos_from_j2000(j2000_pos, ref_state)?;
    Ok(vvlh_from_lvlh(&lvlh_pos))
}
