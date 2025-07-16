use isomorphic::{impl_mut, impl_ref, impl_transmute};

use crate::{Matrix3, MatrixError, RealField, Scalar, Vector3, Vector6, vector3};

#[cfg(feature = "serde")]
use serde::{Deserialize, Serialize};

/// Satellite state with position in `m` and velocity in `m/s` .
#[repr(C)]
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(Serialize, Deserialize))]
pub struct PosVel<T: Scalar = f64> {
    /// Position in `m`
    pub pos: Vector3<T>,
    /// Velocity in `m/s`
    pub vel: Vector3<T>,
}

impl<T: Scalar + Default> Default for PosVel<T> {
    fn default() -> Self {
        Self {
            pos: Vector3::default(),
            vel: Vector3::default(),
        }
    }
}

impl<T: Scalar> PosVel<T> {
    pub const fn new(pos: Vector3<T>, vel: Vector3<T>) -> Self {
        Self { pos, vel }
    }
}

impl_ref!(PosVel<T> = [T; 6] where T: Scalar);
impl_ref!(PosVel<T> = Vector6<T> where T: Scalar);
impl_mut!(PosVel<T> = [T; 6] where T: Scalar);
impl_mut!(PosVel<T> = Vector6<T> where T: Scalar);
impl_transmute!(PosVel<T> = [T; 6] where T: Scalar);
impl_transmute!(PosVel<T> = Vector6<T> where T: Scalar);

macro_rules! specific {
    ($base: ident => $name:ident) => {
        #[repr(transparent)]
        #[derive(Debug, Clone, PartialEq)]
        #[cfg_attr(feature = "serde", derive(Serialize, Deserialize), serde(transparent))]
        pub struct $name<T: Scalar = f64>(pub $base<T>);

        impl<T: Scalar> From<$base<T>> for $name<T> {
            fn from(state: $base<T>) -> Self {
                Self(state)
            }
        }

        impl<T: Scalar> From<$name<T>> for $base<T> {
            fn from(state: $name<T>) -> Self {
                state.0
            }
        }

        impl<T: Scalar> std::ops::Deref for $name<T> {
            type Target = $base<T>;
            fn deref(&self) -> &Self::Target {
                &self.0
            }
        }

        impl<T: Scalar> std::ops::DerefMut for $name<T> {
            fn deref_mut(&mut self) -> &mut Self::Target {
                &mut self.0
            }
        }

        impl<T: Scalar> $name<T> {
            pub const fn new(pos: Vector3<T>, vel: Vector3<T>) -> Self {
                Self($base::new(pos, vel))
            }
        }

        impl<T: Scalar + Default> Default for $name<T> {
            fn default() -> Self {
                Self($base::default())
            }
        }

        impl_ref!($name<T> = Vector6<T> where T: Scalar);
        impl_ref!($name<T> = [T; 6] where T: Scalar);
        impl_ref!($name<T> = $base<T> where T: Scalar);
        impl_mut!($name<T> = Vector6<T> where T: Scalar);
        impl_mut!($name<T> = [T; 6] where T: Scalar);
        impl_mut!($name<T> = $base<T> where T: Scalar);
        impl_transmute!($name<T> = Vector6<T> where T: Scalar);
        impl_transmute!($name<T> = [T; 6] where T: Scalar);
    };
    ($base: ident => $($name:ident)|+) => { $(specific!($base => $name);)* };
}

specific!(PosVel => J2000 | LVLH | VVLH);

impl<T: RealField + Copy> J2000<T> {
    pub fn try_from_lvlh(lvlh: LVLH<T>, ref_state: &Self) -> Result<Self, MatrixError> {
        let rou_norm = vector3::norm_not_zero(&ref_state.pos)?;
        let rou_hat = ref_state.pos / rou_norm;
        let h = ref_state.pos.cross(&ref_state.vel);
        let h_hat = vector3::normalize(&h)?;
        let a = Matrix3::from_columns(&[rou_hat, h_hat.cross(&rou_hat), h_hat]).transpose();
        let inv_a = a.try_inverse().ok_or(MatrixError::InverseSingular)?;
        let lvlh_pos_j2000 = inv_a * lvlh.pos;
        let [x, y, z] = vector3::to_array(h / (rou_norm * rou_norm));
        let zero = T::zero();
        let omega_w = Matrix3::from_row_slice(&[zero, -z, y, z, zero, -x, -y, x, zero]);
        let vel_j2000 = ref_state.vel + inv_a * (lvlh.vel + a * omega_w * lvlh_pos_j2000);
        Ok(J2000::new(ref_state.pos + lvlh_pos_j2000, vel_j2000))
    }

    pub fn try_from_vvlh(vvlh: VVLH<T>, ref_state: &Self) -> Result<Self, MatrixError> {
        J2000::try_from_lvlh(vvlh.into(), ref_state)
    }
}

impl<T: RealField + Copy> LVLH<T> {
    pub fn try_from_j2000(j2000: J2000<T>, ref_state: &J2000<T>) -> Result<Self, MatrixError> {
        let rou_norm = vector3::norm_not_zero(&ref_state.pos)?;
        let rou_hat = ref_state.pos / rou_norm;
        let h = &ref_state.pos.cross(&ref_state.vel);
        let h_hat = vector3::normalize(h)?;
        let a = Matrix3::from_columns(&[rou_hat, h_hat.cross(&rou_hat), h_hat]).transpose();
        let [x, y, z] = vector3::to_array(h / (rou_norm * rou_norm));
        let zero = T::zero();
        let omega_w = Matrix3::from_row_slice(&[zero, -z, y, z, zero, -x, -y, x, zero]);
        let rou_rel = j2000.pos - ref_state.pos;
        let lvlh_vel = a * (j2000.vel - ref_state.vel - omega_w * rou_rel);
        Ok(LVLH::new(a * rou_rel, lvlh_vel))
    }
}

impl<T: RealField + Copy> VVLH<T> {
    pub fn try_from_j2000(j2000: J2000<T>, ref_state: &J2000<T>) -> Result<Self, MatrixError> {
        Ok(LVLH::try_from_j2000(j2000, ref_state)?.into())
    }
}

impl<T: Scalar + std::ops::Neg<Output = T>> From<LVLH<T>> for VVLH<T> {
    fn from(lvlh: LVLH<T>) -> Self {
        use vector3::vvlh_from_lvlh as to_vvlh;
        VVLH::new(to_vvlh(&lvlh.pos), to_vvlh(&lvlh.vel))
    }
}

impl<T: Scalar + std::ops::Neg<Output = T>> From<VVLH<T>> for LVLH<T> {
    fn from(vvlh: VVLH<T>) -> Self {
        use vector3::lvlh_from_vvlh as to_lvlh;
        LVLH::new(to_lvlh(&vvlh.pos), to_lvlh(&vvlh.vel))
    }
}
