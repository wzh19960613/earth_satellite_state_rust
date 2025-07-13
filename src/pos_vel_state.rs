use crate::{InverseOfSingularMatrix, Matrix3, Vector3, Vector6, vector3};
use crate::{RealField, Scalar, impl_isomorphic};

/// Satellite state with position in `m` and velocity in `m/s` .
#[repr(C)]
#[derive(Debug, Clone, PartialEq)]
pub struct PosVelState<T> {
    /// Position in `m`
    pub pos: Vector3<T>,
    /// Velocity in `m/s`
    pub vel: Vector3<T>,
}

impl<T: Scalar + Default> Default for PosVelState<T> {
    fn default() -> Self {
        Self {
            pos: Vector3::default(),
            vel: Vector3::default(),
        }
    }
}

impl<T> PosVelState<T> {
    pub const fn new(pos: Vector3<T>, vel: Vector3<T>) -> Self {
        Self { pos, vel }
    }
}

impl_isomorphic!(PosVelState<T> = Vector6<T> | [T; 6]);

macro_rules! specific_pos_vel_state {
    ($name:ident) => {
        #[repr(transparent)]
        #[derive(Debug, Clone, PartialEq)]
        pub struct $name<T>(pub PosVelState<T>);

        impl<T> From<PosVelState<T>> for $name<T> {
            fn from(state: PosVelState<T>) -> Self {
                Self(state)
            }
        }

        impl<T: Clone> From<[T; 6]> for $name<T> {
            fn from(slice: [T; 6]) -> Self {
                let pos_vel: PosVelState<T> = slice.into();
                Self(pos_vel)
            }
        }

        impl<T> std::ops::Deref for $name<T> {
            type Target = PosVelState<T>;
            fn deref(&self) -> &Self::Target {
                &self.0
            }
        }

        impl<T> std::ops::DerefMut for $name<T> {
            fn deref_mut(&mut self) -> &mut Self::Target {
                &mut self.0
            }
        }

        impl<T> $name<T> {
            pub const fn new(pos: Vector3<T>, vel: Vector3<T>) -> Self {
                Self(PosVelState::new(pos, vel))
            }
        }

        impl<T: Scalar + Default> Default for $name<T> {
            fn default() -> Self {
                Self(PosVelState::default())
            }
        }
    };
    ($($name:ident),*) => { $(specific_pos_vel_state!($name);)* };
}

specific_pos_vel_state!(J2000, LVLH, VVLH);

impl<T: RealField + Copy> J2000<T> {
    pub fn try_from_lvlh(lvlh: LVLH<T>, ref_state: &Self) -> Result<Self, InverseOfSingularMatrix> {
        let Self(PosVelState { pos, vel }) = ref_state;
        let (rou_ref, v_ref) = (pos, vel);
        let LVLH(PosVelState { pos, vel }) = lvlh;
        let rou_norm = rou_ref.norm();
        let rou_hat = rou_ref / rou_norm;
        let h = rou_ref.cross(v_ref);
        let h_hat = h.normalize();
        let theta_hat = h_hat.cross(&rou_hat);
        let a = Matrix3::from_columns(&[rou_hat, theta_hat, h_hat]).transpose();
        let inv_a = a.try_inverse().ok_or(InverseOfSingularMatrix)?;
        let lvlh_pos_j2000 = inv_a * pos;
        let rou_j2000 = rou_ref + lvlh_pos_j2000;
        let [x, y, z] = vector3::to_array(h / (rou_norm * rou_norm)); // omega
        let zero = T::zero();
        let omega_w = Matrix3::from_row_slice(&[zero, -z, y, z, zero, -x, -y, x, zero]);
        let v_j2000 = v_ref + inv_a * (vel + (a * omega_w) * lvlh_pos_j2000);
        Ok(J2000::new(rou_j2000, v_j2000))
    }

    pub fn try_from_vvlh(vvlh: VVLH<T>, ref_state: &Self) -> Result<Self, InverseOfSingularMatrix> {
        J2000::try_from_lvlh(vvlh.into(), ref_state)
    }
}

impl<T: RealField + Copy> LVLH<T> {
    pub fn from_j2000(j2000: J2000<T>, ref_state: &J2000<T>) -> Self {
        let vals = vector3::lvlh_vel_from_j2000_impl(j2000.vel, ref_state);
        LVLH::new(vals.a * vals.rou_rel, vals.lvlh_vel)
    }
}

impl<T: RealField + Copy> VVLH<T> {
    pub fn from_j2000(j2000: J2000<T>, ref_state: &J2000<T>) -> Self {
        LVLH::from_j2000(j2000, ref_state).into()
    }
}

impl<T: Scalar + std::ops::Neg<Output = T>> From<LVLH<T>> for VVLH<T> {
    fn from(lvlh: LVLH<T>) -> Self {
        use vector3::vvlh_from_lvlh as to_vvlh;
        VVLH::new(to_vvlh(lvlh.0.pos), to_vvlh(lvlh.0.vel))
    }
}

impl<T: Scalar + std::ops::Neg<Output = T>> From<VVLH<T>> for LVLH<T> {
    fn from(vvlh: VVLH<T>) -> Self {
        use vector3::lvlh_from_vvlh as to_lvlh;
        LVLH::new(to_lvlh(vvlh.0.pos), to_lvlh(vvlh.0.vel))
    }
}
