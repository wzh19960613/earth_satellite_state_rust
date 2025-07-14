macro_rules! isomorphic_ref {
    ($t_from: ty => $t_to: ty) => {
        impl<T: Scalar> std::convert::AsRef<$t_to> for $t_from {
            fn as_ref(&self) -> &$t_to {
                unsafe { &*(self as *const Self as *const $t_to) }
            }
        }

        impl<T: Scalar> std::convert::AsMut<$t_to> for $t_from {
            fn as_mut(&mut self) -> &mut $t_to {
                unsafe { &mut *(self as *mut Self as *mut $t_to) }
            }
        }
    };
}

macro_rules! isomorphic_transmute {
    ($t_from: ty => $t_to: ty) => {
        impl<T: Scalar> From<$t_from> for $t_to {
            fn from(t: $t_from) -> Self {
                let r: &$t_to = t.as_ref();
                r.clone()
            }
        }
    };
}

macro_rules! impls {
    ($macro: ident ! ($t_from: ty => $($t_to: ty)|+)) => { $($macro!($t_from => $t_to);)+ };
    ($macro: ident ! ($t_from: ty = $($t_to: ty)|+)) => { $(
        $macro!($t_from => $t_to);
        $macro!($t_to => $t_from);
    )+ };
}

pub(crate) use impls;
pub(crate) use isomorphic_ref;
pub(crate) use isomorphic_transmute;
