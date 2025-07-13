macro_rules! impl_isomorphic {
    ($t_from: ty => $t_to: ty) => {
        impl<T> std::convert::AsRef<$t_to> for $t_from {
            fn as_ref(&self) -> &$t_to {
                unsafe { &*(self as *const Self as *const $t_to) }
            }
        }

        impl<T> std::convert::AsMut<$t_to> for $t_from {
            fn as_mut(&mut self) -> &mut $t_to {
                unsafe { &mut *(self as *mut Self as *mut $t_to) }
            }
        }

        impl<T: Clone> From<$t_from> for $t_to {
            fn from(t: $t_from) -> Self {
                let r: &$t_to = t.as_ref();
                r.clone()
            }
        }
    };

    ($($t_from: ty => $t_to: ty),*) => { $(impl_isomorphic!($t_from => $t_to);)* };

    ($t_from: ty = $t_to: ty) => { impl_isomorphic!($t_from => $t_to, $t_to => $t_from); };

    ($t_from: ty = $($t_to: ty)|*) => { $(impl_isomorphic!($t_from = $t_to);)* };
}

pub(crate) use impl_isomorphic;
