use crate::*;

#[test]
fn vec_to_slice() {
    let mut mut_vec = Vector3::new(1.0, 2.0, 3.0);
    let vec = mut_vec.clone();
    let mut_slice = mut_vec.as_mut();
    let slice: &[f64; 3] = vec.as_ref();
    assert_eq!(mut_slice, slice);
}

#[test]
fn sat_state_to_vec6() {
    let mut mut_sat_state: PosVelState<_> = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0].into();
    let sat_state = mut_sat_state.clone();
    let mut_vec6 = mut_sat_state.as_mut();
    let vec6: &Vector6<_> = sat_state.as_ref();
    assert_eq!(vec6, mut_vec6);
}

#[test]
fn sat_state_to_slice() {
    let mut mut_sat_state =
        PosVelState::new(Vector3::new(1.0, 2.0, 3.0), Vector3::new(4.0, 5.0, 6.0));
    let sat_state = mut_sat_state.clone();
    let mut_slice: &mut [f64; 6] = mut_sat_state.as_mut();
    let slice: &[f64; 6] = sat_state.as_ref();
    assert_eq!(slice, mut_slice);
}

#[test]
fn vector3_lvlh_and_vvlh_trans() {
    let vvlh = Vector3::new(1.0, 2.0, 3.0);
    let lvlh = vector3::lvlh_from_vvlh(vvlh);
    assert_eq!(lvlh, Vector3::new(-3.0, 1.0, -2.0));
    let vvlh_back = vector3::vvlh_from_lvlh(lvlh);
    assert_eq!(vvlh_back, vvlh);
}

#[test]
fn from_lvlh_to_j2000() {
    let ref_state = [
        -5042137.452,
        42166300.0,
        14286.07219,
        -2732.96475,
        -456.767294,
        7.459819,
    ]
    .into();
    for (lvlh, j2000) in [
        (
            [-100000.0, 200000.0, -300000.0, 0.0, -5.0, 10.0],
            [
                -5229668.299157037,
                42043264.60687324,
                -285204.28391968715,
                -2719.9534663356803,
                -468.343192497209,
                17.424173822617966,
            ],
        ),
        (
            [500000.0, -200000.0, 400000.0, 0.0, 15.0, -5.0],
            [
                -4901826.073964696,
                42686504.709447816,
                413910.6178101943,
                -2781.7781127925523,
                -449.47424296750046,
                2.5932278090834116,
            ],
        ),
    ] {
        assert_eq!(
            J2000::try_from_lvlh(lvlh.into(), &ref_state).unwrap(),
            j2000.into()
        );
    }
}

#[test]
fn from_lvlh_vel_to_j2000() {
    let ref_state = [
        -5042137.452,
        42166300.0,
        14286.07219,
        -2732.96475,
        -456.767294,
        7.459819,
    ]
    .into();
    for (lvlh, j2000) in [
        (
            [0.0, -5.0, 10.0],
            [4.991929089449267, 0.5935385856962423, 9.986408559296978],
        ),
        (
            [0.0, 15.0, -5.0],
            [
                -14.907497440247987,
                -1.7809198808180138,
                -4.9593189499271855,
            ],
        ),
    ] {
        assert_eq!(
            vector3::try_j2000_vel_from_lvlh(lvlh.into(), &ref_state).unwrap(),
            Vector3::from(j2000)
        );
    }
}
