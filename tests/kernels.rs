use num_bigint::BigInt;
use ratmath::dyadic::{self, E, Iv};
use ratmath::{matrix, rat, wire};

fn p(n: i64, d: i64) -> (BigInt, BigInt) {
    rat::reduce(n.into(), d.into()).unwrap()
}

#[test]
fn determinant_small_and_modular_agree() {
    let m = vec![
        vec![p(2, 1), p(1, 3), p(0, 1), p(5, 1)],
        vec![p(-1, 2), p(4, 1), p(7, 1), p(1, 1)],
        vec![p(3, 1), p(0, 1), p(1, 5), p(-2, 1)],
        vec![p(1, 1), p(1, 1), p(1, 1), p(1, 1)],
    ];
    let a = matrix::determinant_bareiss(&m);
    let b = matrix::determinant_modular(&m);
    assert_eq!(a, b);
    assert_eq!(
        matrix::determinant(&[vec![p(1, 2), p(1, 3)], vec![p(1, 4), p(1, 5)]]),
        p(1, 60)
    );
}

#[test]
fn hilbert_inverse_is_integral() {
    let n = 6;
    let h: Vec<Vec<_>> = (0..n)
        .map(|i| (0..n).map(|j| p(1, (i + j + 1) as i64)).collect())
        .collect();
    let inv = matrix::inverse_pairs(&h).unwrap();
    assert!(inv.iter().flatten().all(|(_, d)| *d == BigInt::from(1)));
    assert_eq!(inv[0][0], p(36, 1));
}

#[test]
fn singular_and_solve_statuses() {
    let a = vec![vec![p(1, 1), p(2, 1)], vec![p(2, 1), p(4, 1)]];
    assert!(matrix::inverse_pairs(&a).is_none());
    assert_eq!(matrix::determinant(&a), p(0, 1));
    let s = matrix::solve_pairs(&a, &[vec![p(1, 1)], vec![p(2, 1)]]);
    assert_eq!(s.status, matrix::SolveStatus::Infinite);
    let s = matrix::solve_pairs(&a, &[vec![p(1, 1)], vec![p(3, 1)]]);
    assert_eq!(s.status, matrix::SolveStatus::None);
}

#[test]
fn pi_and_exp_enclose_the_truth() {
    let ctx = dyadic::context(60);
    let pi = dyadic::pi_interval(&ctx);
    let f = |e: &E| rat::to_f64(e.fin().unwrap()) / 2f64.powi(60);
    assert!(f(&pi.lo) <= std::f64::consts::PI && std::f64::consts::PI <= f(&pi.hi));
    let one = Iv {
        lo: E::Fin(ctx.one.clone()),
        hi: E::Fin(ctx.one.clone()),
    };
    let e = dyadic::exp(&ctx, &Some(one)).unwrap();
    assert!(f(&e.lo) <= std::f64::consts::E && std::f64::consts::E <= f(&e.hi));
}

#[test]
fn to_f64_rounds_like_js() {
    // 2^53 + 1 rounds to even (2^53); 2^53 + 3 rounds up to 2^53 + 4
    let two53 = BigInt::from(1u64 << 53);
    assert_eq!(rat::to_f64(&(&two53 + 1u8)), 9007199254740992.0);
    assert_eq!(rat::to_f64(&(&two53 + 3u8)), 9007199254740996.0);
    assert_eq!(rat::to_f64(&(-(&two53 + 3u8))), -9007199254740996.0);
    let big: BigInt = BigInt::from(3) << 2000u32;
    assert_eq!(rat::to_f64(&big), f64::INFINITY);
}

#[test]
fn wire_det_round_trip() {
    // 2x2 [[1/2, 1/3], [1/4, 1/5]] → det = 1/10 - 1/12 = 1/60
    let out = wire::call(wire::ops::DET, "2 2 1 2 1 3 1 4 1 5").unwrap();
    assert_eq!(out, "1 3c");
    assert!(wire::call(wire::ops::DET, "2 2 1 0 1 1 1 1 1 1").is_err());
    assert!(wire::call(99, "").is_err());
}
