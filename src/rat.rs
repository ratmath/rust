//! Exact rationals as `(numerator, denominator)` pairs of `BigInt`, the
//! denominator positive and the pair reduced. Bare pairs (not a class) so the
//! matrix kernels don't renormalize on every step; this mirrors `src/q.js` in
//! ratmath/linalg.

use num_bigint::{BigInt, Sign};
use num_integer::Integer;
use num_traits::{One, Signed, Zero};

/// A reduced rational `n / d` with `d > 0`.
pub type Pair = (BigInt, BigInt);

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ZeroDenominator;

/// Reduce `n / d`, moving the sign to the numerator.
pub fn reduce(mut n: BigInt, mut d: BigInt) -> Result<Pair, ZeroDenominator> {
    if d.is_zero() {
        return Err(ZeroDenominator);
    }
    if d.is_negative() {
        n = -n;
        d = -d;
    }
    let g = n.gcd(&d);
    if g > BigInt::one() {
        n /= &g;
        d /= &g;
    }
    Ok((n, d))
}

/// Reduce when the denominator is known to be nonzero.
pub fn reduce_nz(n: BigInt, d: BigInt) -> Pair {
    reduce(n, d).expect("nonzero denominator")
}

pub fn lcm(a: &BigInt, b: &BigInt) -> BigInt {
    if a.is_zero() || b.is_zero() {
        return BigInt::zero();
    }
    (a / a.gcd(b) * b).abs()
}

/// Number of bits in `|n|` (0 for 0), like `n.toString(2).length` in JS.
pub fn bit_length(n: &BigInt) -> u64 {
    n.bits()
}

/// Nearest double to a BigInt, ties to even — the same rounding as JS `Number(bigint)`.
pub fn to_f64(n: &BigInt) -> f64 {
    let bits = n.bits();
    let neg = n.sign() == Sign::Minus;
    let mag = n.magnitude();
    let v = if bits <= 53 {
        let mut x = 0u64;
        for (i, d) in mag.iter_u64_digits().enumerate() {
            x |= d << (64 * i);
        }
        x as f64
    } else {
        let shift = bits - 54; // keep 54 bits: 53 + one rounding bit
        let q: num_bigint::BigUint = mag >> shift;
        let q = q.iter_u64_digits().next().unwrap_or(0);
        let sticky = mag.trailing_zeros().map_or(false, |tz| tz < shift);
        let round_bit = q & 1;
        let mut m = q >> 1;
        if round_bit == 1 && (sticky || m & 1 == 1) {
            m += 1;
        }
        let e = (shift + 1) as i32;
        if e > 1100 {
            f64::INFINITY
        } else {
            (m as f64) * 2f64.powi(e)
        }
    };
    if neg { -v } else { v }
}
