//! Exact matrix kernels over the rationals. Each routine is a faithful port of
//! the matching function in ratmath/linalg `src/exact.js` (same pivot choice,
//! same row order), so the outputs are identical, not merely equal in value.
//!
//! Matrices are `Vec<Vec<Pair>>` (rows of reduced rationals).

use crate::rat::{Pair, lcm, reduce_nz};
use num_bigint::BigInt;
use num_traits::{One, ToPrimitive, Zero};

pub type PairMatrix = Vec<Vec<Pair>>;

/// Clear denominators row by row. Returns integer rows and each row's scale
/// `L_i` (row i was multiplied by `L_i`).
pub fn integer_rows(pairs: &[Vec<Pair>]) -> (Vec<Vec<BigInt>>, Vec<BigInt>) {
    let mut scales = Vec::with_capacity(pairs.len());
    let rows = pairs
        .iter()
        .map(|row| {
            let mut l = BigInt::one();
            for (_, d) in row {
                if !d.is_one() {
                    l = lcm(&l, d);
                }
            }
            let out = row
                .iter()
                .map(|(n, d)| if &l == d { n.clone() } else { n * (&l / d) })
                .collect();
            scales.push(l);
            out
        })
        .collect();
    (rows, scales)
}

pub struct Bareiss {
    pub rows: Vec<Vec<BigInt>>,
    pub pivots: Vec<usize>,
    pub last: BigInt,
    pub sign: i32,
}

/// Fraction-free Gauss–Jordan (Bareiss), pivoting on the first `columns`
/// columns with the smallest nonzero pivot (as linalg does). After the run every
/// pivot entry equals `last`, so the RREF is `rows[k] / last`.
pub fn bareiss(rows: Vec<Vec<BigInt>>, columns: usize) -> Bareiss {
    if let Some(small) = bareiss_i128(&rows, columns) {
        return small;
    }
    bareiss_big(rows, columns)
}

/// The same elimination in machine integers. Bareiss entries are minors of
/// the input, so for the small matrices a calculator sees they fit in i128;
/// on the first overflow this gives up and the BigInt version runs instead.
fn bareiss_i128(rows: &[Vec<BigInt>], columns: usize) -> Option<Bareiss> {
    let m = rows.len();
    let width = rows.first().map_or(0, |r| r.len());
    let mut a: Vec<i128> = Vec::with_capacity(m * width);
    for row in rows {
        for v in row {
            a.push(v.to_i64()? as i128);
        }
    }
    let mut pivots = Vec::new();
    let mut prev: i128 = 1;
    let mut sign = 1;
    let mut r = 0;
    let mut col = 0;
    while col < columns && r < m {
        let mut found: Option<usize> = None;
        for i in r..m {
            let v = a[i * width + col];
            if v != 0 && found.is_none_or(|f| v.unsigned_abs() < a[f * width + col].unsigned_abs())
            {
                found = Some(i);
            }
        }
        let Some(found) = found else {
            col += 1;
            continue;
        };
        if found != r {
            for j in 0..width {
                a.swap(r * width + j, found * width + j);
            }
            sign = -sign;
        }
        let pbase = r * width;
        let p = a[pbase + col];
        for i in 0..m {
            if i == r {
                continue;
            }
            let base = i * width;
            let f = a[base + col];
            for j in 0..width {
                if j == col {
                    continue;
                }
                let x = a[base + j];
                let b = a[pbase + j];
                if x == 0 && (f == 0 || b == 0) {
                    continue;
                }
                let t = p.checked_mul(x)?.checked_sub(f.checked_mul(b)?)?;
                a[base + j] = t / prev;
            }
            a[base + col] = 0;
        }
        prev = p;
        pivots.push(col);
        r += 1;
        col += 1;
    }
    let rows = a
        .chunks(width.max(1))
        .take(m)
        .map(|row| row.iter().map(|&v| BigInt::from(v)).collect())
        .collect();
    let rows = if width == 0 {
        vec![Vec::new(); m]
    } else {
        rows
    };
    Some(Bareiss {
        rows,
        pivots,
        last: BigInt::from(prev),
        sign,
    })
}

fn bareiss_big(mut rows: Vec<Vec<BigInt>>, columns: usize) -> Bareiss {
    let m = rows.len();
    let width = rows.first().map_or(0, |r| r.len());
    let mut pivots = Vec::new();
    let mut prev = BigInt::one();
    let mut sign = 1;
    let mut r = 0;
    let mut col = 0;
    let mut t = BigInt::zero();
    while col < columns && r < m {
        let mut found: Option<usize> = None;
        for i in r..m {
            let v = &rows[i][col];
            if !v.is_zero() {
                match found {
                    None => found = Some(i),
                    Some(f) => {
                        if v.magnitude() < rows[f][col].magnitude() {
                            found = Some(i)
                        }
                    }
                }
            }
        }
        let Some(found) = found else {
            col += 1;
            continue;
        };
        if found != r {
            rows.swap(r, found);
            sign = -sign;
        }
        let pivot_row = std::mem::take(&mut rows[r]);
        let p = pivot_row[col].clone();
        for (i, row) in rows.iter_mut().enumerate() {
            if i == r {
                continue;
            }
            let f = std::mem::take(&mut row[col]);
            let f_zero = f.is_zero();
            for j in 0..width {
                if j == col {
                    continue;
                }
                let b = &pivot_row[j];
                if row[j].is_zero() && (f_zero || b.is_zero()) {
                    continue;
                }
                // row[j] = (p * a - f * b) / prev, exact
                t.clone_from(&row[j]);
                t *= &p;
                if !f_zero && !b.is_zero() {
                    t -= &f * b;
                }
                if !prev.is_one() {
                    t /= &prev;
                }
                std::mem::swap(&mut row[j], &mut t);
            }
            row[col] = BigInt::zero();
        }
        rows[r] = pivot_row;
        prev = p;
        pivots.push(col);
        r += 1;
        col += 1;
    }
    Bareiss {
        rows,
        pivots,
        last: prev,
        sign,
    }
}

/// Reduced row echelon form as rows of pairs plus the pivot columns.
pub fn rref_pairs(pairs: &[Vec<Pair>], columns: usize) -> (PairMatrix, Vec<usize>) {
    let (ints, _) = integer_rows(pairs);
    let Bareiss {
        rows, pivots, last, ..
    } = bareiss(ints, columns);
    let out = rows
        .into_iter()
        .map(|row| {
            row.into_iter()
                .map(|v| reduce_nz(v, last.clone()))
                .collect()
        })
        .collect();
    (out, pivots)
}

pub fn rank(pairs: &[Vec<Pair>], columns: usize) -> usize {
    let (ints, _) = integer_rows(pairs);
    bareiss(ints, columns).pivots.len()
}

pub fn determinant_bareiss(pairs: &[Vec<Pair>]) -> Pair {
    let n = pairs.len();
    if n == 0 {
        return (BigInt::one(), BigInt::one());
    }
    let (ints, scales) = integer_rows(pairs);
    let Bareiss {
        pivots, last, sign, ..
    } = bareiss(ints, n);
    if pivots.len() < n {
        return (BigInt::zero(), BigInt::one());
    }
    let denom: BigInt = scales.iter().product();
    reduce_nz(if sign < 0 { -last } else { last }, denom)
}

#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub enum SolveStatus {
    Unique = 0,
    Infinite = 1,
    None = 2,
}

pub struct Solution {
    pub status: SolveStatus,
    /// One solution (n rows × k columns; free unknowns set to 0), or None.
    pub x: Option<PairMatrix>,
    /// Null-space vectors of A, each of length n.
    pub basis: Vec<Vec<Pair>>,
}

fn zero_pair() -> Pair {
    (BigInt::zero(), BigInt::one())
}

/// Solve A x = B (B may have several columns), like linalg's `solvePairs`.
pub fn solve_pairs(a: &[Vec<Pair>], b: &[Vec<Pair>]) -> Solution {
    let n = a.first().map_or(0, |r| r.len());
    let k = b.first().map_or(0, |r| r.len());
    let wide: Vec<Vec<Pair>> = a
        .iter()
        .zip(b)
        .map(|(ra, rb)| ra.iter().chain(rb.iter()).cloned().collect())
        .collect();
    let (rows, pivots) = rref_pairs(&wide, n);
    for row in rows.iter().skip(pivots.len()) {
        if row[n..].iter().any(|(v, _)| !v.is_zero()) {
            return Solution {
                status: SolveStatus::None,
                x: None,
                basis: null_basis(&rows, &pivots, n),
            };
        }
    }
    let mut x = vec![vec![zero_pair(); k]; n];
    for (r, &col) in pivots.iter().enumerate() {
        for j in 0..k {
            x[col][j] = rows[r][n + j].clone();
        }
    }
    let basis = null_basis(&rows, &pivots, n);
    let status = if basis.is_empty() {
        SolveStatus::Unique
    } else {
        SolveStatus::Infinite
    };
    Solution {
        status,
        x: Some(x),
        basis,
    }
}

pub fn null_basis(rows: &[Vec<Pair>], pivots: &[usize], n: usize) -> Vec<Vec<Pair>> {
    let mut is_pivot = vec![false; n];
    for &p in pivots {
        is_pivot[p] = true;
    }
    let mut basis = Vec::new();
    for f in 0..n {
        if is_pivot[f] {
            continue;
        }
        let mut v = vec![zero_pair(); n];
        v[f] = (BigInt::one(), BigInt::one());
        for (r, &col) in pivots.iter().enumerate() {
            let (a, b) = &rows[r][f];
            v[col] = (-a, b.clone());
        }
        basis.push(v);
    }
    basis
}

pub fn inverse_pairs(pairs: &[Vec<Pair>]) -> Option<PairMatrix> {
    let n = pairs.len();
    let wide: Vec<Vec<Pair>> = pairs
        .iter()
        .enumerate()
        .map(|(i, row)| {
            let mut w = row.clone();
            w.extend((0..n).map(|j| (BigInt::from((i == j) as u8), BigInt::one())));
            w
        })
        .collect();
    let (rows, pivots) = rref_pairs(&wide, n);
    if pivots.len() < n {
        return None;
    }
    Some(rows.into_iter().map(|row| row[n..].to_vec()).collect())
}

// ---------------------------------------------------------------------------
// Multi-modular determinant. Same idea as linalg's `determinantModular`
// (determinant mod many word-size primes, Chinese remaindering past twice the
// Hadamard bound) but with 31-bit primes: all the O(n³) work is one u64
// multiply-add and remainder per entry (a·b + c < 2^63), native in x86, ARM
// and WASM alike. The answer is exact, so it equals the JS result.

const PRIME_TOP: u32 = 1 << 31;

fn is_prime_u32(p: u32) -> bool {
    if p < 2 {
        return false;
    }
    for &q in &[2u32, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37] {
        if p % q == 0 {
            return p == q;
        }
    }
    // deterministic Miller–Rabin for 32-bit n with bases 2, 7, 61
    let n = p as u64;
    let mut d = n - 1;
    let mut s = 0;
    while d % 2 == 0 {
        d /= 2;
        s += 1;
    }
    'outer: for &a in &[2u64, 7, 61] {
        let mut x = pow_mod(a % n, d, n);
        if x == 1 || x == n - 1 {
            continue;
        }
        for _ in 1..s {
            x = x * x % n;
            if x == n - 1 {
                continue 'outer;
            }
        }
        return false;
    }
    true
}

fn pow_mod(mut b: u64, mut e: u64, m: u64) -> u64 {
    let mut r = 1u64;
    while e > 0 {
        if e & 1 == 1 {
            r = r * b % m;
        }
        b = b * b % m;
        e >>= 1;
    }
    r
}

thread_local! {
    static PRIMES: std::cell::RefCell<Vec<u32>> = const { std::cell::RefCell::new(Vec::new()) };
}

/// The i-th prime below 2^31, counting down.
pub fn prime(i: usize) -> u32 {
    PRIMES.with(|ps| {
        let mut ps = ps.borrow_mut();
        let mut p = ps.last().map_or(PRIME_TOP - 1, |&l| l - 2);
        while ps.len() <= i {
            if is_prime_u32(p) {
                ps.push(p);
            }
            p -= 2;
        }
        ps[i]
    })
}

fn inv_mod(a: u64, p: u64) -> u64 {
    pow_mod(a % p, p - 2, p)
}

fn residue(v: &BigInt, p: u32) -> u64 {
    if let Some(s) = v.to_i64() {
        return s.rem_euclid(p as i64) as u64;
    }
    let r = (v % p).to_i64().unwrap();
    (if r < 0 { r + p as i64 } else { r }) as u64
}

/// x mod p for x < 2^63 without a hardware divide (Barrett reduction: the
/// estimated quotient is short by at most one, which the fix-up corrects). A
/// `%` by a runtime divisor is a slow 64-bit divide on x86 and in WASM.
#[derive(Clone, Copy)]
struct Reducer {
    p: u64,
    #[cfg(not(target_arch = "wasm32"))]
    m: u64, // floor(2^64 / p), Barrett constant
    #[cfg(target_arch = "wasm32")]
    pinv: f64,
}

impl Reducer {
    fn new(p: u32) -> Self {
        Reducer {
            p: p as u64,
            #[cfg(not(target_arch = "wasm32"))]
            m: u64::MAX / p as u64,
            #[cfg(target_arch = "wasm32")]
            pinv: 1.0 / p as f64,
        }
    }
    /// Native: Barrett with a 64×64→128 multiply (one instruction on x86-64 and arm64).
    #[cfg(not(target_arch = "wasm32"))]
    #[inline(always)]
    fn reduce(&self, x: u64) -> u64 {
        let q = ((x as u128 * self.m as u128) >> 64) as u64;
        let r = x - q * self.p;
        if r >= self.p { r - self.p } else { r }
    }
    /// WASM has no wide multiply; a double-precision quotient is off by at most one.
    #[cfg(target_arch = "wasm32")]
    #[inline(always)]
    fn reduce(&self, x: u64) -> u64 {
        let p = self.p as i64;
        let q = (x as f64 * self.pinv) as i64;
        let mut r = x as i64 - q * p;
        if r < 0 {
            r += p;
        } else if r >= p {
            r -= p;
        }
        r as u64
    }
}

/// Determinant mod p of an integer matrix.
pub fn det_mod(rows: &[Vec<BigInt>], p: u32) -> u64 {
    let n = rows.len();
    let pp = p as u64;
    let red = Reducer::new(p);
    let mut a = vec![0u64; n * n];
    for i in 0..n {
        for j in 0..n {
            a[i * n + j] = residue(&rows[i][j], p);
        }
    }
    let mut det = 1u64;
    let mut neg = false;
    for col in 0..n {
        let Some(found) = (col..n).find(|&i| a[i * n + col] != 0) else {
            return 0;
        };
        if found != col {
            for j in col..n {
                a.swap(col * n + j, found * n + j);
            }
            neg = !neg;
        }
        let base = col * n;
        let pivot = a[base + col];
        det = det * pivot % pp;
        let inv = inv_mod(pivot, pp);
        let (head, tail) = a.split_at_mut((col + 1) * n);
        let prow = &head[base + col + 1..base + n];
        for ri in tail.chunks_exact_mut(n) {
            let f = ri[col] * inv % pp;
            if f == 0 {
                continue;
            }
            let g = pp - f; // add g·pivotRow instead of subtracting f·pivotRow
            for (x, &b) in ri[col + 1..].iter_mut().zip(prow) {
                *x = red.reduce(*x + g * b);
            }
        }
    }
    if neg && det != 0 { pp - det } else { det }
}

/// log2 of the Hadamard bound Π‖row‖₂, rounded up per row (0 if a row is zero).
pub fn hadamard_bits(rows: &[Vec<BigInt>]) -> u64 {
    let mut bits = 0;
    for row in rows {
        let s: BigInt = row.iter().map(|v| v * v).sum();
        if s.is_zero() {
            return 0;
        }
        bits += s.bits().div_ceil(2);
    }
    bits
}

pub fn determinant_modular(pairs: &[Vec<Pair>]) -> Pair {
    let n = pairs.len();
    if n == 0 {
        return (BigInt::one(), BigInt::one());
    }
    let (rows, scales) = integer_rows(pairs);
    let need = (hadamard_bits(&rows) + 2) as f64;
    let mut x = BigInt::zero();
    let mut m = BigInt::one();
    let mut bits = 0f64;
    let mut i = 0;
    while bits < need {
        let p = prime(i);
        let r = det_mod(&rows, p);
        let pp = p as u64;
        let xm = residue(&x, p);
        let t = ((r + pp - xm) % pp) * inv_mod(residue(&m, p), pp) % pp;
        x += &m * t;
        m *= p;
        bits += (p as f64).log2();
        i += 1;
    }
    if x > (&m >> 1) {
        x -= &m;
    }
    let denom: BigInt = scales.iter().product();
    reduce_nz(x, denom)
}

/// Determinant: Bareiss for n ≤ 3 (as linalg does), multi-modular otherwise.
pub fn determinant(pairs: &[Vec<Pair>]) -> Pair {
    if pairs.len() <= 3 {
        determinant_bareiss(pairs)
    } else {
        determinant_modular(pairs)
    }
}

/// True when A is certainly nonsingular (det ≠ 0 mod one of `tries` primes).
pub fn provably_nonsingular(pairs: &[Vec<Pair>], tries: usize) -> bool {
    let (rows, _) = integer_rows(pairs);
    (0..tries).any(|i| det_mod(&rows, prime(i)) != 0)
}
