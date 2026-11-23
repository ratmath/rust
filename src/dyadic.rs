//! Certified interval arithmetic on dyadic fixed-point numbers: a faithful
//! port of ratmath/graphics `src/interval.js`. An endpoint is `Fin(m)` for the
//! exact rational m / 2^k, or ±∞. Lower endpoints round down, upper endpoints
//! round up, so the true result is always inside. Every algorithm (series,
//! guard bits, range reduction, rounding) follows the JS step for step so the
//! two produce bit-identical intervals.
//!
//! `None` is the empty interval (e.g. sqrt of a negative interval).

use crate::rat::{bit_length, to_f64};
use num_bigint::BigInt;
use num_integer::Integer;
use num_traits::{One, Pow, Signed, Zero};
use std::cell::OnceCell;
use std::collections::HashMap;
use std::rc::Rc;

#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub enum E {
    NegInf,
    Fin(BigInt),
    PosInf,
}
use E::{Fin, NegInf, PosInf};

impl E {
    pub fn is_fin(&self) -> bool {
        matches!(self, Fin(_))
    }
    fn is_zero(&self) -> bool {
        matches!(self, Fin(v) if v.is_zero())
    }
    fn sgn(&self) -> i32 {
        match self {
            NegInf => -1,
            PosInf => 1,
            Fin(v) => match v.sign() {
                num_bigint::Sign::Minus => -1,
                num_bigint::Sign::NoSign => 0,
                num_bigint::Sign::Plus => 1,
            },
        }
    }
    fn neg(&self) -> E {
        match self {
            NegInf => PosInf,
            PosInf => NegInf,
            Fin(v) => Fin(-v),
        }
    }
    pub fn fin(&self) -> Option<&BigInt> {
        match self {
            Fin(v) => Some(v),
            _ => None,
        }
    }
    fn gt0(&self) -> bool {
        match self {
            PosInf => true,
            NegInf => false,
            Fin(v) => v.is_positive(),
        }
    }
    fn lt0(&self) -> bool {
        match self {
            PosInf => false,
            NegInf => true,
            Fin(v) => v.is_negative(),
        }
    }
    fn ge0(&self) -> bool {
        !self.lt0()
    }
    fn le0(&self) -> bool {
        !self.gt0()
    }
}

fn zero() -> E {
    Fin(BigInt::zero())
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Iv {
    pub lo: E,
    pub hi: E,
}

pub type I = Option<Iv>;

/// Raised where the JS kernel throws a RangeError (callers treat it as "anything").
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Overflow;

pub fn whole() -> Iv {
    Iv {
        lo: NegInf,
        hi: PosInf,
    }
}

pub fn floor_div(a: &BigInt, b: &BigInt) -> BigInt {
    a.div_floor(b)
}

pub fn ceil_div(a: &BigInt, b: &BigInt) -> BigInt {
    -((-a).div_floor(b))
}

/// Arithmetic shift right (floor), as BigInt `>>` in JS.
fn shr(a: &BigInt, s: u64) -> BigInt {
    if a.is_negative() {
        -((-a - 1u8) >> s) - 1u8
    } else {
        a >> s
    }
}

fn shr_up(a: &BigInt, s: u64) -> BigInt {
    -shr(&-a, s)
}

pub fn isqrt(n: &BigInt) -> BigInt {
    n.sqrt()
}

/// A precision context: `k` fractional bits, `one` = 2^k. Cached per k, with
/// lazily computed constants (pi, ln 2), like `context(k)` in the JS.
pub struct Ctx {
    pub k: i64,
    pub one: BigInt,
    pi: OnceCell<(BigInt, BigInt)>,
    ln2: OnceCell<(BigInt, BigInt)>,
}

thread_local! {
    static CONTEXTS: std::cell::RefCell<HashMap<i64, Rc<Ctx>>> = std::cell::RefCell::new(HashMap::new());
}

pub fn context(k: i64) -> Rc<Ctx> {
    let k = k.clamp(8, 4000);
    CONTEXTS.with(|c| {
        c.borrow_mut()
            .entry(k)
            .or_insert_with(|| {
                Rc::new(Ctx {
                    k,
                    one: BigInt::one() << k,
                    pi: OnceCell::new(),
                    ln2: OnceCell::new(),
                })
            })
            .clone()
    })
}

impl Ctx {
    fn kk(&self) -> u64 {
        self.k as u64
    }
}

/// Outward-rounded enclosure of the exact rational n/d (d ≠ 0).
pub fn from_rational(ctx: &Ctx, n: &BigInt, d: &BigInt) -> Iv {
    let s = n << ctx.kk();
    Iv {
        lo: Fin(floor_div(&s, d)),
        hi: Fin(ceil_div(&s, d)),
    }
}

pub fn hull(a: I, b: I) -> I {
    match (a, b) {
        (None, b) => b,
        (a, None) => a,
        (Some(a), Some(b)) => Some(Iv {
            lo: a.lo.min(b.lo),
            hi: a.hi.max(b.hi),
        }),
    }
}

fn add_e(a: &E, b: &E) -> E {
    match (a, b) {
        (Fin(x), Fin(y)) => Fin(x + y),
        (a, _) if !a.is_fin() => a.clone(),
        (_, b) => b.clone(),
    }
}

fn mul_e(ctx: &Ctx, a: &E, b: &E, down: bool) -> E {
    if a.is_zero() || b.is_zero() {
        return zero();
    }
    match (a, b) {
        (Fin(x), Fin(y)) => {
            let p = x * y;
            Fin(if down {
                shr(&p, ctx.kk())
            } else {
                shr_up(&p, ctx.kk())
            })
        }
        _ => {
            if a.sgn() * b.sgn() > 0 {
                PosInf
            } else {
                NegInf
            }
        }
    }
}

pub fn add(a: &I, b: &I) -> I {
    let (a, b) = (a.as_ref()?, b.as_ref()?);
    Some(Iv {
        lo: add_e(&a.lo, &b.lo),
        hi: add_e(&a.hi, &b.hi),
    })
}

pub fn neg(a: &I) -> I {
    let a = a.as_ref()?;
    Some(Iv {
        lo: a.hi.neg(),
        hi: a.lo.neg(),
    })
}

pub fn sub(a: &I, b: &I) -> I {
    if a.is_none() {
        return None;
    }
    add(a, &neg(b))
}

pub fn mul(ctx: &Ctx, a: &I, b: &I) -> I {
    let (a, b) = (a.as_ref()?, b.as_ref()?);
    let mut lo = PosInf;
    let mut hi = NegInf;
    for x in [&a.lo, &a.hi] {
        for y in [&b.lo, &b.hi] {
            let l = mul_e(ctx, x, y, true);
            let h = mul_e(ctx, x, y, false);
            if l < lo {
                lo = l;
            }
            if h > hi {
                hi = h;
            }
        }
    }
    Some(Iv { lo, hi })
}

fn recip_pos(ctx: &Ctx, lo: &E, hi: &E) -> Iv {
    let one2 = &ctx.one * &ctx.one;
    Iv {
        lo: hi.fin().map_or(zero(), |h| Fin(floor_div(&one2, h))),
        hi: lo.fin().map_or(zero(), |l| Fin(ceil_div(&one2, l))),
    }
}

pub fn recip(ctx: &Ctx, a: &I) -> I {
    let a = a.as_ref()?;
    if a.lo.gt0() {
        return Some(recip_pos(ctx, &a.lo, &a.hi));
    }
    if a.hi.lt0() {
        let lo2 = if a.lo.is_fin() { a.lo.neg() } else { PosInf };
        return neg(&Some(recip_pos(ctx, &a.hi.neg(), &lo2)));
    }
    if a.lo.is_zero() && a.hi.is_zero() {
        return None;
    }
    let one2 = || &ctx.one * &ctx.one;
    if a.lo.is_zero() {
        return Some(Iv {
            lo: a.hi.fin().map_or(zero(), |h| Fin(floor_div(&one2(), h))),
            hi: PosInf,
        });
    }
    if a.hi.is_zero() {
        return Some(Iv {
            lo: NegInf,
            hi: a.lo.fin().map_or(zero(), |l| Fin(ceil_div(&one2(), l))),
        });
    }
    Some(whole())
}

pub fn div(ctx: &Ctx, a: &I, b: &I) -> I {
    if a.is_none() || b.is_none() {
        return None;
    }
    mul(ctx, a, &recip(ctx, b))
}

fn pow_e(ctx: &Ctx, v: &E, n: u64, down: bool) -> E {
    match v {
        Fin(x) => {
            let p: BigInt = Pow::pow(x, n);
            let shift = ctx.kk() * (n - 1);
            Fin(if down {
                shr(&p, shift)
            } else {
                shr_up(&p, shift)
            })
        }
        _ => {
            if n % 2 == 0 {
                PosInf
            } else {
                v.clone()
            }
        }
    }
}

pub fn pow_int(ctx: &Ctx, a: &I, n: i64) -> I {
    let iv = a.as_ref()?;
    if n == 0 {
        return Some(Iv {
            lo: Fin(ctx.one.clone()),
            hi: Fin(ctx.one.clone()),
        });
    }
    if n < 0 {
        return recip(ctx, &pow_int(ctx, a, -n));
    }
    if n == 1 {
        return a.clone();
    }
    let u = n as u64;
    if u % 2 == 1 || iv.lo.ge0() {
        return Some(Iv {
            lo: pow_e(ctx, &iv.lo, u, true),
            hi: pow_e(ctx, &iv.hi, u, false),
        });
    }
    if iv.hi.le0() {
        return Some(Iv {
            lo: pow_e(ctx, &iv.hi, u, true),
            hi: pow_e(ctx, &iv.lo, u, false),
        });
    }
    let m = if iv.lo.neg() > iv.hi { &iv.lo } else { &iv.hi };
    Some(Iv {
        lo: zero(),
        hi: pow_e(ctx, m, u, false),
    })
}

pub fn abs(a: &I) -> I {
    let iv = a.as_ref()?;
    if iv.lo.ge0() {
        return a.clone();
    }
    if iv.hi.le0() {
        return neg(a);
    }
    let nl = iv.lo.neg();
    Some(Iv {
        lo: zero(),
        hi: if nl > iv.hi { nl } else { iv.hi.clone() },
    })
}

pub fn min(a: &I, b: &I) -> I {
    let (a, b) = (a.as_ref()?, b.as_ref()?);
    Some(Iv {
        lo: a.lo.clone().min(b.lo.clone()),
        hi: a.hi.clone().min(b.hi.clone()),
    })
}

pub fn max(a: &I, b: &I) -> I {
    let (a, b) = (a.as_ref()?, b.as_ref()?);
    Some(Iv {
        lo: a.lo.clone().max(b.lo.clone()),
        hi: a.hi.clone().max(b.hi.clone()),
    })
}

pub fn sqrt(ctx: &Ctx, a: &I) -> I {
    let iv = a.as_ref()?;
    if iv.hi.lt0() {
        return None;
    }
    let lo = match &iv.lo {
        Fin(l) if l.is_positive() => Fin(isqrt(&(l << ctx.kk()))),
        PosInf => PosInf,
        _ => zero(),
    };
    let hi = match &iv.hi {
        Fin(h) => {
            let s = h << ctx.kk();
            let r = isqrt(&s);
            Fin(if &r * &r == s { r } else { r + 1u8 })
        }
        _ => PosInf,
    };
    Some(Iv { lo, hi })
}

// ---- transcendental functions ------------------------------------------------

fn guard(ctx: &Ctx, extra: i64) -> Rc<Ctx> {
    context(ctx.k + extra)
}

fn rescale(from: &Ctx, to: &Ctx, iv: &I) -> I {
    let iv = iv.as_ref()?;
    let s = (from.k - to.k) as u64;
    Some(Iv {
        lo: match &iv.lo {
            Fin(v) => Fin(shr(v, s)),
            e => e.clone(),
        },
        hi: match &iv.hi {
            Fin(v) => Fin(shr_up(v, s)),
            e => e.clone(),
        },
    })
}

fn lift(from: &Ctx, to: &Ctx, iv: &Iv) -> Iv {
    let s = (to.k - from.k) as u64;
    let up = |e: &E| match e {
        Fin(v) => Fin(v << s),
        e => e.clone(),
    };
    Iv {
        lo: up(&iv.lo),
        hi: up(&iv.hi),
    }
}

/// atan(1/n) at context g (n ≥ 2).
fn atan_inv(g: &Ctx, n: u32) -> (BigInt, BigInt) {
    let mut sum = BigInt::zero();
    let mut err = BigInt::one();
    let n = BigInt::from(n);
    let n2 = &n * &n;
    let mut j = 0u64;
    let mut denom = n.clone();
    loop {
        let term = &g.one / (BigInt::from(2 * j + 1) * &denom);
        if term.is_zero() {
            break;
        }
        if j % 2 == 0 {
            sum += term;
        } else {
            sum -= term;
        }
        err += 1u8;
        j += 1;
        denom *= &n2;
    }
    (&sum - &err - 1u8, &sum + &err + 1u8)
}

fn pi_pair(ctx: &Ctx) -> &(BigInt, BigInt) {
    ctx.pi.get_or_init(|| {
        let g = guard(ctx, 16);
        let a = atan_inv(&g, 5);
        let b = atan_inv(&g, 239);
        let lo = BigInt::from(16) * &a.0 - BigInt::from(4) * &b.1;
        let hi = BigInt::from(16) * &a.1 - BigInt::from(4) * &b.0;
        let r = rescale(
            &g,
            ctx,
            &Some(Iv {
                lo: Fin(lo),
                hi: Fin(hi),
            }),
        )
        .unwrap();
        (r.lo.fin().unwrap().clone(), r.hi.fin().unwrap().clone())
    })
}

pub fn pi_interval(ctx: &Ctx) -> Iv {
    let (lo, hi) = pi_pair(ctx);
    Iv {
        lo: Fin(lo.clone()),
        hi: Fin(hi.clone()),
    }
}

/// 2·atanh(tn/td) for 0 ≤ t ≤ 1/3.
fn two_atanh(g: &Ctx, tn: &BigInt, td: &BigInt) -> (BigInt, BigInt) {
    let mut sum = BigInt::zero();
    let mut err = BigInt::one();
    let mut num = tn.clone();
    let mut den = td.clone();
    let tn2 = tn * tn;
    let td2 = td * td;
    let mut j = 0u64;
    loop {
        let term = (&g.one * &num) / (&den * BigInt::from(2 * j + 1));
        if term.is_zero() {
            break;
        }
        sum += term;
        err += 1u8;
        num *= &tn2;
        den *= &td2;
        if j > 100000 {
            break;
        }
        j += 1;
    }
    (
        BigInt::from(2) * &sum - BigInt::from(2) * &err,
        BigInt::from(2) * (&sum + &err + 2u8),
    )
}

fn ln2(ctx: &Ctx) -> &(BigInt, BigInt) {
    ctx.ln2
        .get_or_init(|| two_atanh(ctx, &BigInt::one(), &BigInt::from(3)))
}

fn log_point(g: &Ctx, m: &BigInt) -> (BigInt, BigInt) {
    let bl = bit_length(m) as i64;
    let e = BigInt::from(bl - 1 - g.k);
    let fd = BigInt::one() << (bl - 1) as u64;
    let t = two_atanh(g, &(m - &fd), &(m + &fd));
    let l2 = ln2(g);
    let (elo, ehi) = if !e.is_negative() {
        (&e * &l2.0, &e * &l2.1)
    } else {
        (&e * &l2.1, &e * &l2.0)
    };
    (t.0 + elo, t.1 + ehi)
}

pub fn log(ctx: &Ctx, a: &I) -> I {
    let iv = a.as_ref()?;
    if iv.hi.le0() {
        return None;
    }
    let g = guard(ctx, 24);
    let up = lift(ctx, &g, iv);
    let lo = match &up.lo {
        Fin(v) if v.is_positive() => Fin(log_point(&g, v).0),
        PosInf => PosInf,
        _ => NegInf,
    };
    let hi = match &up.hi {
        Fin(v) => Fin(log_point(&g, v).1),
        _ => PosInf,
    };
    rescale(&g, ctx, &Some(Iv { lo, hi }))
}

fn exp_point(g: &Ctx, m: &BigInt) -> (BigInt, E) {
    let limit = BigInt::from(2000) * &g.one;
    if m > &limit {
        return (&g.one << 1000u32, PosInf);
    }
    if m < &-&limit {
        return (BigInt::zero(), Fin(BigInt::one()));
    }
    let mut s = 0u64;
    let half = &g.one >> 1u32;
    let mut am = m.abs();
    while am > half {
        am >>= 1u32;
        s += 1;
    }
    let den_shift = g.kk() + s;
    let mut sum = BigInt::zero();
    let mut err = BigInt::one();
    let mut num = BigInt::one();
    let mut fact = BigInt::one();
    let mut shift = 0u64;
    let mut j = 0u64;
    loop {
        if j > 0 {
            num *= m;
            fact *= j;
            shift += den_shift;
        }
        let full = shr(&(&num << g.kk()), shift);
        let term = &full / &fact; // truncates toward zero, like JS
        if term.is_zero() && j > 2 {
            break;
        }
        sum += term;
        err += 2u8;
        if j > 4000 {
            break;
        }
        j += 1;
    }
    let mut lo = &sum - &err - 2u8;
    let mut hi = &sum + &err + 2u8;
    if lo.is_negative() {
        lo = BigInt::zero();
    }
    for _ in 0..s {
        lo = shr(&(&lo * &lo), g.kk());
        hi = shr_up(&(&hi * &hi), g.kk());
    }
    (lo, Fin(hi))
}

fn bit_length_of(v: &E) -> i64 {
    v.fin().map_or(0, |x| bit_length(x) as i64)
}

pub fn exp(ctx: &Ctx, a: &I) -> I {
    let iv = a.as_ref()?;
    let mag = bit_length_of(&iv.lo).max(bit_length_of(&iv.hi)) - ctx.k;
    let g = guard(ctx, 24 + mag.max(0) * 2);
    let up = lift(ctx, &g, iv);
    let lo = match &up.lo {
        Fin(v) => Fin(exp_point(&g, v).0),
        _ => zero(),
    };
    let hi = match &up.hi {
        Fin(v) => exp_point(&g, v).1,
        _ => PosInf,
    };
    rescale(&g, ctx, &Some(Iv { lo, hi }))
}

/// Taylor enclosure of sin or cos on a narrow interval t with |t| ≤ 1.
fn sin_cos_small(g: &Ctx, lo: &BigInt, hi: &BigInt, sin: bool) -> (BigInt, BigInt) {
    let mid = shr(&(lo + hi), 1);
    let rad = hi - &mid + 1u8;
    let mut sum = BigInt::zero();
    let mut err = BigInt::one();
    let mut term = if sin { mid.clone() } else { g.one.clone() };
    let mut j: u64 = if sin { 1 } else { 0 };
    let m2 = &mid * &mid;
    let one2 = &g.one * &g.one;
    let mut i = 0u64;
    loop {
        if term.is_zero() {
            break;
        }
        if i % 2 == 0 {
            sum += &term;
        } else {
            sum -= &term;
        }
        err += 1u8;
        term = (&term * &m2) / (&one2 * BigInt::from(j + 1) * BigInt::from(j + 2));
        j += 2;
        if i > 2000 {
            break;
        }
        i += 1;
    }
    (&sum - &err - 1u8 - &rad, &sum + &err + 1u8 + &rad)
}

/// JS `Math.round`: halves round toward +∞.
fn js_round(v: f64) -> f64 {
    let f = v.floor();
    if v - f >= 0.5 { f + 1.0 } else { f }
}

fn sin_narrow(g: &Ctx, xlo: &BigInt, xhi: &BigInt) -> Result<(BigInt, BigInt), Overflow> {
    let (plo, phi) = pi_pair(g);
    let sh = if g.k > 40 { (g.k - 40) as u64 } else { 0 };
    let xf = to_f64(&shr(xlo, sh)) / 2f64.powi(g.k.min(40) as i32);
    let qf = js_round(xf / (std::f64::consts::PI / 2.0));
    if !qf.is_finite() {
        return Err(Overflow); // BigInt(Infinity) throws in JS
    }
    let q = <BigInt as num_traits::FromPrimitive>::from_f64(qf).ok_or(Overflow)?;
    let (qp_lo, qp_hi) = if !q.is_negative() {
        (&q * plo, &q * phi)
    } else {
        (&q * phi, &q * plo)
    };
    let two = BigInt::from(2);
    let tlo = xlo - ceil_div(&qp_hi, &two);
    let thi = xhi - floor_div(&qp_lo, &two);
    let r = q.mod_floor(&BigInt::from(4));
    let r: u8 = r.try_into().unwrap();
    Ok(match r {
        0 => sin_cos_small(g, &tlo, &thi, true),
        1 => sin_cos_small(g, &tlo, &thi, false),
        _ => {
            let v = sin_cos_small(g, &tlo, &thi, r == 2);
            (-v.1, -v.0)
        }
    })
}

pub fn sin(ctx: &Ctx, a: &I) -> Result<I, Overflow> {
    sin_shift(ctx, a, false)
}

pub fn cos(ctx: &Ctx, a: &I) -> Result<I, Overflow> {
    sin_shift(ctx, a, true)
}

fn sin_shift(ctx: &Ctx, a: &I, quarter: bool) -> Result<I, Overflow> {
    let Some(iv) = a.as_ref() else {
        return Ok(None);
    };
    let unit = || {
        Some(Iv {
            lo: Fin(-&ctx.one),
            hi: Fin(ctx.one.clone()),
        })
    };
    let (Fin(alo), Fin(ahi)) = (&iv.lo, &iv.hi) else {
        return Ok(unit());
    };
    if ahi - alo >= BigInt::from(7) * &ctx.one {
        return Ok(unit());
    }
    let g = guard(ctx, 24);
    let (plo, phi) = pi_pair(&g).clone();
    let s = (g.k - ctx.k) as u64;
    let (mut xlo, mut xhi) = (alo << s, ahi << s);
    if quarter {
        xlo += shr(&plo, 1);
        xhi += ceil_div(&phi, &BigInt::from(2));
    }
    let lo = sin_narrow(&g, &xlo, &xlo)?;
    let hi = sin_narrow(&g, &xhi, &xhi)?;
    let mut rlo = lo.0.min(hi.0.clone());
    let mut rhi = lo.1.max(hi.1);
    let scale = 2f64.powi(g.k as i32);
    let xl = to_f64(&xlo) / scale;
    let xh = to_f64(&xhi) / scale;
    let hp = std::f64::consts::PI / 2.0;
    let j0 = ((xl - hp) / std::f64::consts::PI).floor() - 1.0;
    let j1 = ((xh - hp) / std::f64::consts::PI).ceil() + 1.0;
    if j0.abs() > 9.0e15 || j1.abs() > 9.0e15 {
        return Err(Overflow); // the JS loop cannot step past 2^53 either
    }
    if j0.is_finite() && j1.is_finite() {
        let two = BigInt::from(2);
        let mut j = j0 as i64;
        while (j as f64) <= j1 {
            let jj = BigInt::from(2 * j + 1);
            let (c_lo, c_hi) = if !jj.is_negative() {
                (
                    floor_div(&(&jj * &plo), &two),
                    ceil_div(&(&jj * &phi), &two),
                )
            } else {
                (
                    floor_div(&(&jj * &phi), &two),
                    ceil_div(&(&jj * &plo), &two),
                )
            };
            if !(c_hi < xlo || c_lo > xhi) {
                let v = if j.rem_euclid(2) == 0 {
                    g.one.clone()
                } else {
                    -&g.one
                };
                if v < rlo {
                    rlo = v.clone();
                }
                if v > rhi {
                    rhi = v;
                }
            }
            j += 1;
        }
    }
    // clamp to [-1, 1]
    let m1 = -&g.one;
    if rlo < m1 {
        rlo = m1;
    }
    if rhi > g.one {
        rhi = g.one.clone();
    }
    Ok(rescale(
        &g,
        ctx,
        &Some(Iv {
            lo: Fin(rlo),
            hi: Fin(rhi),
        }),
    ))
}

/// a^c for an exponent interval c > 0 on the domain a ≥ 0 (0^c = 0); the
/// ends of a give the ends, each as exp(c · ln m) (graphics `powPos`).
pub fn pow_pos(ctx: &Ctx, a: &I, c: &I) -> I {
    let (iv, _) = (a.as_ref()?, c.as_ref()?);
    if iv.hi.lt0() {
        return None;
    }
    let at = |m: &E| {
        exp(
            ctx,
            &mul(
                ctx,
                c,
                &log(
                    ctx,
                    &Some(Iv {
                        lo: m.clone(),
                        hi: m.clone(),
                    }),
                ),
            ),
        )
    };
    let lo = if iv.lo.gt0() { at(&iv.lo)?.lo } else { zero() };
    let hi = if !iv.hi.is_fin() {
        PosInf
    } else if iv.hi.is_zero() {
        zero()
    } else {
        at(&iv.hi)?.hi
    };
    Some(Iv { lo, hi })
}

pub fn tan(ctx: &Ctx, a: &I) -> Result<I, Overflow> {
    let s = sin(ctx, a)?;
    let c = cos(ctx, a)?;
    Ok(div(ctx, &s, &c))
}
