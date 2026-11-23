//! Certified graphs of y = f(x): the per-pixel-column bisection of
//! ratmath/graphics `certifyGraph` / `certifyColumn`, evaluating a small stack
//! program (compiled from the graphics AST by `js/bytecode.js`) over
//! `dyadic` intervals. Output matches the JS frame column for column.

use crate::dyadic::{self as d, Ctx, E, I, Iv, Overflow};
use num_bigint::BigInt;
use num_traits::ToPrimitive;

/// One stack-machine instruction. Operands are popped right-to-left.
#[derive(Debug, Clone, PartialEq)]
pub enum Op {
    X,
    Param(usize),
    Const(usize),
    Pi,
    E,
    Add,
    Sub,
    Mul,
    Div,
    Neg,
    PowInt(i64),
    Pow,
    PowPos,
    Recip,
    Sin,
    Cos,
    Tan,
    Exp,
    Ln,
    Sqrt,
    Abs,
    Min,
    Max,
}

impl Op {
    /// Numeric opcodes used on the wire (see WIRE.md).
    pub fn decode(code: i64, arg: impl FnOnce() -> Option<i64>) -> Option<Op> {
        Some(match code {
            0 => Op::X,
            1 => Op::Param(arg()? as usize),
            2 => Op::Const(arg()? as usize),
            3 => Op::Pi,
            4 => Op::E,
            10 => Op::Add,
            11 => Op::Sub,
            12 => Op::Mul,
            13 => Op::Div,
            14 => Op::Neg,
            15 => Op::PowInt(arg()?),
            16 => Op::Pow,
            17 => Op::PowPos,
            18 => Op::Recip,
            20 => Op::Sin,
            21 => Op::Cos,
            22 => Op::Tan,
            23 => Op::Exp,
            24 => Op::Ln,
            25 => Op::Sqrt,
            26 => Op::Abs,
            27 => Op::Min,
            28 => Op::Max,
            _ => return None,
        })
    }
}

pub struct Program {
    pub code: Vec<Op>,
    pub consts: Vec<(BigInt, BigInt)>,
}

/// A program bound to one precision context, with its literals, parameters
/// and `e` converted once.
pub struct Bound<'a> {
    prog: &'a Program,
    ctx: std::rc::Rc<Ctx>,
    consts: Vec<Iv>,
    params: Vec<Iv>,
    e: std::cell::OnceCell<I>,
}

impl<'a> Bound<'a> {
    pub fn new(prog: &'a Program, k: i64, params: &[(BigInt, BigInt)]) -> Self {
        let ctx = d::context(k);
        let consts = prog
            .consts
            .iter()
            .map(|(n, q)| d::from_rational(&ctx, n, q))
            .collect();
        let params = params
            .iter()
            .map(|(n, q)| d::from_rational(&ctx, n, q))
            .collect();
        Bound {
            prog,
            ctx,
            consts,
            params,
            e: std::cell::OnceCell::new(),
        }
    }

    pub fn ctx(&self) -> &Ctx {
        &self.ctx
    }

    /// Enclosure of f over X. Err where the JS kernel would throw a RangeError.
    pub fn eval(&self, x: &Iv) -> Result<I, Overflow> {
        let c = &*self.ctx;
        let mut st: Vec<I> = Vec::with_capacity(8);
        for op in &self.prog.code {
            let v = match op {
                Op::X => Some(x.clone()),
                Op::Param(i) => Some(self.params.get(*i).ok_or(Overflow)?.clone()),
                Op::Const(i) => Some(self.consts.get(*i).ok_or(Overflow)?.clone()),
                Op::Pi => Some(d::pi_interval(c)),
                Op::E => self
                    .e
                    .get_or_init(|| {
                        d::exp(
                            c,
                            &Some(Iv {
                                lo: E::Fin(c.one.clone()),
                                hi: E::Fin(c.one.clone()),
                            }),
                        )
                    })
                    .clone(),
                Op::Neg => d::neg(&pop(&mut st)?),
                Op::PowInt(k) => d::pow_int(c, &pop(&mut st)?, *k),
                Op::Sin => d::sin(c, &pop(&mut st)?)?,
                Op::Cos => d::cos(c, &pop(&mut st)?)?,
                Op::Tan => d::tan(c, &pop(&mut st)?)?,
                Op::Exp => d::exp(c, &pop(&mut st)?),
                Op::Ln => d::log(c, &pop(&mut st)?),
                Op::Sqrt => d::sqrt(c, &pop(&mut st)?),
                Op::Abs => d::abs(&pop(&mut st)?),
                Op::Recip => d::recip(c, &pop(&mut st)?),
                Op::Add
                | Op::Sub
                | Op::Mul
                | Op::Div
                | Op::Pow
                | Op::PowPos
                | Op::Min
                | Op::Max => {
                    let b = pop(&mut st)?;
                    let a = pop(&mut st)?;
                    match op {
                        Op::Add => d::add(&a, &b),
                        Op::Sub => d::sub(&a, &b),
                        Op::Mul => d::mul(c, &a, &b),
                        Op::Div => d::div(c, &a, &b),
                        // base^e = exp(e · ln(base)), as the JS compiler does
                        Op::Pow => d::exp(c, &d::mul(c, &b, &d::log(c, &a))),
                        Op::PowPos => d::pow_pos(c, &a, &b),
                        Op::Min => d::min(&a, &b),
                        _ => d::max(&a, &b),
                    }
                }
            };
            st.push(v);
        }
        if st.len() != 1 {
            return Err(Overflow);
        }
        Ok(st.pop().unwrap())
    }
}

fn pop(st: &mut Vec<I>) -> Result<I, Overflow> {
    st.pop().ok_or(Overflow)
}

/// The exact viewport, as rationals (any representation; only floors are taken).
pub struct Viewport {
    pub xmin: (BigInt, BigInt),
    pub xmax: (BigInt, BigInt),
    pub ymin: (BigInt, BigInt),
    pub ymax: (BigInt, BigInt),
    pub width: u32,
    pub height: u32,
}

/// Exact x-interval of pixel column i at ctx, rounded outward (`columnInterval`).
pub fn column_interval(ctx: &Ctx, vp: &Viewport, i: u32) -> Iv {
    let (an, ad) = &vp.xmin;
    let (bn, bd) = &vp.xmax;
    let w = BigInt::from(vp.width);
    // xmin + (xmax - xmin) * i / width, as num / den with den > 0
    let den = ad * bd * &w;
    let at = |i: u32| an * bd * &w + (bn * ad - an * bd) * BigInt::from(i);
    let k = ctx.k as u64;
    Iv {
        lo: E::Fin(d::floor_div(&(at(i) << k), &den)),
        hi: E::Fin(d::ceil_div(&(at(i + 1) << k), &den)),
    }
}

/// Pixel rows [top, bottom] touched by the y-interval, clipped, or None (`rowSpan`).
pub fn row_span(ctx: &Ctx, vp: &Viewport, iv: &I) -> Option<(i64, i64)> {
    let iv = iv.as_ref()?;
    let h = vp.height as i64;
    let (a, b) = &vp.ymax;
    let (yn, yd) = &vp.ymin;
    // span = ymax - ymin = c / dd
    let c = a * yd - yn * b;
    let dd = b * yd;
    let hb = BigInt::from(h);
    let row_of = |y: &E| -> i64 {
        match y {
            E::PosInf => -1,
            E::NegInf => h,
            E::Fin(y) => {
                let num = (a * &ctx.one - y * b) * &hb * &dd;
                let den = b * &c * &ctx.one;
                let r = d::floor_div(&num, &den);
                if r < BigInt::from(-1) {
                    -1
                } else if r > hb {
                    h
                } else {
                    r.to_i64().unwrap()
                }
            }
        }
    };
    let mut top = row_of(&iv.hi);
    let mut bot = row_of(&iv.lo);
    if bot < 0 || top >= h {
        return None;
    }
    if top < 0 {
        top = 0;
    }
    if bot >= h {
        bot = h - 1;
    }
    Some((top, bot))
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Column {
    /// Merged inclusive row ranges, flattened: top, bottom, top, bottom, …
    pub spans: Vec<i64>,
    pub tight: bool,
    pub evaluations: u64,
}

/// Certified enclosure of the graph over one pixel column (`certifyColumn`).
pub fn certify_column(f: &Bound, vp: &Viewport, x: &Iv, max_depth: u32) -> Column {
    let ctx = f.ctx();
    let mut spans: Vec<i64> = Vec::new();
    let mut tight = true;
    let mut evaluations = 0u64;
    let span = |iv: &Iv, evaluations: &mut u64| -> Option<(i64, i64)> {
        *evaluations += 1;
        let y = f.eval(iv).unwrap_or_else(|_| Some(d::whole()));
        row_span(ctx, vp, &y)
    };
    fn lohi(iv: &Iv) -> (&BigInt, &BigInt) {
        (iv.lo.fin().unwrap(), iv.hi.fin().unwrap())
    }
    // explicit stack in the same depth-first order as the JS recursion
    let first = span(x, &mut evaluations);
    let mut stack: Vec<(Iv, Option<(i64, i64)>, u32)> = vec![(x.clone(), first, 0)];
    while let Some((iv, rows, depth)) = stack.pop() {
        let Some(rows) = rows else { continue };
        let count = rows.1 - rows.0 + 1;
        if count <= 2 {
            spans.extend([rows.0, rows.1]);
            continue;
        }
        let (lo, hi) = lohi(&iv);
        if depth >= max_depth || hi - lo < BigInt::from(2) {
            tight = false;
            spans.extend([rows.0, rows.1]);
            continue;
        }
        let mid = d::floor_div(&(lo + hi), &BigInt::from(2));
        let l = Iv {
            lo: iv.lo.clone(),
            hi: E::Fin(mid.clone()),
        };
        let r = Iv {
            lo: E::Fin(mid),
            hi: iv.hi.clone(),
        };
        let rl = span(&l, &mut evaluations);
        let rr = span(&r, &mut evaluations);
        let covered = covered_rows(rl, rr);
        if covered >= count - 1 && depth >= 1 {
            if let Some(a) = rl {
                spans.extend([a.0, a.1]);
            }
            if let Some(b) = rr {
                spans.extend([b.0, b.1]);
            }
            continue;
        }
        // visit(L) then visit(R): push R first
        stack.push((r, rr, depth + 1));
        stack.push((l, rl, depth + 1));
    }
    Column {
        spans: merge_spans(spans),
        tight,
        evaluations,
    }
}

fn covered_rows(a: Option<(i64, i64)>, b: Option<(i64, i64)>) -> i64 {
    match (a, b) {
        (None, None) => 0,
        (None, Some(b)) => b.1 - b.0 + 1,
        (Some(a), None) => a.1 - a.0 + 1,
        (Some(a), Some(b)) => {
            let lo = a.0.min(b.0);
            let hi = a.1.max(b.1);
            let gap = 0.max(a.0.max(b.0) - a.1.min(b.1) - 1);
            hi - lo + 1 - gap
        }
    }
}

fn merge_spans(flat: Vec<i64>) -> Vec<i64> {
    if flat.len() <= 2 {
        return flat;
    }
    let mut pairs: Vec<(i64, i64)> = flat.chunks_exact(2).map(|p| (p[0], p[1])).collect();
    pairs.sort_by_key(|p| p.0); // stable, like Array.prototype.sort
    let mut out = vec![pairs[0].0, pairs[0].1];
    for &(t, b) in &pairs[1..] {
        let last = out.len() - 1;
        if t <= out[last] + 1 {
            if b > out[last] {
                out[last] = b;
            }
        } else {
            out.extend([t, b]);
        }
    }
    out
}

/// Certify pixel columns `from..to` of the graph (`certifyGraph` without the
/// progress yields: the caller picks the slice size).
pub fn certify_columns(
    prog: &Program,
    params: &[(BigInt, BigInt)],
    vp: &Viewport,
    k: i64,
    max_depth: u32,
    from: u32,
    to: u32,
) -> Vec<Column> {
    let f = Bound::new(prog, k, params);
    (from..to.min(vp.width))
        .map(|i| {
            let x = column_interval(f.ctx(), vp, i);
            certify_column(&f, vp, &x, max_depth)
        })
        .collect()
}
