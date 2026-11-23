//! The text wire format shared by the C ABI (WASM and Swift). Tokens are
//! separated by single spaces; every integer, big or small, is written in
//! hexadecimal with an optional leading `-` (what JS `BigInt#toString(16)`
//! and `Number#toString(16)` produce). See WIRE.md for each operation's layout.

use crate::graph::{self, Op, Program, Viewport};
use crate::matrix::{self, SolveStatus};
use crate::rat::{Pair, reduce};
use num_bigint::BigInt;
use num_traits::{ToPrimitive, Zero};
use std::fmt::Write;

pub const VERSION: &str = env!("CARGO_PKG_VERSION");

/// Operation codes for `rm_call`.
pub mod ops {
    pub const VERSION: u32 = 0;
    pub const DET: u32 = 1;
    pub const DET_BAREISS: u32 = 2;
    pub const DET_MODULAR: u32 = 3;
    pub const RANK: u32 = 4;
    pub const RREF: u32 = 5;
    pub const SOLVE: u32 = 6;
    pub const INVERSE: u32 = 7;
    pub const NONSINGULAR: u32 = 8;
    pub const GRAPH_COLUMNS: u32 = 20;
}

pub struct Reader<'a> {
    toks: std::str::SplitAsciiWhitespace<'a>,
}

impl<'a> Reader<'a> {
    pub fn new(s: &'a str) -> Self {
        Reader {
            toks: s.split_ascii_whitespace(),
        }
    }
    fn tok(&mut self) -> Result<&'a str, String> {
        self.toks
            .next()
            .ok_or_else(|| "unexpected end of input".to_string())
    }
    pub fn big(&mut self) -> Result<BigInt, String> {
        let t = self.tok()?;
        let (neg, digits) = match t.strip_prefix('-') {
            Some(r) => (true, r),
            None => (false, t),
        };
        let v = BigInt::parse_bytes(digits.as_bytes(), 16)
            .ok_or_else(|| format!("bad hex integer {t:?}"))?;
        Ok(if neg { -v } else { v })
    }
    pub fn int(&mut self) -> Result<i64, String> {
        self.big()?
            .to_i64()
            .ok_or_else(|| "integer out of range".to_string())
    }
    pub fn usize(&mut self) -> Result<usize, String> {
        let v = self.int()?;
        usize::try_from(v).map_err(|_| format!("expected a count, got {v}"))
    }
    pub fn pair(&mut self) -> Result<Pair, String> {
        let n = self.big()?;
        let d = self.big()?;
        reduce(n, d).map_err(|_| "zero denominator".to_string())
    }
    pub fn matrix(&mut self, rows: usize, cols: usize) -> Result<Vec<Vec<Pair>>, String> {
        (0..rows)
            .map(|_| (0..cols).map(|_| self.pair()).collect())
            .collect()
    }
    /// `rows cols` followed by the entries.
    pub fn sized_matrix(&mut self) -> Result<Vec<Vec<Pair>>, String> {
        let r = self.usize()?;
        let c = self.usize()?;
        self.matrix(r, c)
    }
}

#[derive(Default)]
pub struct Writer {
    pub out: String,
}

impl Writer {
    fn sep(&mut self) {
        if !self.out.is_empty() {
            self.out.push(' ');
        }
    }
    pub fn big(&mut self, v: &BigInt) {
        self.sep();
        let _ = write!(self.out, "{v:x}");
    }
    pub fn int(&mut self, v: i64) {
        self.sep();
        if v < 0 {
            let _ = write!(self.out, "-{:x}", v.unsigned_abs());
        } else {
            let _ = write!(self.out, "{v:x}");
        }
    }
    pub fn pair(&mut self, p: &Pair) {
        self.big(&p.0);
        self.big(&p.1);
    }
    pub fn pairs(&mut self, m: &[Vec<Pair>]) {
        for row in m {
            for p in row {
                self.pair(p);
            }
        }
    }
}

fn square(m: &[Vec<Pair>], what: &str) -> Result<(), String> {
    if m.iter().any(|r| r.len() != m.len()) {
        return Err(format!("{what} needs a square matrix"));
    }
    Ok(())
}

/// Run one operation on wire-format input; returns wire-format output.
pub fn call(op: u32, input: &str) -> Result<String, String> {
    let mut r = Reader::new(input);
    let mut w = Writer::default();
    match op {
        ops::VERSION => w.out.push_str(VERSION),
        ops::DET | ops::DET_BAREISS | ops::DET_MODULAR => {
            let m = r.sized_matrix()?;
            square(&m, "det")?;
            let p = match op {
                ops::DET => matrix::determinant(&m),
                ops::DET_BAREISS => matrix::determinant_bareiss(&m),
                _ => matrix::determinant_modular(&m),
            };
            w.pair(&p);
        }
        ops::RANK => {
            let m = r.sized_matrix()?;
            let cols = m.first().map_or(0, |r| r.len());
            w.int(matrix::rank(&m, cols) as i64);
        }
        ops::RREF => {
            // rows cols columns entries -> npivots pivots… entries…
            let rows = r.usize()?;
            let cols = r.usize()?;
            let columns = r.usize()?.min(cols);
            let m = r.matrix(rows, cols)?;
            let (out, pivots) = matrix::rref_pairs(&m, columns);
            w.int(pivots.len() as i64);
            for p in &pivots {
                w.int(*p as i64);
            }
            w.pairs(&out);
        }
        ops::SOLVE => {
            // m n k A… B… -> status hasX X… nbasis basis…
            let m = r.usize()?;
            let n = r.usize()?;
            let k = r.usize()?;
            let a = r.matrix(m, n)?;
            let b = r.matrix(m, k)?;
            let s = matrix::solve_pairs(&a, &b);
            w.int(match s.status {
                SolveStatus::Unique => 0,
                SolveStatus::Infinite => 1,
                SolveStatus::None => 2,
            });
            match &s.x {
                Some(x) => {
                    w.int(1);
                    w.pairs(x);
                }
                None => w.int(0),
            }
            w.int(s.basis.len() as i64);
            for v in &s.basis {
                for p in v {
                    w.pair(p);
                }
            }
        }
        ops::INVERSE => {
            let m = r.sized_matrix()?;
            square(&m, "inverse")?;
            match matrix::inverse_pairs(&m) {
                Some(inv) => {
                    w.int(1);
                    w.pairs(&inv);
                }
                None => w.int(0),
            }
        }
        ops::NONSINGULAR => {
            let tries = r.usize()?;
            let m = r.sized_matrix()?;
            square(&m, "nonsingular")?;
            w.int(matrix::provably_nonsingular(&m, tries) as i64);
        }
        ops::GRAPH_COLUMNS => {
            let k = r.int()?;
            let max_depth = r.int()? as u32;
            let width = r.int()? as u32;
            let height = r.int()? as u32;
            let from = r.int()? as u32;
            let to = r.int()? as u32;
            let mut rat = || -> Result<(BigInt, BigInt), String> {
                let n = r.big()?;
                let d = r.big()?;
                if d.is_zero() {
                    return Err("zero denominator".into());
                }
                Ok((n, d))
            };
            let xmin = rat()?;
            let xmax = rat()?;
            let ymin = rat()?;
            let ymax = rat()?;
            let np = r.usize()?;
            let params = (0..np).map(|_| r.pair()).collect::<Result<Vec<_>, _>>()?;
            let nc = r.usize()?;
            let consts = (0..nc).map(|_| r.pair()).collect::<Result<Vec<_>, _>>()?;
            let ncode = r.usize()?;
            let mut code = Vec::with_capacity(ncode);
            for _ in 0..ncode {
                let c = r.int()?;
                let mut arg_err = None;
                let op = Op::decode(c, || match r.int() {
                    Ok(v) => Some(v),
                    Err(e) => {
                        arg_err = Some(e);
                        None
                    }
                });
                if let Some(e) = arg_err {
                    return Err(e);
                }
                code.push(op.ok_or_else(|| format!("unknown opcode {c}"))?);
            }
            if width == 0 || height == 0 {
                return Err("viewport width and height must be positive".into());
            }
            let vp = Viewport {
                xmin,
                xmax,
                ymin,
                ymax,
                width,
                height,
            };
            let prog = Program { code, consts };
            let cols = graph::certify_columns(&prog, &params, &vp, k, max_depth, from, to);
            w.int(cols.len() as i64);
            for c in &cols {
                w.int(c.tight as i64);
                w.int(c.evaluations as i64);
                w.int(c.spans.len() as i64);
                for s in &c.spans {
                    w.int(*s);
                }
            }
        }
        _ => return Err(format!("unknown operation {op}")),
    }
    Ok(w.out)
}
