//! # ratmath (Rust)
//!
//! Rust kernels for the hot paths of the RatMath JavaScript libraries. The JS
//! libraries stay the reference implementation; each kernel here is a faithful
//! port whose results are identical (checked by `tests/crosscheck.test.js`), so
//! a library may hand an operation to Rust whenever that is faster and get the
//! same answer.
//!
//! * [`matrix`]: exact rational linear algebra (ratmath/linalg `exact.js`):
//!   Bareiss RREF/solve/inverse/rank and the multi-modular determinant.
//! * [`dyadic`]: certified interval arithmetic on dyadic fixed point
//!   (ratmath/graphics `interval.js`), including exp, ln, sin, cos, tan.
//! * [`graph`]: certified pixel-column enclosures of y = f(x)
//!   (ratmath/graphics `certifyGraph`).
//! * [`ffi`] / [`wire`]: one C ABI for WASM and native (Swift, later Android).

pub mod dyadic;
pub mod ffi;
pub mod graph;
pub mod matrix;
pub mod rat;
pub mod wire;
