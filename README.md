# ratmath/rust

Rust kernels for the hot exact arithmetic in the RatMath JavaScript libraries,
callable from JS through WebAssembly and from Swift (and later Kotlin) through
one C ABI.

The JS libraries ([linalg](https://github.com/ratmath/linalg),
[graphics](https://github.com/ratmath/graphics), analytic) remain the
reference implementation. Each kernel here is a faithful port whose results are
**identical** to the JS (same pivots, same row order, same certified pixel
spans), so a library can hand an operation to Rust whenever the backend says it
pays, and nothing downstream can tell the difference except the clock.

## What is in it

| module | ports | what |
|---|---|---|
| `matrix` | linalg `src/exact.js` | Bareiss RREF / rank / solve / inverse / null space; Bareiss and multi-modular determinant; modular nonsingularity proof. Small matrices run in `i128` and fall back to BigInt on overflow. |
| `dyadic` | graphics `src/interval.js` | certified interval arithmetic on dyadic fixed point: + − × ÷, integer powers, sqrt, abs, min, max, exp, ln, sin, cos, tan, π |
| `graph` | graphics `certifyGraph` | per-pixel-column bisection producing the certified cell frame of y = f(x), from a small stack program |
| `ffi`, `wire` | — | the C ABI (`include/ratmath.h`) and its text wire format ([WIRE.md](WIRE.md)) |
| `js/` | — | `@ratmath/rust`: WASM loader, the opt-in backend object, the graphics-AST → program compiler ([BACKEND.md](BACKEND.md)) |

## Use from JS

```js
import { loadRustBackend } from "@ratmath/rust";
const rust = await loadRustBackend();
rust.det([[[1n, 2n], [1n, 3n]], [[1n, 4n], [1n, 5n]]]);   // [1n, 60n]
if (rust.wants("det", n)) …                                  // policy check
```

See [BACKEND.md](BACKEND.md) for how linalg and graphics opt in with one
`setBackend` slot.

## Use from Swift (Apple first)

`scripts/build-apple.sh` (on a Mac) builds `RatMathCore.xcframework` for iOS
devices, the simulator and macOS; `swift/RatMath.swift` sketches the wrapper.
This has not been built on a Mac yet; the C ABI it wraps is the same one the
WASM tests exercise.

## Speed

`bun bench/bench.js` times the JS reference (Bun runs JavaScriptCore, the
engine in Safari and iOS) against the Rust kernels through WASM, including the
cost of crossing. `cargo run --release --example bench` then times the same
inputs natively, which is what a Swift app calling the C ABI gets. Random
matrices with entries a/b, |a| ≤ 10, b ≤ 6; graphs are 300×200 certified frames.
One run on a 2.1 GHz Xeon (cloud container, noisy):

| operation | size | JS | Rust WASM | Rust native |
|---|---|---|---|---|
| det | 16 | 409 µs | 286 µs | 176 µs |
| det | 32 | 1.8 ms | 1.4 ms | 0.77 ms |
| det | 64 | 17.4 ms | 11.5 ms | 5.4 ms |
| det | 100 | 77.5 ms | 54.1 ms | 23.1 ms |
| solve | 16 | 454 µs | 941 µs | 480 µs |
| solve | 48 | 14.4 ms | 31.4 ms | 14.9 ms |
| inverse | 32 | 17.0 ms | 28.9 ms | 18.8 ms |
| rref | 48 | 17.6 ms | 47.3 ms | 17.3 ms |
| graph `x^2 - 2` | 300×200 | 2.6 ms | 1.1 ms | 0.56 ms |
| graph `sin(x)` | 300×200 | 3.3 ms | 3.2 ms | 1.9 ms |
| graph `a*sin(3x)/x + exp(-x^2)` | 300×200 | 17.7 ms | 16.4 ms | 13.2 ms |

What this says:

* **Word-size algorithms win.** The multi-modular determinant (31-bit primes,
  Barrett reduction) is 3–4× faster natively and about 1.4× in WASM. Certified
  graphs are 1.3–4× faster natively.
* **BigInt elimination does not.** Fraction-free elimination spends its time
  in big-integer multiply and divide, where JavaScriptCore's BigInt is as fast
  as Rust's `num-bigint` natively and faster than it in WASM (32-bit limbs
  there). So `DEFAULT_POLICY` sends only `det` (n ≥ 32) to WASM; rref, solve,
  inverse and graphs stay in JS unless a caller opts in.
* The way to make solve/inverse fast is algorithmic, not a language change:
  modular/p-adic solving (Dixon lifting with rational reconstruction) moves
  them onto machine words as the determinant already is. That is the next
  kernel to add here.

## Tests

```
cargo test --release            # Rust unit tests
bun run test                    # builds the WASM, then cross-checks against JS
```

The cross-check needs checkouts of ratmath/linalg and ratmath/graphics next to
this repo (`../linalg`, `../graphics`, with `bun install`), or set
`RATMATH_LINALG` / `RATMATH_GRAPHICS`. It compares 120 random matrices
(square, rectangular, singular, huge entries) across det/rref/rank/solve/
inverse, 17 hand-picked graphs, and 40 random expressions on four viewports
(including tiny and far-from-origin ones). A library whose code is missing is
skipped.

`js/ratmath.wasm` is committed so JS users need no Rust toolchain;
`bun run build:wasm` rebuilds it.

## License

MIT
