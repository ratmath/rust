# Opt-in backend interface

The RatMath JS libraries are the reference implementation. A library may hand
an operation to the Rust kernels when that is faster, and must get **the
identical result** (not just an equal value: same pivots, same row order, same
pixel spans). `tests/crosscheck.test.js` checks this against the JS code.

## Contract

```js
import { loadRustBackend } from "@ratmath/rust";
const rust = await loadRustBackend();       // WASM next to js/index.js, or pass bytes/URL/Module
```

A backend object has:

| member | meaning |
|---|---|
| `name`, `version` | `"rust-wasm"`, crate version |
| `wants(op, size)` | true when the policy says Rust pays for this op at this size |
| `det(pairs)`, `detBareiss`, `detModular` | `[n, d]` |
| `rank(pairs)` | number |
| `rref(pairs, columns)` | `{ rows, pivots }` = linalg `rrefPairs` |
| `solve(a, b)` | `{ status, x, basis }` = linalg `solvePairs` |
| `inverse(pairs)` | rows or `null` = linalg `inversePairs` |
| `provablyNonsingular(pairs, tries)` | boolean, never wrong when true |
| `graphColumns({ program, params, viewport, k, maxDepth, from, to })` | columns = graphics `certifyGraph` frame columns |

Matrices are rows of `[numerator, denominator]` BigInt pairs, reduced with a
positive denominator: the shape linalg already uses internally.

Calls are synchronous. They give no progress ticks, so a library should call
them inside its exact-stage generator on slices it chooses (for graphs, a
range of columns per call via `from`/`to`), or from a worker.

## How a library opts in

Each library keeps one module-level slot and checks it inside the exact stage;
nothing changes for callers that never set it.

```js
// linalg src/backend.js (suggested)
let backend = null;
export function setBackend(b) { backend = b; }   // e.g. setBackend(await loadRustBackend())
export function viaBackend(op, size) { return backend?.wants(op, size) ? backend : null; }

// linalg ops.js, det's exact stage
function* () {
  const fast = viaBackend("det", A.rows);
  if (fast) return { value: toRational(fast.det(A.pairs())), via: "rust" };
  ...existing JS path...
}
```

graphics does the same in `certifyGraph`: when `viaBackend("graphColumns",
vp.width)` is set, compile once with `compileProgram(ast, variable)` and fetch
columns in slices of `partialEvery`, yielding the same progress ticks.

`policy` thresholds can be overridden per app:
`loadRustBackend(src, { policy: { det: 16, graphColumns: 1 } })`.

## Native (Swift, later Android)

The same operations are exported through the C ABI in `include/ratmath.h`;
`swift/RatMath.swift` sketches the wrapper. A native app that runs the JS
libraries in JavaScriptCore can route the same calls to the native kernels
instead of WASM.
