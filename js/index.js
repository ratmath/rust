// @ratmath/rust: load the Rust kernels (WebAssembly) and expose them as an
// opt-in backend for the RatMath JS libraries. Results are identical to the JS
// reference implementations; the backend only changes how fast they arrive.
//
//   import { loadRustBackend } from "@ratmath/rust";
//   const rust = await loadRustBackend();          // or pass bytes / a URL / a Module
//   if (rust.wants("det", n)) value = rust.det(pairs);
//
// All matrices are rows of [numerator, denominator] BigInt pairs (reduced,
// denominator > 0), the same shape ratmath/linalg uses internally.

import { Tokens, encodePairs, hex } from "./wire.js";

export { compileProgram, OP } from "./bytecode.js";

const OPS = {
  version: 0, det: 1, detBareiss: 2, detModular: 3, rank: 4, rref: 5,
  solve: 6, inverse: 7, nonsingular: 8, graphColumns: 20,
};

/**
 * Default dispatch policy: the smallest size at which handing an operation to
 * the WASM kernels pays, including the cost of crossing (BigInts to text and
 * back). From `bun bench/bench.js` (see README): the multi-modular
 * determinant wins from about 32×32; BigInt elimination (rref, solve,
 * inverse) does not beat the engine's own BigInt in WASM, and graphs are about
 * even, so those stay in JS unless a caller opts in. `size` is the matrix
 * dimension, or the number of pixel columns for graphColumns.
 */
export const DEFAULT_POLICY = {
  det: 32,
  rank: Infinity,
  rref: Infinity,
  solve: Infinity,
  inverse: Infinity,
  graphColumns: Infinity,
};

async function readSource(source) {
  if (source instanceof WebAssembly.Module) return source;
  if (source instanceof ArrayBuffer || ArrayBuffer.isView(source)) return source;
  if (typeof source === "string" || source instanceof URL) {
    const url = source instanceof URL ? source : new URL(source, import.meta.url);
    if (url.protocol === "file:") {
      const { readFile } = await import("node:fs/promises");
      return readFile(url);
    }
    const res = await fetch(url);
    return res.arrayBuffer();
  }
  throw new TypeError("loadRustBackend: pass wasm bytes, a URL, or a WebAssembly.Module");
}

/** Instantiate the WASM kernels. Default source: ./ratmath.wasm next to this file. */
export async function loadRustBackend(source = new URL("./ratmath.wasm", import.meta.url), { policy = {} } = {}) {
  const bytes = await readSource(source);
  const { instance } = bytes instanceof WebAssembly.Module
    ? { instance: await WebAssembly.instantiate(bytes, {}) }
    : await WebAssembly.instantiate(bytes, {});
  return makeBackend(instance.exports, { ...DEFAULT_POLICY, ...policy });
}

const encoder = new TextEncoder();
const decoder = new TextDecoder();

function makeBackend(x, policy) {
  function call(op, text) {
    if (backend._onCall) backend._onCall(op, text);
    const input = encoder.encode(text);
    const len = input.length;
    const ptr = x.rm_alloc(len);
    new Uint8Array(x.memory.buffer, ptr, len).set(input);
    const status = x.rm_call(op, ptr, len);
    x.rm_free(ptr, len);
    const out = decoder.decode(new Uint8Array(x.memory.buffer, x.rm_out_ptr(), x.rm_out_len()));
    if (status !== 0) throw new Error(`ratmath-rust: ${out}`);
    return out;
  }

  const sized = (rows) => {
    const out = [hex(rows.length), hex(rows[0]?.length ?? 0)];
    encodePairs(out, rows);
    return out.join(" ");
  };

  const backend = {
    name: "rust-wasm",
    version: null,
    policy,

    /** Should operation `op` of this size go to Rust? */
    wants(op, size) {
      const t = policy[op];
      return t !== undefined && size >= t;
    },

    /** Exact determinant as [n, d] (Bareiss for n ≤ 3, multi-modular above, as linalg). */
    det(pairs) {
      return new Tokens(call(OPS.det, sized(pairs))).pair();
    },
    detBareiss(pairs) {
      return new Tokens(call(OPS.detBareiss, sized(pairs))).pair();
    },
    detModular(pairs) {
      return new Tokens(call(OPS.detModular, sized(pairs))).pair();
    },
    rank(pairs) {
      return new Tokens(call(OPS.rank, sized(pairs))).int();
    },
    /** Same result as linalg's rrefPairs(pairs, columns): { rows, pivots }. */
    rref(pairs, columns = pairs[0]?.length ?? 0) {
      const m = pairs.length;
      const c = pairs[0]?.length ?? 0;
      const out = [hex(m), hex(c), hex(columns)];
      encodePairs(out, pairs);
      const t = new Tokens(call(OPS.rref, out.join(" ")));
      const np = t.int();
      const pivots = [];
      for (let i = 0; i < np; i++) pivots.push(t.int());
      return { rows: t.pairs(m, c), pivots };
    },
    /** Same result as linalg's solvePairs(a, b): { status, x, basis }. */
    solve(a, b) {
      const m = a.length;
      const n = a[0]?.length ?? 0;
      const k = b[0]?.length ?? 0;
      const out = [hex(m), hex(n), hex(k)];
      encodePairs(out, a);
      encodePairs(out, b);
      const t = new Tokens(call(OPS.solve, out.join(" ")));
      const status = ["unique", "infinite", "none"][t.int()];
      const x = t.int() ? t.pairs(n, k) : null;
      const nb = t.int();
      const basis = [];
      for (let i = 0; i < nb; i++) basis.push(t.pairs(n, 1).map((r) => r[0]));
      return { status, x, basis };
    },
    /** Inverse as rows of pairs, or null when singular. */
    inverse(pairs) {
      const n = pairs.length;
      const t = new Tokens(call(OPS.inverse, sized(pairs)));
      return t.int() ? t.pairs(n, n) : null;
    },
    provablyNonsingular(pairs, tries = 2) {
      return new Tokens(call(OPS.nonsingular, hex(tries) + " " + sized(pairs))).int() === 1;
    },

    /**
     * Certified graph columns, equal to graphics' certifyGraph frame columns.
     *   program   from compileProgram(ast, variable)
     *   params    { name: [n, d] } exact parameter values (names per program.params)
     *   viewport  { xmin, xmax, ymin, ymax: [n, d], width, height }
     *   k         precision bits (graphics' precisionFor(vp, extraBits))
     *   maxDepth  bisection depth (graphics default 8)
     *   from, to  column range (default all)
     * Returns [{ spans: number[], tight, evaluations }] for columns from..to-1.
     */
    graphColumns({ program, params = {}, viewport: vp, k, maxDepth = 8, from = 0, to = vp.width }) {
      const out = [hex(k), hex(maxDepth), hex(vp.width), hex(vp.height), hex(from), hex(to)];
      for (const r of [vp.xmin, vp.xmax, vp.ymin, vp.ymax]) out.push(hex(r[0]), hex(r[1]));
      out.push(hex(program.params.length));
      for (const name of program.params) {
        const v = params[name];
        if (!v) throw new ReferenceError(`no value for ${name}`);
        out.push(hex(v[0]), hex(v[1]));
      }
      out.push(hex(program.consts.length));
      for (const [n, d] of program.consts) out.push(hex(n), hex(d));
      out.push(hex(program.code.length - countArgs(program.code)));
      for (const c of program.code) out.push(hex(c));
      const t = new Tokens(call(OPS.graphColumns, out.join(" ")));
      const n = t.int();
      const cols = new Array(n);
      for (let i = 0; i < n; i++) {
        const tight = t.int() === 1;
        const evaluations = t.int();
        const ns = t.int();
        const spans = new Array(ns);
        for (let j = 0; j < ns; j++) spans[j] = t.int();
        cols[i] = { spans, tight, evaluations };
      }
      return cols;
    },
  };
  backend.version = call(OPS.version, "");
  return backend;
}

// number of inline arguments in a code array (PARAM, CONST, POWI carry one)
function countArgs(code) {
  let n = 0;
  for (let i = 0; i < code.length; i++) {
    if (code[i] === 1 || code[i] === 2 || code[i] === 15) {
      n++;
      i++;
    }
  }
  return n;
}

