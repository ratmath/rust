// Bit-for-bit cross-checks of the Rust kernels (through WebAssembly) against
// the JS reference implementations in ratmath/linalg and ratmath/graphics.
import { describe, expect, test } from "bun:test";
import { resolve } from "node:path";
import { loadRustBackend, compileProgram } from "../js/index.js";
import { linalgDir, graphicsDir, haveLinalg, haveGraphics } from "./refs.js";

const rust = await loadRustBackend(new URL("../js/ratmath.wasm", import.meta.url));

// deterministic PRNG (mulberry32)
function rng(seed) {
  return () => {
    seed |= 0;
    seed = (seed + 0x6d2b79f5) | 0;
    let t = Math.imul(seed ^ (seed >>> 15), 1 | seed);
    t = (t + Math.imul(t ^ (t >>> 7), 61 | t)) ^ t;
    return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
  };
}

function gcd(a, b) {
  a = a < 0n ? -a : a;
  b = b < 0n ? -b : b;
  while (b) [a, b] = [b, a % b];
  return a;
}
const red = (n, d) => {
  const g = gcd(n, d) || 1n;
  return [n / g, d / g];
};

function randomMatrix(r, rows, cols, { big = false, zeros = 0.2, fractions = 0.4 } = {}) {
  const m = [];
  for (let i = 0; i < rows; i++) {
    const row = [];
    for (let j = 0; j < cols; j++) {
      if (r() < zeros) {
        row.push([0n, 1n]);
        continue;
      }
      let n = BigInt(Math.floor(r() * 41) - 20);
      if (big) n = n * 10n ** BigInt(Math.floor(r() * 30)) + BigInt(Math.floor(r() * 1e9));
      const d = r() < fractions ? BigInt(1 + Math.floor(r() * 12)) : 1n;
      row.push(red(n, d));
    }
    m.push(row);
  }
  return m;
}

function withDependentRow(r, m) {
  // make the last row a combination of the first two (forces rank deficiency)
  if (m.length < 3) return m;
  const a = BigInt(Math.floor(r() * 5) - 2);
  const out = m.map((row) => row.map((p) => [...p]));
  out[out.length - 1] = m[0].map(([n1, d1], j) => {
    const [n2, d2] = m[1][j];
    return red(a * n1 * d2 + n2 * d1, d1 * d2);
  });
  return out;
}

const shapes = [];
{
  const r = rng(12345);
  for (let t = 0; t < 120; t++) {
    const rows = 1 + Math.floor(r() * 9);
    const cols = t % 3 === 0 ? rows : 1 + Math.floor(r() * 9);
    let m = randomMatrix(r, rows, cols, { big: t % 5 === 0, zeros: t % 7 === 0 ? 0.6 : 0.15 });
    if (t % 4 === 0) m = withDependentRow(r, m);
    shapes.push(m);
  }
}

describe.skipIf(!haveLinalg)("linalg cross-check", async () => {
  if (!haveLinalg) return; // bun still runs a skipped describe body
  const X = await import(resolve(linalgDir, "src/exact.js"));
  const squares = shapes.filter((m) => m.length === m[0].length);

  test("determinant (Bareiss, modular, dispatching) equals linalg", () => {
    for (const m of squares) {
      const ref = X.run(X.determinantBareiss(m));
      expect(X.run(X.determinantModular(m))).toEqual(ref);
      expect(rust.detBareiss(m)).toEqual(ref);
      expect(rust.detModular(m)).toEqual(ref);
      expect(rust.det(m)).toEqual(ref);
    }
  });

  test("rref equals rrefPairs exactly (rows and pivots)", () => {
    for (const m of shapes) {
      for (const columns of [m[0].length, Math.max(1, m[0].length - 1)]) {
        const ref = X.run(X.rrefPairs(m, columns));
        expect(rust.rref(m, columns)).toEqual(ref);
      }
    }
  });

  test("rank", () => {
    for (const m of shapes) expect(rust.rank(m)).toBe(X.run(X.rrefPairs(m, m[0].length)).pivots.length);
  });

  test("solve equals solvePairs (unique, infinite and none)", () => {
    const r = rng(7);
    const seen = new Set();
    for (const a of shapes) {
      const b = randomMatrix(r, a.length, 1 + Math.floor(r() * 2));
      const ref = X.run(X.solvePairs(a, b));
      seen.add(ref.status);
      expect(rust.solve(a, b)).toEqual(ref);
    }
    expect([...seen].sort()).toEqual(["infinite", "none", "unique"]);
  });

  test("inverse equals inversePairs (including singular → null)", () => {
    let singular = 0;
    for (const m of squares) {
      const ref = X.run(X.inversePairs(m));
      if (ref === null) singular++;
      expect(rust.inverse(m)).toEqual(ref);
    }
    expect(singular).toBeGreaterThan(0);
  });

  test("provablyNonsingular is never wrong when true", () => {
    for (const m of squares) {
      const nz = X.run(X.determinantBareiss(m))[0] !== 0n;
      if (rust.provablyNonsingular(m)) expect(nz).toBe(true);
      if (!nz) expect(rust.provablyNonsingular(m)).toBe(false);
    }
  });

  test("larger matrices: 40×40 determinant and 25×25 inverse", () => {
    const r = rng(99);
    const m = randomMatrix(r, 40, 40, { zeros: 0 });
    expect(rust.det(m)).toEqual(X.run(X.determinantModular(m)));
    const s = randomMatrix(r, 25, 25, { zeros: 0 });
    expect(rust.inverse(s)).toEqual(X.run(X.inversePairs(s)));
  });
});

describe.skipIf(!haveGraphics)("graphics cross-check", async () => {
  if (!haveGraphics) return;
  const { certifyGraph } = await import(resolve(graphicsDir, "src/plot.js"));
  const { parse, paramsToRational } = await import(resolve(graphicsDir, "src/expr.js"));
  const { viewport, precisionFor } = await import(resolve(graphicsDir, "src/viewport.js"));

  const drive = (gen) => {
    for (;;) {
      const { value, done } = gen.next();
      if (done) return value;
    }
  };

  const cases = [
    ["sin(x)", {}, { x: [-7, 7], y: [-2, 2] }],
    ["x^2 - 2", {}, { x: ["-3/2", "5/2"], y: [-3, 5] }],
    ["1/x", {}, { x: [-2, 2], y: [-5, 5] }],
    ["tan(x)", {}, { x: [-4, 4], y: [-6, 6] }],
    ["exp(x)", {}, { x: [-3, 3], y: [-1, 20] }],
    ["ln(x)", {}, { x: ["-1/3", 5], y: [-4, 2] }],
    ["sqrt(x) - 1", {}, { x: [-1, 4], y: [-2, 2] }],
    ["a*sin(b*x) + c", { a: "3/2", b: "2.5", c: "-1/7" }, { x: [-5, 5], y: [-3, 3] }],
    ["x^x", {}, { x: ["1/100", 3], y: [0, 4] }],
    ["abs(x) - 1", {}, { x: [-2, 2], y: [-2, 2] }],
    ["min(x, 1) + max(x, -1)/3", {}, { x: [-3, 3], y: [-3, 3] }],
    ["e^x - pi", {}, { x: [-2, 2], y: [-4, 4] }],
    ["cos(3x)/x", {}, { x: [-3, 3], y: [-10, 10] }],
    ["x^-2", {}, { x: [-2, 2], y: [-1, 10] }],
    ["sin(1/x)", {}, { x: ["-1/2", "1/2"], y: [-2, 2] }],
    ["x^3 - x", {}, { x: [-1000000, 1000000], y: [-1e15, 1e15] }],
    ["cos(x)", {}, { x: [99990, 100010], y: [-2, 2] }],
    ["x^(1/2) + x^(2/3)", {}, { x: [-1, 4], y: [-1, 5] }],
    ["x^(-1/2) - (x+1)^(3/2 - 1)", {}, { x: [-1, 3], y: [-3, 5] }],
    ["(x^2)^(1/3) - 2^x", {}, { x: [-2, 2], y: [-4, 4] }],
  ];

  for (const [f, params, box] of cases) {
    test(`certified columns equal graphics for ${f}`, () => {
      const vp = viewport({ ...box, width: 120, height: 80 });
      const ast = parse(f);
      const ref = drive(certifyGraph(ast, params, vp, { asymptotes: false })).value;
      const program = compileProgram(ast, "x");
      const rp = paramsToRational(params);
      const p = Object.fromEntries(Object.entries(rp).map(([k, v]) => [k, [v.numerator, v.denominator]]));
      const pair = (q) => [q.numerator, q.denominator];
      const cols = rust.graphColumns({
        program,
        params: p,
        viewport: { xmin: pair(vp.xmin), xmax: pair(vp.xmax), ymin: pair(vp.ymin), ymax: pair(vp.ymax), width: vp.width, height: vp.height },
        k: precisionFor(vp, 12),
        maxDepth: 8,
      });
      expect(cols.map((c) => c.spans)).toEqual(ref.columns);
      expect(cols.every((c) => c.tight)).toBe(ref.tight);
      expect(cols.reduce((s, c) => s + c.evaluations, 0)).toBe(ref.evaluations);
    });
  }
});

describe.skipIf(!haveGraphics)("graphics fuzz", async () => {
  if (!haveGraphics) return;
  const { certifyGraph } = await import(resolve(graphicsDir, "src/plot.js"));
  const { parse } = await import(resolve(graphicsDir, "src/expr.js"));
  const { viewport, precisionFor } = await import(resolve(graphicsDir, "src/viewport.js"));
  const drive = (gen) => {
    for (;;) {
      const { value, done } = gen.next();
      if (done) return value;
    }
  };
  const r = rng(2026);
  const pick = (a) => a[Math.floor(r() * a.length)];
  const atom = () => pick(["x", "x", "2", "1/3", "pi", "e", "0.7", "x"]);
  function expr(depth) {
    if (depth <= 0 || r() < 0.25) return atom();
    switch (Math.floor(r() * 4)) {
      case 0: return `(${expr(depth - 1)} ${pick(["+", "-", "*", "/"])} ${expr(depth - 1)})`;
      case 1: return `${pick(["sin", "cos", "tan", "exp", "ln", "sqrt", "abs"])}(${expr(depth - 1)})`;
      case 2: return `(${expr(depth - 1)})^${pick(["2", "3", "-1", "(1/2)", "(-3/2)", "(1/2+1/3)", "x"])}`;
      default: return `${pick(["min", "max"])}(${expr(depth - 1)}, ${expr(depth - 1)})`;
    }
  }
  const boxes = [
    { x: [-6, 6], y: [-4, 4] },
    { x: ["-1/1000", "1/1000"], y: ["-1/1000", "1/1000"] },
    { x: ["314159/100000", "314160/100000"], y: [-1, 1] },
    { x: [-50, 50], y: [-1000, 1000] },
  ];
  for (let t = 0; t < 40; t++) {
    const f = expr(3);
    const box = boxes[t % boxes.length];
    test(`fuzz ${t}: ${f}`, () => {
      const vp = viewport({ ...box, width: 48, height: 40 });
      const ast = parse(f);
      let ref;
      try {
        ref = drive(certifyGraph(ast, {}, vp, { asymptotes: false })).value;
      } catch (e) {
        return; // the JS kernel itself throws here (e.g. a crash path); nothing to compare
      }
      const pair = (q) => [q.numerator, q.denominator];
      const cols = rust.graphColumns({
        program: compileProgram(ast, "x"),
        viewport: { xmin: pair(vp.xmin), xmax: pair(vp.xmax), ymin: pair(vp.ymin), ymax: pair(vp.ymax), width: vp.width, height: vp.height },
        k: precisionFor(vp, 12),
      });
      expect(cols.map((c) => c.spans)).toEqual(ref.columns);
      expect(cols.reduce((s, c) => s + c.evaluations, 0)).toBe(ref.evaluations);
    });
  }
});
