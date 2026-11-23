// Timing of the JS reference vs the Rust kernels through WebAssembly
// (including the cost of crossing: encoding BigInts to text and back).
//   bun bench/bench.js            (needs ../linalg and ../graphics checkouts)
import { resolve } from "node:path";
import { loadRustBackend, compileProgram } from "../js/index.js";
import { linalgDir, graphicsDir, haveLinalg, haveGraphics } from "../tests/refs.js";

const rust = await loadRustBackend(new URL("../js/ratmath.wasm", import.meta.url));
// Each case's wire input is saved so `cargo run --release --example bench`
// can time the same work natively (what an iOS app calling the C ABI gets).
const { writeFileSync, mkdirSync } = await import("node:fs");
const dump = [];
let label = null;
rust._onCall = (op, text) => { if (label && !dump.find((d) => d.label === label)) dump.push({ label, op, text }); };

function time(fn, minMs = 300) {
  fn();
  let n = 0;
  const t0 = performance.now();
  let t;
  do {
    fn();
    n++;
    t = performance.now() - t0;
  } while (t < minMs && n < 1000);
  return t / n;
}
const fmt = (ms) => (ms < 1 ? `${(ms * 1000).toFixed(0)} µs` : `${ms.toFixed(1)} ms`);

let seed = 42;
const rnd = () => ((seed = (seed * 1103515245 + 12345) % 2147483648) / 2147483648);
function matrix(n, m = n, den = 6) {
  return Array.from({ length: n }, () =>
    Array.from({ length: m }, () => {
      const a = BigInt(Math.floor(rnd() * 21) - 10);
      const d = BigInt(1 + Math.floor(rnd() * den));
      let x = a, y = d;
      while (y) [x, y] = [y, x % y];
      const g = x < 0n ? -x : x || 1n;
      return [a / g, d / g];
    }),
  );
}

const rows = [];
if (haveLinalg) {
  const X = await import(resolve(linalgDir, "src/exact.js"));
  for (const n of [4, 8, 16, 32, 64, 100]) {
    const m = matrix(n);
    const js = time(() => X.run(n <= 3 ? X.determinantBareiss(m) : X.determinantModular(m)));
    label = `det ${n}`;
    const rs = time(() => rust.det(m));
    rows.push(["det", n, js, rs]);
  }
  for (const n of [4, 8, 16, 32, 48]) {
    const m = matrix(n);
    const b = matrix(n, 1);
    label = `solve ${n}`;
    rows.push(["solve", n, time(() => X.run(X.solvePairs(m, b))), time(() => rust.solve(m, b))]);
    label = `inverse ${n}`;
    rows.push(["inverse", n, time(() => X.run(X.inversePairs(m))), time(() => rust.inverse(m))]);
    label = `rref ${n}`;
    rows.push(["rref", n, time(() => X.run(X.rrefPairs(m, n))), time(() => rust.rref(m, n))]);
  }
}
if (haveGraphics) {
  const { certifyGraph } = await import(resolve(graphicsDir, "src/plot.js"));
  const { parse } = await import(resolve(graphicsDir, "src/expr.js"));
  const { viewport, precisionFor } = await import(resolve(graphicsDir, "src/viewport.js"));
  const drive = (g) => { for (;;) { const s = g.next(); if (s.done) return s.value; } };
  const pair = (q) => [q.numerator, q.denominator];
  for (const f of ["x^2 - 2", "sin(x)", "a*sin(3x)/x + exp(-x^2)"]) {
    const vp = viewport({ x: [-6, 6], y: [-3, 3], width: 300, height: 200 });
    const ast = parse(f);
    const program = compileProgram(ast);
    const params = { a: [3n, 2n] };
    const js = time(() => drive(certifyGraph(ast, { a: "3/2" }, vp, { asymptotes: false })), 1000);
    label = `graph ${f}`;
    const rs = time(() => rust.graphColumns({
      program, params, k: precisionFor(vp, 12),
      viewport: { xmin: pair(vp.xmin), xmax: pair(vp.xmax), ymin: pair(vp.ymin), ymax: pair(vp.ymax), width: 300, height: 200 },
    }), 1000);
    rows.push([`graph ${f}`, "300×200", js, rs]);
  }
}

label = null;
mkdirSync(new URL("./cases/", import.meta.url), { recursive: true });
writeFileSync(new URL("./cases/cases.txt", import.meta.url), dump.map((d) => `${d.label}\t${d.op}\t${d.text}`).join("\n") + "\n");
console.log(`ratmath/rust ${rust.version}\n`);
console.log("| operation | size | JS | Rust (WASM) | speed-up |\n|---|---|---|---|---|");
for (const [op, n, js, rs] of rows) console.log(`| ${op} | ${n} | ${fmt(js)} | ${fmt(rs)} | ${(js / rs).toFixed(1)}× |`);
