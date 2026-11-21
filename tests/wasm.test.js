// Smoke tests of the WASM backend that need no JS reference checkout.
import { expect, test } from "bun:test";
import { loadRustBackend } from "../js/index.js";

const rust = await loadRustBackend(new URL("../js/ratmath.wasm", import.meta.url));

test("version and policy", () => {
  expect(rust.version).toMatch(/^\d+\.\d+\.\d+$/);
  expect(rust.wants("det", 40)).toBe(true);
  expect(rust.wants("det", 3)).toBe(false);
  expect(rust.wants("solve", 1000)).toBe(false);
});

test("determinant, inverse, solve", () => {
  const p = (n, d = 1n) => [n, d];
  expect(rust.det([[p(1n, 2n), p(1n, 3n)], [p(1n, 4n), p(1n, 5n)]])).toEqual([1n, 60n]);
  expect(rust.inverse([[p(2n), p(0n)], [p(0n), p(4n)]])).toEqual([[p(1n, 2n), p(0n)], [p(0n), p(1n, 4n)]]);
  expect(rust.inverse([[p(1n), p(2n)], [p(2n), p(4n)]])).toBeNull();
  const s = rust.solve([[p(1n), p(1n)], [p(1n), p(-1n)]], [[p(3n)], [p(1n)]]);
  expect(s).toEqual({ status: "unique", x: [[p(2n)], [p(1n)]], basis: [] });
});

test("errors come back as exceptions", () => {
  expect(() => rust.det([[[1n, 0n]]])).toThrow(/zero denominator/);
});
