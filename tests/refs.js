// Locate the JS reference libraries for cross-checking. Defaults to sibling
// checkouts (../linalg, ../graphics); override with RATMATH_LINALG /
// RATMATH_GRAPHICS. Each needs its dependencies installed (bun install).
import { existsSync } from "node:fs";
import { resolve } from "node:path";

const root = resolve(import.meta.dir, "..");
export const linalgDir = resolve(process.env.RATMATH_LINALG ?? resolve(root, "../linalg"));
export const graphicsDir = resolve(process.env.RATMATH_GRAPHICS ?? resolve(root, "../graphics"));
export const haveLinalg = existsSync(resolve(linalgDir, "src/exact.js"));
export const haveGraphics = existsSync(resolve(graphicsDir, "src/plot.js"));
