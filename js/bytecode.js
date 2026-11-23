// Compile a ratmath/graphics expression AST (src/expr.js: num, const, var,
// op, call nodes) into the stack program the Rust `graph` module runs. The
// translation mirrors graphics' `compileInterval` node for node, so the Rust
// enclosure equals the JS one.

export const OP = {
  X: 0, PARAM: 1, CONST: 2, PI: 3, E: 4,
  ADD: 10, SUB: 11, MUL: 12, DIV: 13, NEG: 14, POWI: 15, POW: 16, POWPOS: 17, RECIP: 18,
  SIN: 20, COS: 21, TAN: 22, EXP: 23, LN: 24, SQRT: 25, ABS: 26, MIN: 27, MAX: 28,
};

const UNARY = { sin: OP.SIN, cos: OP.COS, tan: OP.TAN, exp: OP.EXP, ln: OP.LN, log: OP.LN, sqrt: OP.SQRT, abs: OP.ABS };
const BINARY = { "+": OP.ADD, "-": OP.SUB, "*": OP.MUL, "/": OP.DIV };

/**
 * The exact value of a constant, rational-only subtree, or null (graphics
 * `constValue`): `1/2`, `-3`, `(1/2 - 1)`. Uses the AST's Core Rationals.
 */
export function constValue(ast) {
  switch (ast.type) {
    case "num":
      return ast.value;
    case "op": {
      const vals = ast.args.map(constValue);
      if (vals.some((x) => x === null)) return null;
      const [a, b] = vals;
      switch (ast.op) {
        case "+": return a.add(b);
        case "-": return a.subtract(b);
        case "*": return a.multiply(b);
        case "/": return b.numerator === 0n ? null : a.divide(b);
        case "neg": return a.negate();
        case "^": {
          if (b.denominator !== 1n || b.numerator > 4096n || b.numerator < -4096n) return null;
          if (a.numerator === 0n && b.numerator <= 0n) return null;
          return a.pow(Number(b.numerator));
        }
      }
    }
  }
  return null;
}

/** The exponent as a JS integer when it is a constant integer (graphics `integerExponent`). */
function integerExponent(ast) {
  const c = constValue(ast);
  if (c && c.denominator === 1n && c.numerator >= -4096n && c.numerator <= 4096n) return Number(c.numerator);
  return null;
}

/**
 * Returns { code: number[] (opcodes with inline arguments), consts: [n, d][],
 * params: string[] (names, in PARAM index order) }.
 */
export function compileProgram(ast, variable = "x") {
  const code = [];
  const consts = [];
  const params = [];
  const go = (n) => {
    switch (n.type) {
      case "num":
        consts.push([n.value.numerator, n.value.denominator]);
        code.push(OP.CONST, consts.length - 1);
        return;
      case "const":
        code.push(n.name === "pi" ? OP.PI : OP.E);
        return;
      case "var": {
        if (n.name === variable) {
          code.push(OP.X);
          return;
        }
        let i = params.indexOf(n.name);
        if (i < 0) i = params.push(n.name) - 1;
        code.push(OP.PARAM, i);
        return;
      }
      case "op": {
        if (n.op === "neg") {
          go(n.args[0]);
          code.push(OP.NEG);
          return;
        }
        if (n.op === "^") {
          const k = integerExponent(n.args[1]);
          go(n.args[0]);
          if (k !== null) {
            code.push(OP.POWI, k);
            return;
          }
          const cv = constValue(n.args[1]);
          go(n.args[1]);
          if (cv === null) code.push(OP.POW); // exp(b · ln a)
          else if (cv.numerator > 0n) code.push(OP.POWPOS);
          else code.push(OP.NEG, OP.POWPOS, OP.RECIP); // 1 / a^(-c)
          return;
        }
        const o = BINARY[n.op];
        if (o === undefined) break;
        go(n.args[0]);
        go(n.args[1]);
        code.push(o);
        return;
      }
      case "call": {
        if (UNARY[n.name] !== undefined) {
          go(n.args[0]);
          code.push(UNARY[n.name]);
          return;
        }
        if (n.name === "min" || n.name === "max") {
          go(n.args[0]);
          go(n.args[1]);
          code.push(n.name === "min" ? OP.MIN : OP.MAX);
          return;
        }
        break;
      }
    }
    throw new Error(`cannot compile ${n.type} ${n.op ?? n.name ?? ""}`);
  };
  go(ast);
  return { code, consts, params };
}
