// Encoders and decoders for the wire format (see ../WIRE.md): space-separated
// tokens, every integer in hexadecimal with an optional leading "-".

export const hex = (v) => v.toString(16);

export function big(tok) {
  return tok[0] === "-" ? -BigInt("0x" + tok.slice(1)) : BigInt("0x" + tok);
}

export const int = (tok) => Number(big(tok));

export function encodePairs(out, rows) {
  for (const row of rows) for (const [n, d] of row) out.push(hex(n), hex(d));
}

export class Tokens {
  constructor(text) {
    this.t = text.length ? text.split(" ") : [];
    this.i = 0;
  }
  big() {
    return big(this.t[this.i++]);
  }
  int() {
    return int(this.t[this.i++]);
  }
  pair() {
    const n = this.big();
    return [n, this.big()];
  }
  pairs(rows, cols) {
    const out = new Array(rows);
    for (let r = 0; r < rows; r++) {
      const row = new Array(cols);
      for (let c = 0; c < cols; c++) row[c] = this.pair();
      out[r] = row;
    }
    return out;
  }
}
