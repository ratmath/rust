# Wire format

Every operation takes one text input and returns one text output through
`rm_call(op, input, len)` (see `include/ratmath.h`). Tokens are separated by
single spaces. **Every integer is hexadecimal**, with an optional leading `-`:
exactly what JavaScript's `BigInt#toString(16)` / `Number#toString(16)` print and
what `BigInt("0x" + s)` reads back (both linear time, unlike decimal).

A rational is two integers `n d` (d ≠ 0; inputs are reduced on read, outputs
are reduced with d > 0). A matrix is its entries row by row.

| op | code | input | output |
|---|---|---|---|
| version | 0 | (empty) | version string, e.g. `0.1.0` |
| det | 1 | `rows cols` entries | `n d` (Bareiss for n ≤ 3, multi-modular above) |
| detBareiss | 2 | same | `n d` |
| detModular | 3 | same | `n d` |
| rank | 4 | `rows cols` entries | rank |
| rref | 5 | `rows cols columns` entries | `npivots` pivot columns… then the rows×cols entries |
| solve | 6 | `m n k` A (m×n) B (m×k) | `status` (0 unique, 1 infinite, 2 none), `hasX`, X (n×k if hasX), `nbasis`, basis vectors (n entries each) |
| inverse | 7 | `n n` entries | `0` (singular) or `1` followed by the n×n entries |
| nonsingular | 8 | `tries` then `n n` entries | `1` if provably nonsingular mod one of `tries` primes, else `0` |
| graphColumns | 20 | see below | see below |

`rref`, `solve` and `inverse` return exactly what ratmath/linalg's
`rrefPairs`, `solvePairs` and `inversePairs` return (same pivots, same row
order, same basis vectors).

## graphColumns

Input:

```
k maxDepth width height from to
xmin_n xmin_d xmax_n xmax_d ymin_n ymin_d ymax_n ymax_d
nparams (n d)…            parameter values, in program.params order
nconsts (n d)…            literals referenced by CONST
ncode  instructions…      ncode counts instructions, not inline arguments
```

Instructions (stack machine; binary operators pop b, then a, and push a∘b):

| code | instruction | | code | instruction |
|---|---|---|---|---|
| 0 | X (the plotted variable) | | 15 `k` | integer power a^k |
| 1 `i` | parameter i | | 16 | a^b = exp(b·ln a) |
| 2 `i` | constant i | | 20–22 | sin, cos, tan |
| 3 | π | | 23, 24 | exp, ln |
| 4 | e | | 25, 26 | sqrt, abs |
| 17 | a^b for constant b > 0 (graphics `powPos`) | | 18 | 1/a |
| 10–13 | + − × ÷ | | 27, 28 | min, max |
| 14 | negate | | | |

`js/bytecode.js` compiles a ratmath/graphics AST to this form.

Output: `ncolumns`, then for each column `tight evaluations nspans spans…`,
where spans are the merged inclusive pixel-row ranges top, bottom, top, bottom…
— the same as a column of graphics' `certifyGraph` frame.
