//! Native timings of the same work `bun bench/bench.js` measures, read from
//! the wire inputs it saves in bench/cases/cases.txt (run that first):
//!     bun bench/bench.js && cargo run --release --example bench
//! Includes parsing the input and writing the output, as a Swift caller pays.
use std::time::Instant;

fn main() {
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/bench/cases/cases.txt");
    let text =
        std::fs::read_to_string(path).expect("run `bun bench/bench.js` first to save the cases");
    println!("| operation | native |\n|---|---|");
    for line in text.lines().filter(|l| !l.is_empty()) {
        let mut parts = line.splitn(3, '\t');
        let (label, op, input) = (
            parts.next().unwrap(),
            parts.next().unwrap(),
            parts.next().unwrap_or(""),
        );
        let op: u32 = op.parse().unwrap();
        ratmath::wire::call(op, input).unwrap();
        let t0 = Instant::now();
        let mut n = 0;
        while (t0.elapsed().as_secs_f64() < 0.3 && n < 1000) || n == 0 {
            ratmath::wire::call(op, input).unwrap();
            n += 1;
        }
        let ms = t0.elapsed().as_secs_f64() * 1000.0 / n as f64;
        let shown = if ms < 1.0 {
            format!("{:.0} µs", ms * 1000.0)
        } else {
            format!("{ms:.1} ms")
        };
        println!("| {label} | {shown} |");
    }
}
