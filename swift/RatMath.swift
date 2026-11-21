// A thin Swift wrapper over the C ABI (sketch; add RatMathCore.xcframework to
// the target). Values cross as wire-format text (WIRE.md); a Swift-side
// big-integer type can encode/decode hex the same way js/wire.js does.
import RatMathCore

public enum RatMathError: Error { case kernel(String) }

public enum RatMath {
    /// Run one kernel operation on wire-format input.
    public static func call(_ op: UInt32, _ input: String) throws -> String {
        let bytes = Array(input.utf8)
        let status = bytes.withUnsafeBufferPointer { buf -> Int32 in
            rm_call(op, buf.baseAddress, buf.count)
        }
        let out = String(decoding: UnsafeBufferPointer(start: rm_out_ptr(), count: rm_out_len()), as: UTF8.self)
        if status != 0 { throw RatMathError.kernel(out) }
        return out
    }

    public static var version: String { (try? call(UInt32(RM_OP_VERSION), "")) ?? "?" }

    /// Determinant of a square matrix of (numerator, denominator) hex strings.
    public static func det(_ m: [[(String, String)]]) throws -> (String, String) {
        let n = m.count
        var toks = [String(n, radix: 16), String(n, radix: 16)]
        for row in m { for (a, b) in row { toks.append(a); toks.append(b) } }
        let out = try call(UInt32(RM_OP_DET), toks.joined(separator: " ")).split(separator: " ")
        return (String(out[0]), String(out[1]))
    }
}
