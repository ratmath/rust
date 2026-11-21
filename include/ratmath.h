/*
 * ratmath.h — C ABI of the RatMath Rust kernels (github.com/ratmath/rust).
 *
 * One entry point for every caller: Swift (iOS/macOS, via a module map or a
 * bridging header), Kotlin/JNI (Android, later), and WebAssembly (js/index.js
 * does the same calls). Inputs and outputs are wire-format text; see WIRE.md.
 *
 *   uint8_t *buf = rm_alloc(len);  memcpy(buf, input, len);
 *   int32_t status = rm_call(RM_OP_DET, buf, len);
 *   rm_free(buf, len);
 *   // rm_out_len() bytes at rm_out_ptr(): the result (status 0) or an error message (status 1)
 *
 * The output buffer is owned by the library and stays valid until the next
 * rm_call on the same thread. Calls on different threads are independent.
 */
#ifndef RATMATH_H
#define RATMATH_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

enum {
  RM_OP_VERSION = 0,
  RM_OP_DET = 1,           /* Bareiss for n <= 3, multi-modular above */
  RM_OP_DET_BAREISS = 2,
  RM_OP_DET_MODULAR = 3,
  RM_OP_RANK = 4,
  RM_OP_RREF = 5,
  RM_OP_SOLVE = 6,
  RM_OP_INVERSE = 7,
  RM_OP_NONSINGULAR = 8,
  RM_OP_GRAPH_COLUMNS = 20
};

uint8_t *rm_alloc(size_t len);
void rm_free(uint8_t *ptr, size_t len);
int32_t rm_call(uint32_t op, const uint8_t *input, size_t len);
const uint8_t *rm_out_ptr(void);
size_t rm_out_len(void);

#ifdef __cplusplus
}
#endif

#endif /* RATMATH_H */
