//! C ABI, the one entry point for both WASM (JS) and native (Swift on iOS,
//! Kotlin/JNI on Android later). Input and output are wire-format text
//! (WIRE.md); the declarations are in `include/ratmath.h`.
//!
//! Protocol: `rm_alloc(len)` a buffer, write the input into it, call
//! `rm_call(op, ptr, len)`, `rm_free(ptr, len)`; then read `rm_out_len()`
//! bytes at `rm_out_ptr()`. The output buffer belongs to the library and stays
//! valid until the next `rm_call` on the same thread. `rm_call` returns 0 on
//! success, 1 when the output is an error message.

use std::cell::RefCell;

thread_local! {
    static OUT: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
}

#[unsafe(no_mangle)]
pub extern "C" fn rm_alloc(len: usize) -> *mut u8 {
    let mut v = Vec::<u8>::with_capacity(len.max(1));
    let p = v.as_mut_ptr();
    std::mem::forget(v);
    p
}

/// # Safety
/// `ptr` must come from `rm_alloc(len)` with the same `len`, freed once.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rm_free(ptr: *mut u8, len: usize) {
    if !ptr.is_null() {
        drop(unsafe { Vec::from_raw_parts(ptr, 0, len.max(1)) });
    }
}

/// # Safety
/// `ptr` must point to `len` readable bytes.
#[unsafe(no_mangle)]
pub unsafe extern "C" fn rm_call(op: u32, ptr: *const u8, len: usize) -> i32 {
    let input = if len == 0 {
        &[][..]
    } else {
        unsafe { std::slice::from_raw_parts(ptr, len) }
    };
    let result = match std::str::from_utf8(input) {
        Ok(s) => std::panic::catch_unwind(|| crate::wire::call(op, s))
            .unwrap_or_else(|_| Err("internal error (panic)".to_string())),
        Err(_) => Err("input is not UTF-8".to_string()),
    };
    let (code, text) = match result {
        Ok(t) => (0, t),
        Err(e) => (1, e),
    };
    OUT.with(|o| {
        let mut o = o.borrow_mut();
        o.clear();
        o.extend_from_slice(text.as_bytes());
    });
    code
}

#[unsafe(no_mangle)]
pub extern "C" fn rm_out_ptr() -> *const u8 {
    OUT.with(|o| o.borrow().as_ptr())
}

#[unsafe(no_mangle)]
pub extern "C" fn rm_out_len() -> usize {
    OUT.with(|o| o.borrow().len())
}
