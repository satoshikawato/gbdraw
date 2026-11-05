//! The drawing engine for the diagram Worker: a plain C ABI, so the module needs
//! no Emscripten build matched to Pyodide. `gbdraw/web/js/workers/gbdraw-core.js`
//! is its only caller.

use std::cell::RefCell;

thread_local! {
    static RESULT: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
}

#[no_mangle]
pub extern "C" fn gbdraw_alloc(length: usize) -> *mut u8 {
    let mut buffer = std::mem::ManuallyDrop::new(Vec::<u8>::with_capacity(length));
    buffer.as_mut_ptr()
}

/// # Safety
/// `pointer` and `length` come from one `gbdraw_alloc` call.
#[no_mangle]
pub unsafe extern "C" fn gbdraw_free(pointer: *mut u8, length: usize) {
    drop(Vec::from_raw_parts(pointer, 0, length));
}

/// 0 on success, with the text at `gbdraw_result_pointer`; 1 with the reason there.
///
/// # Safety
/// The pointers and lengths name valid bytes; `params` is UTF-8.
#[no_mangle]
pub unsafe extern "C" fn gbdraw_pairwise_ribbons(params: *const u8, params_length: usize, columns: *const u8, columns_length: usize) -> u32 {
    let params = std::str::from_utf8_unchecked(std::slice::from_raw_parts(params, params_length));
    let columns = std::slice::from_raw_parts(columns, columns_length);
    let (status, text) = match gbdraw_engine::ribbons::render(params, columns) {
        Ok(text) => (0, text),
        Err(reason) => (1, reason.0),
    };
    RESULT.with(|result| *result.borrow_mut() = text.into_bytes());
    status
}

#[no_mangle]
pub extern "C" fn gbdraw_result_pointer() -> *const u8 {
    RESULT.with(|result| result.borrow().as_ptr())
}

#[no_mangle]
pub extern "C" fn gbdraw_result_length() -> usize {
    RESULT.with(|result| result.borrow().len())
}
