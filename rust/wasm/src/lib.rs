//! The drawing engine for the diagram Worker: a plain C ABI, so the module needs
//! no Emscripten build matched to Pyodide. `gbdraw/web/js/workers/gbdraw-core.js`
//! is its only caller.

use std::cell::RefCell;

thread_local! {
    static RESULT: RefCell<Vec<u8>> = const { RefCell::new(Vec::new()) };
}

#[no_mangle]
pub extern "C" fn gbdraw_contract_version() -> u32 {
    gbdraw_engine::CONTRACT_VERSION
}

#[no_mangle]
pub extern "C" fn gbdraw_alloc(length: usize) -> *mut u8 {
    Box::into_raw(vec![0u8; length].into_boxed_slice()).cast()
}

/// # Safety
/// `pointer` and `length` come from one `gbdraw_alloc` call.
#[no_mangle]
pub unsafe extern "C" fn gbdraw_free(pointer: *mut u8, length: usize) {
    drop(Box::from_raw(std::ptr::slice_from_raw_parts_mut(
        pointer, length,
    )));
}

/// 0 on success, with the text at `gbdraw_result_pointer`; 1 with the reason there.
///
/// # Safety
/// The pointers and lengths name valid bytes.
#[no_mangle]
pub unsafe extern "C" fn gbdraw_pairwise_ribbons(
    params: *const u8,
    params_length: usize,
    columns: *const u8,
    columns_length: usize,
) -> u32 {
    let params = std::slice::from_raw_parts(params, params_length);
    let columns = std::slice::from_raw_parts(columns, columns_length);
    let drawn = std::str::from_utf8(params)
        .map_err(|error| gbdraw_engine::Unsupported(error.to_string()))
        .and_then(|params| gbdraw_engine::ribbons::render(params, columns));
    let (status, text) = match drawn {
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
