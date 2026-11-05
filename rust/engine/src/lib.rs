//! gbdraw's drawing engine: the SVG text of diagram groups, byte for byte as
//! the Python drawers write it.

pub mod ribbons;
pub mod svg;

use std::panic::{catch_unwind, AssertUnwindSafe};

/// What this engine draws and how it is called. Raise it, with `ENGINE_CONTRACT`
/// in `gbdraw/render/native.py`, whenever the drawn bytes or the call change:
/// Python then leaves an older installed or cached engine unused.
pub const CONTRACT_VERSION: u32 = 1;

/// Why a group cannot be drawn here. The caller draws it in Python, which
/// reports any input error in gbdraw's own terms.
#[derive(Debug, PartialEq, Eq)]
pub struct Unsupported(pub String);

impl std::fmt::Display for Unsupported {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}

/// `draw`, with a panic reported as a refusal: a bug here leaves the group to
/// Python instead of stopping the render. Wasm aborts instead of unwinding, and
/// its caller treats the trap as the same refusal.
pub fn guarded<T>(draw: impl FnOnce() -> Result<T, Unsupported>) -> Result<T, Unsupported> {
    catch_unwind(AssertUnwindSafe(draw))
        .unwrap_or_else(|_| Err(Unsupported("the engine panicked".to_string())))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_panic_is_a_refusal() {
        assert_eq!(guarded(|| Ok(1)), Ok(1));
        assert_eq!(
            guarded(|| -> Result<(), Unsupported> { panic!("a bug") }),
            Err(Unsupported("the engine panicked".to_string()))
        );
    }
}
