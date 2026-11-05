//! gbdraw's drawing engine: the SVG text of diagram groups, byte for byte as
//! the Python drawers write it.

pub mod ribbons;
pub mod svg;

/// Why a group cannot be drawn here. The caller draws it in Python, which
/// reports any input error in gbdraw's own terms.
#[derive(Debug, PartialEq, Eq)]
pub struct Unsupported(pub String);

impl std::fmt::Display for Unsupported {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(&self.0)
    }
}
