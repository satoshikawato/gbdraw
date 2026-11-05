//! `gbdraw_core`: the drawing engine for native Python.

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyBytes;

/// The children text of a pairwise comparison group (`gbdraw_engine::ribbons`), as UTF-8 bytes.
#[pyfunction]
fn pairwise_ribbons<'py>(py: Python<'py>, params: &str, columns: &[u8]) -> PyResult<Bound<'py, PyBytes>> {
    py.detach(|| gbdraw_engine::ribbons::render(params, columns))
        .map(|text| PyBytes::new(py, text.as_bytes()))
        .map_err(|reason| PyValueError::new_err(reason.to_string()))
}

#[pymodule]
fn gbdraw_core(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(pairwise_ribbons, module)?)
}
