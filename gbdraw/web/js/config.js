export const GBDRAW_WHEEL_NAME = "gbdraw-0.9.0b0-py3-none-any.whl";
export const GBDRAW_WHEEL_CACHE_BUST = "20260222-1600";

export const PYODIDE_INDEX_URL = "./vendor/pyodide/v0.29.0/full/";
export const PYODIDE_WHEELS_BASE_URL = "./vendor/pyodide-wheels/";
export const PYODIDE_REQUIRED_WHEELS = [
  "six-1.17.0-py2.py3-none-any.whl",
  "pytz-2025.2-py2.py3-none-any.whl",
  "python_dateutil-2.9.0.post0-py2.py3-none-any.whl",
  "numpy-2.2.5-cp313-cp313-pyodide_2025_0_wasm32.whl",
  "svgwrite-1.4.3-py3-none-any.whl",
  "fonttools-4.56.0-py3-none-any.whl",
  "biopython-1.85-cp313-cp313-pyodide_2025_0_wasm32.whl",
  "pandas-2.3.2-cp313-cp313-pyodide_2025_0_wasm32.whl",
  "bcbio_gff-0.7.1-py3-none-any.whl"
];

export const DEBUG = false;

export const debugLog = (...args) => {
  if (DEBUG) {
    console.log('[DEBUG]', ...args);
  }
};
