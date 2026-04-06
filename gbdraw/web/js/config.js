export const GBDRAW_WHEEL_NAME = "gbdraw-0.9.0-py3-none-any.whl";
export const GBDRAW_WHEEL_CACHE_BUST = "a464c04d3161";
export const PYODIDE_INDEX_URL = "./vendor/pyodide/v0.29.0/full/";
export const PYODIDE_LOCAL_WHEELS = [
  "./vendor/pyodide-wheels/six-1.17.0-py2.py3-none-any.whl",
  "./vendor/pyodide-wheels/python_dateutil-2.9.0.post0-py2.py3-none-any.whl",
  "./vendor/pyodide-wheels/pytz-2025.2-py2.py3-none-any.whl",
  "./vendor/pyodide-wheels/numpy-2.2.5-cp313-cp313-pyodide_2025_0_wasm32.whl",
  "./vendor/pyodide-wheels/biopython-1.85-cp313-cp313-pyodide_2025_0_wasm32.whl",
  "./vendor/pyodide-wheels/fonttools-4.56.0-py3-none-any.whl",
  "./vendor/pyodide-wheels/svgwrite-1.4.3-py3-none-any.whl",
  "./vendor/pyodide-wheels/bcbio_gff-0.7.1-py3-none-any.whl",
  "./vendor/pyodide-wheels/pandas-2.3.2-cp313-cp313-pyodide_2025_0_wasm32.whl"
];
export const WASI_SHIM_URL = "./vendor/browser_wasi_shim/dist/index.js";

export const DEBUG = false;

export const debugLog = (...args) => {
  if (DEBUG) {
    console.log('[DEBUG]', ...args);
  }
};
