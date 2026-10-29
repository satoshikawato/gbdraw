// @ts-check
export const GBDRAW_WHEEL_NAME = "gbdraw-0.14.4-py3-none-any.whl";
// The gbdraw version: the version field of the wheel's file name (PEP 427),
// which tools/prepare_browser_wheel.py writes from pyproject.toml.
export const GBDRAW_VERSION = GBDRAW_WHEEL_NAME.split('-')[1];
// The `?v=` query of the wheel URL: the wheel's content hash, which
// tools/prepare_cloudflare_pages.py writes into the hosted bundle's copy of
// this file. Null here: the package (`gbdraw gui`) and the test servers load
// the wheel by its plain name.
/** @type {string | null} */
export const GBDRAW_WHEEL_CACHE_BUST = null;
// The folder of a Pyodide snapshot the diagram Worker restores instead of
// installing its packages (workers/pyodide-snapshot.js). Null here: only a test
// server sets it (docs/internal/SELECTIVE_CI.md, "Pyodide snapshot").
/** @type {string | null} */
export const PYODIDE_SNAPSHOT_URL = null;
// tools/prepare_browser_wheel.py writes the four PYODIDE_ values below from
// gbdraw/_web_assets.py (PYODIDE_RUNTIME_DIR, PYODIDE_LOADER_FILES,
// PYODIDE_LOCAL_WHEELS, PYODIDE_RUNTIME_WHEELS); tests/test_web_packaging.py
// fails when they differ.
export const PYODIDE_INDEX_URL = "./vendor/pyodide/v0.29.0/full/";
// The files in PYODIDE_INDEX_URL that Pyodide's loader downloads at start.
export const PYODIDE_RUNTIME_FILES = [
  "pyodide.asm.js",
  "pyodide.asm.wasm",
  "python_stdlib.zip",
  "pyodide-lock.json"
];
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
// The Pyodide distribution's wheels the diagram Worker loads by name before
// micropip installs PYODIDE_LOCAL_WHEELS: micropip and what those wheels
// require from the distribution (pandas' tzdata).
export const PYODIDE_RUNTIME_WHEELS = [
  "./vendor/pyodide/v0.29.0/full/micropip-0.11.0-py3-none-any.whl",
  "./vendor/pyodide/v0.29.0/full/tzdata-2025.2-py2.py3-none-any.whl"
];
export const LOSAT_THREADED_WASM_URL = "./wasm/losat/losat-threaded.wasm";
export const WASI_SHIM_URL = "./vendor/browser_wasi_shim/dist/index.js";
