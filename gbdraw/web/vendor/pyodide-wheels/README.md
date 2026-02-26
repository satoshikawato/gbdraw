# Pyodide Wheels (Offline)

Place all Python dependency wheels used by `gbdraw gui` in this directory.

The app reads required wheel filenames from:
`gbdraw/web/js/config.js` -> `PYODIDE_REQUIRED_WHEELS`.

Use versioned wheel filenames (no short aliases), and keep this directory
strictly synchronized with `PYODIDE_REQUIRED_WHEELS`.

Each listed file must exist here and be a valid `.whl` archive (`PK` header).
