# Web Vendor Assets (Offline Bundle)

This directory contains third-party frontend assets used by `gbdraw gui`.

All runtime assets required by the web UI are vendored here for offline execution.

## Required assets

- `vue/vue.global.js`
- `tailwind/app.min.css` (generated from `../css/tailwind.source.css` via `python tools/build_web_css.py`)
- `pyodide/v0.29.0/full/*` (Pyodide core runtime files + `micropip` wheel)
- `pyodide-wheels/*.whl` listed in `gbdraw/web/js/config.js`
- `phosphor/index.js` and `phosphor/regular/*` (offline icon CSS/fonts)
- `jspdf/jspdf.umd.min.js`
- `svg2pdf/svg2pdf.umd.min.js`
- `dompurify/purify.min.js`
- `browser_wasi_shim/*.js` (all `dist` modules required by `index.js`)
- `fonts/Inter-Variable.ttf` and `fonts/NotoSansJP-Variable.ttf`

## Verification checklist

1. `rg -n \"Placeholder vendor asset\" gbdraw/web/vendor` returns no matches.
2. Every filename in `PYODIDE_REQUIRED_WHEELS` exists in `pyodide-wheels/`.
3. `python tools/build_web_css.py --check` succeeds.
4. `pytest tests/test_web_offline_assets.py -v` succeeds.

See `manifest.json` for exact source URLs, versions, SHA-256 values, and license tags.
