# Third-Party Licenses (Web Vendor Assets)

This file tracks licenses for third-party assets vendored under `gbdraw/web/vendor`.

## License summary

- Vue.js (`vue@3.5.25`): MIT
- Tailwind CSS (generated utility stylesheet): MIT
- Pyodide core/runtime (`0.29.0`): MPL-2.0
- Pyodide package wheels (from Pyodide lock): package-specific OSS licenses (see `manifest.json`)
- `bcbio_gff-0.7.1`: Biopython License
- @phosphor-icons/web (`2.1.2`, regular assets + local loader adapter): MIT
- jsPDF (`3.0.3`): MIT
- svg2pdf.js (`2.6.0`): MIT
- DOMPurify (`3.2.7`): Apache-2.0 OR MPL-2.0
- @bjorn3/browser_wasi_shim (`0.4.2`): MIT OR Apache-2.0
- Inter (`google/fonts`): SIL Open Font License 1.1
- Noto Sans JP (`google/fonts`): SIL Open Font License 1.1

## Source of truth

- Asset inventory, source URL, version, and SHA-256:
  - `gbdraw/web/vendor/manifest.json`
- Wheel set consumed by Pyodide:
  - `gbdraw/web/js/config.js` (`PYODIDE_REQUIRED_WHEELS`)
