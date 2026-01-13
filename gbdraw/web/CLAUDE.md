# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with the gbdraw web application.

## Overview

`index.html` is a **self-contained single-page application (SPA)** that runs gbdraw entirely in the browser using WebAssembly. No server is required - all genome data processing happens client-side via Pyodide (Python compiled to WebAssembly).

- **Location:** `gbdraw/web/index.html`
- **Served by:** `gbdraw gui` command or hosted at https://gbdraw.app/
- **Constraint:** Must remain a single file with no build step

## Quick Reference

```bash
# Run locally
gbdraw gui                    # Opens browser at http://localhost:<port>

# Build wheel for web deployment
python -m build
# Wheel version in GBDRAW_WHEEL_NAME constant must match pyproject.toml

# Debug in browser
# Open DevTools Console (F12) to see Pyodide output and errors
```

## Technology Stack

| Library | Version | Purpose |
|---------|---------|---------|
| Vue.js 3 | 3.5.25 | Reactive UI framework |
| Pyodide | 0.29.0 | Python WebAssembly runtime |
| TailwindCSS | CDN | Utility-first CSS styling |
| jsPDF | 3.0.3 | PDF generation |
| svg2pdf.js | 2.6.0 | SVG to PDF conversion |
| DOMPurify | 3.2.7 | XSS protection for SVG output |
| Phosphor Icons | 2.1.2 | Icon library |

## Architecture

```
index.html
├── <head>
│   ├── Content Security Policy (CSP) - lines 5-17
│   ├── CDN script imports with SRI hashes - lines 20-32
│   └── TailwindCSS styles & custom CSS - lines 33-73
├── <body>
│   ├── Vue app container (#app) - lines 77-1095
│   │   ├── Loading/Processing overlays
│   │   ├── Header (mode toggle, config save/load)
│   │   ├── Left Panel: Input, Settings, Colors, Advanced
│   │   └── Right Panel: Preview, Controls, Editors
│   ├── Vue component templates (HelpTip, FileUploader) - lines 1097-1127
│   └── Main Vue application script - lines 1129-7629
```

## Key Code Sections

Use these approximate line ranges to navigate the codebase:

| Section | Description |
|---------|-------------|
| **State Management** (~1236-1480) | Reactive refs: `pyodideReady`, `processing`, `results`, `mode`, etc. |
| **Form Data** (~1353-1400) | `form` and `adv` reactive objects for settings |
| **Feature Color Editing** (~1878-2150) | `applyInstantPreview()`, `setFeatureColor()` |
| **Legend Management** (~3162-4500) | `addLegendEntry()`, `removeLegendEntry()`, `extractLegendEntries()` |
| **Drag & Drop** (~5137-5400) | Legend and diagram element repositioning |
| **Pyodide Init** (~6900-7100) | Python environment setup and helper functions |
| **Main Entry** (~7112) | `runAnalysis()` - diagram generation |
| **Export Functions** (~7397-7500) | `downloadSVG()`, `downloadPNG()`, `downloadPDF()` |

## Key Features

### Serverless Architecture
- Pyodide loads Python environment in-browser
- gbdraw wheel installed dynamically (version in `GBDRAW_WHEEL_NAME` constant ~line 1131)
- All processing runs locally - genomic data never leaves user's device

### Dual Mode Support
- **Circular:** Single genome with track type options (tuckin/middle/spreadout)
- **Linear:** Multi-genome comparison with BLAST tracks

### Interactive SVG Editing
- Drag & drop legend and diagram elements
- Click features to change colors with popup picker
- Legend editor for reordering, renaming, deleting entries
- Canvas padding adjustment (top/right/bottom/left)

## Python Integration

### Initialization Flow (~6900-7100)
1. Load Pyodide runtime
2. Install micropip
3. Install Python dependencies (biopython, svgwrite, pandas, fonttools, bcbio-gff)
4. Load gbdraw wheel from same origin
5. Define helper functions
6. Mark `pyodideReady = true`

### Python Helper Functions
```python
get_palettes_json()              # Load color palettes from TOML (~6935)
run_gbdraw_wrapper()             # Execute circular or linear mode (~6941)
generate_legend_entry_svg()      # Generate SVG for dynamic legend entry (~6964)
extract_features_from_genbank()  # Extract features for UI color editor (~7033)
```

### File System Access Pattern
```javascript
// Write file to Pyodide virtual FS
pyodide.FS.writeFile(path, new Uint8Array(buffer));

// Read Python results via JSON
const resultJson = pyodide.runPython("run_gbdraw_wrapper('circular', args)");
const results = JSON.parse(resultJson);
```

## Security

### Content Security Policy (lines 5-17)
Restricts script/style sources to trusted CDNs. Update CSP when adding new dependencies.

### SVG Sanitization
All generated SVG passes through DOMPurify with forbidden tags/attributes for XSS protection.

### Regex Validation
User-provided regex patterns are validated for length and ReDoS patterns.

## Development Notes

### Adding New Settings
1. Add reactive state to `setup()` function (~line 1232)
2. Add form element in appropriate card
3. Include in args array passed to Python (~line 7200)
4. Update Python-side handling if needed

### Modifying the UI
- CSS: Tailwind utility classes inline + custom classes in `<style>` block
- Custom classes: `.card`, `.btn-*`, `.form-input`, etc. (lines 43-72)
- Colors: Slate palette with Blue/Indigo accents

### Debugging
- Browser DevTools Console shows Pyodide output
- Set `DEBUG = true` (~line 1134) for verbose logging
- Check `errorLog.value` for Python exceptions

### Common Issues
1. **"Pyodide not ready"**: Wait for initialization (~5-15 seconds on first load)
2. **Memory errors**: Large genomes may exhaust browser memory
3. **Wheel version mismatch**: Ensure `GBDRAW_WHEEL_NAME` matches `pyproject.toml`
4. **CSP errors**: Check CDN URLs in CSP header if adding new dependencies

## Known Limitations

1. **Memory:** Large genomes may exhaust browser memory
2. **Performance:** Initial Pyodide load takes 5-15 seconds
3. **CairoSVG:** Not available in browser (PNG/PDF use canvas instead)
4. **Fonts:** Custom fonts require paths accessible to Pyodide

## Related Files

- `gbdraw/cli.py` - CLI entry point with `gui` command
- `gbdraw/data/config.toml` - Default configuration values
- `gbdraw/data/color_palettes.toml` - Color palette definitions
- Parent project: See `CLAUDE.md` in project root
