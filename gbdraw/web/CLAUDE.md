# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with the gbdraw web application.

## Overview

`index.html` is the **SPA entry point** and loads ES modules from `gbdraw/web/js` (no build step). The app runs gbdraw entirely in the browser using WebAssembly. No server is required - all genome data processing happens client-side via Pyodide (Python compiled to WebAssembly).

- **File size:** `index.html` contains HTML/CSS/templates; JavaScript lives under `gbdraw/web/js/`
- **Location:** `gbdraw/web/index.html`
- **Served by:** `gbdraw gui` command or hosted at https://gbdraw.app/

## Quick Reference

```bash
# Run locally
gbdraw gui                    # Opens browser at http://localhost:8080

# Build wheel for web deployment (version must match pyproject.toml)
python -m build
# Copy dist/gbdraw-X.X.X-py3-none-any.whl to web server

# Test in browser
# Open DevTools Console (F12) to see Pyodide output and errors
```

### Key File References

| Section | File | Description |
|---------|------|-------------|
| CSP Header | gbdraw/web/index.html | Content Security Policy |
| CSS Styles | gbdraw/web/index.html | TailwindCSS custom classes |
| Vue App Template | gbdraw/web/index.html | HTML structure |
| Core UI Logic | gbdraw/web/js/app.js | Vue setup, handlers, Pyodide integration |
| State Management | gbdraw/web/js/state.js | Reactive refs and computed |
| Components | gbdraw/web/js/components.js | HelpTip / FileUploader |
| Export Functions | gbdraw/web/js/services/export.js | PDF/PNG/SVG download |
| Config I/O | gbdraw/web/js/services/config.js | Save/load settings |
| PNG Helper | gbdraw/web/js/utils/png.js | DPI injection |

## Technology Stack

### Core Libraries (CDN)
| Library | Version | Purpose |
|---------|---------|---------|
| Vue.js 3 | 3.5.25 | Reactive UI framework |
| Pyodide | 0.29.0 | Python WebAssembly runtime |
| TailwindCSS | CDN | Utility-first CSS styling |
| jsPDF | 3.0.3 | PDF generation |
| svg2pdf.js | 2.6.0 | SVG to PDF conversion |
| DOMPurify | 3.2.7 | XSS protection for SVG output |
| Phosphor Icons | 2.1.2 | Icon library |

### Fonts
- Inter (Latin)
- Noto Sans JP (Japanese support)

## Architecture

```
index.html
├── <head>
│   ├── Content Security Policy (CSP)
│   ├── CDN script imports
│   └── TailwindCSS styles & custom CSS
├── <body>
│   ├── Vue app container (#app)
│   │   ├── Loading overlay (Pyodide init)
│   │   ├── Processing overlay
│   │   ├── Header (mode toggle, config save/load)
│   │   ├── Left Panel (4 columns)
│   │   │   ├── Input Genomes card
│   │   │   ├── Basic Settings card
│   │   │   ├── Colors & Filters (collapsible)
│   │   │   ├── Advanced Options (collapsible)
│   │   │   ├── About & Citation (collapsible)
│   │   │   └── Generate button
│   │   └── Right Panel (8 columns)
│   │       ├── Error display
│   │       ├── Result preview (SVG container with zoom)
│   │       ├── Floating controls (zoom, padding, reset)
│   │       ├── Feature Color Picker Popup
│   │       ├── Legend Editor (slide-out drawer)
│   │       └── Feature Color Editor (slide-out drawer)
│   ├── Vue component templates (HelpTip, FileUploader)
│   └── Main Vue application script
```

## Key Features

### 1. Serverless Architecture
- **Pyodide** loads a Python environment in-browser
- gbdraw wheel is installed dynamically (version must match `pyproject.toml`)
- All processing (genome parsing, SVG generation) runs locally
- **Privacy:** Genomic data never leaves the user's device

### 2. Dual Mode Support
- **Circular mode:** Single genome with track type options (tuckin/middle/spreadout)
- **Linear mode:** Multi-genome comparison with BLAST tracks

### 3. Interactive SVG Editing (Post-generation)
- **Drag & Drop:** Reposition legend and diagram elements
- **Feature Click:** Change individual feature colors with popup picker
- **Legend Editor:** Reorder, rename, delete legend entries
- **Canvas Padding:** Adjust whitespace (top/right/bottom/left)

### 4. Color Management
- **Default Colors (-d):** Base palette with 16+ presets
- **Specific Rules (-t):** Regex-based conditional coloring
- **Feature Overrides:** Per-feature color customization
- **Dual Legend Support:** Horizontal and vertical legends synchronized in linear mode

## State Management

### Key Reactive State (Vue refs)
```javascript
// System state (lines ~1339-1360)
pyodideReady            // Pyodide initialization status
processing              // Diagram generation in progress
errorLog                // Error messages

// Results (lines ~1360-1380)
results                 // Array of {name, content} SVG outputs
selectedResultIndex     // Currently viewed result
svgContent              // Computed sanitized SVG

// UI state (lines ~1380-1430)
mode                    // 'circular' | 'linear'
zoom                    // Preview zoom level
showFeaturePanel        // Feature color editor visibility
showLegendPanel         // Legend editor visibility

// Legend state (lines ~1429-1450)
legendEntries           // Extracted legend entries [{caption, color, yPos}]
addedLegendCaptions     // Set of manually added legend captions

// Feature editor state (lines ~1450-1480)
extractedFeatures       // Features from last generation
featureColorOverrides   // {featureKey: {color, caption}}
```

### Form Data Structure
```javascript
form = {
    prefix, species, strain,      // Metadata
    track_type,                   // 'tuckin' | 'middle' | 'spreadout'
    legend,                       // Position: 'right' | 'left' | 'none' | ...
    scale_style,                  // 'bar' | 'ruler' (linear)
    show_labels,                  // boolean (circular)
    show_labels_linear,           // 'none' | 'all' | 'first'
    separate_strands,
    suppress_gc, suppress_skew,   // Circular options
    show_gc, show_skew,           // Linear options
}

adv = {
    features,                     // Array of feature types to draw
    window_size, step_size, nt,   // GC calculation
    def_font_size, label_font_size,
    // Stroke settings, legend settings, etc.
}
```

## Key JavaScript Functions

### Main Entry Points
| Function | Line | Description |
|----------|------|-------------|
| `runAnalysis()` | 6813 | Main diagram generation |
| `downloadSVG()` | ~7100 | SVG export |
| `downloadPNG()` | ~7120 | PNG export via canvas |
| `downloadPDF()` | 7141 | PDF export via jsPDF |

### Legend Management
| Function | Description |
|----------|-------------|
| `getAllFeatureLegendGroups(svg)` | Get all legend groups (handles dual legends) |
| `addLegendEntry(caption, color)` | Add entry to all legend groups |
| `removeLegendEntry(caption)` | Remove entry from all legend groups |
| `extractLegendEntries()` | Extract entries for Legend Editor panel |
| `moveLegendEntryUp/Down(idx)` | Reorder entries |

### Feature Color Editing
| Function | Line | Description |
|----------|------|-------------|
| `setFeatureColor()` | 6259 | Set color and update rules |
| `applyInstantPreview()` | ~6280 | Update SVG instantly |
| `updateClickedFeatureColor()` | 1927 | Handle popup color change |

### Drag & Drop
| Function | Description |
|----------|-------------|
| `startLegendDrag(e)` / `onLegendDrag(e)` / `endLegendDrag()` | Legend positioning |
| `startDiagramDrag(e)` / `onDiagramDrag(e)` / `endDiagramDrag()` | Diagram element positioning |
| `parseTransform(str)` | Extract x,y from transform attribute |
| `resetAllPositions()` | Reset to original positions |

## Python Integration

### Pyodide Initialization Flow
1. Load Pyodide runtime
2. Install micropip
3. Load gbdraw wheel from same origin
4. Install Python dependencies (biopython, svgwrite, pandas, fonttools, bcbio-gff)
5. Initialize helper functions
6. Mark `pyodideReady = true`

### Python Helper Functions (defined in index.html)
```python
get_palettes_json()              # Load color palettes from TOML
run_gbdraw_wrapper()             # Execute circular or linear mode
generate_legend_entry_svg()      # Generate SVG for dynamic legend entry
extract_features_from_genbank()  # Extract features for UI color editor
```

### File System Access Pattern
```javascript
// Write file to Pyodide virtual FS (line 6806)
const writeToFS = async (fileObj, path) => {
    const buffer = await fileObj.arrayBuffer();
    pyodide.FS.writeFile(path, new Uint8Array(buffer));
};

// Read Python results via JSON
const resultJson = pyodide.runPython("run_gbdraw_wrapper('circular', args)");
const results = JSON.parse(resultJson);
```

## Security Considerations

### Content Security Policy (lines 5-17)
```
default-src 'self';
script-src 'self' 'unsafe-inline' 'unsafe-eval' https://unpkg.com https://cdn.tailwindcss.com ...;
style-src 'self' 'unsafe-inline' https://fonts.googleapis.com ...;
connect-src 'self' https://cdn.jsdelivr.net https://pypi.org https://files.pythonhosted.org ...;
frame-ancestors 'none';
```

### SVG Sanitization
All generated SVG passes through DOMPurify:
- `FORBID_TAGS: ['style', 'script', 'foreignObject', 'iframe', 'animate', ...]`
- `FORBID_ATTR: ['onload', 'onclick', 'onerror', ...]`

### Regex Validation
User-provided regex patterns are validated:
- Length check (warn if >50 chars)
- ReDoS pattern detection
- Try-catch around `new RegExp()`

## Development Notes

### Modifying the UI
- CSS: Tailwind utility classes inline + custom classes in `<style>` block
- Custom classes: `.card`, `.btn-*`, `.form-input`, etc. (lines 43-72)
- Colors: Slate palette with Blue/Indigo accents
- Collapsible sections: `<details>` element with custom styling

### Adding New Settings
1. Add reactive state to `setup()` function (~line 1339)
2. Add form element in appropriate card
3. Include in args array passed to Python (~line 6840)
4. Update Python-side handling if needed

### Debugging
- Browser DevTools Console shows Pyodide output
- `[DEBUG]` prefixed logs for incremental edit tracking
- Vue DevTools compatible
- Check `errorLog.value` for Python exceptions

### Common Issues
1. **"Pyodide not ready"**: Wait for initialization (~5-15 seconds on first load)
2. **Memory errors**: Large genomes may exhaust browser memory
3. **Wheel version mismatch**: Ensure wheel version matches `pyproject.toml`
4. **CSP errors**: Check CDN URLs in CSP header if adding new dependencies

## File Dependencies

When deploying:
- `index.html` (this file)
- `gbdraw-X.X.X-py3-none-any.whl` (Python wheel, same origin)
- CDN dependencies (Vue, Pyodide, TailwindCSS, icons, jsPDF, DOMPurify)

## Known Limitations

1. **Memory:** Large genomes may exhaust browser memory
2. **Performance:** Initial Pyodide load takes 5-15 seconds
3. **Fonts:** Custom fonts require paths accessible to Pyodide
4. **CairoSVG:** Not available in browser (PNG/PDF use canvas instead)
5. **Incremental Updates:** Some edits require full regeneration

## Related Files

- `gbdraw/cli.py` - CLI entry point with `gui` command
- `gbdraw/data/config.toml` - Default configuration values
- `gbdraw/data/color_palettes.toml` - Color palette definitions
- Parent project: See `CLAUDE.md` in project root
