# CLAUDE.md - gbdraw Web App (index.html)

This document provides guidance for working with the gbdraw serverless web application.

## Overview

`index.html` is a **self-contained single-page application (SPA)** that runs gbdraw entirely in the browser using WebAssembly. No server is required - all genome data processing happens client-side via Pyodide (Python compiled to WebAssembly).

- **File size:** ~3500+ lines (HTML + embedded JavaScript + CSS)
- **Location:** `gbdraw/web/index.html`
- **Served by:** `gbdraw gui` command or hosted at https://gbdraw.app/

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
│   │   │   └── Generate button
│   │   └── Right Panel (8 columns)
│   │       ├── Error display
│   │       ├── Result preview (SVG container)
│   │       ├── Zoom controls
│   │       ├── Canvas padding controls
│   │       ├── Feature Color Editor (drawer)
│   │       └── Legend Editor (drawer)
│   ├── Vue component templates
│   │   ├── HelpTip (tooltip component)
│   │   └── FileUploader (drag-drop file input)
│   └── Main Vue application script
```

## Key Features

### 1. Serverless Architecture
- **Pyodide** loads a Python environment in-browser
- gbdraw wheel (`gbdraw-0.8.3-py3-none-any.whl`) is installed dynamically
- All processing (genome parsing, SVG generation) runs locally
- **Privacy:** Genomic data never leaves the user's device

### 2. Dual Mode Support
- **Circular mode:** Single genome visualization
- **Linear mode:** Multi-genome comparison with BLAST tracks

### 3. Interactive SVG Editing
Post-generation interactive editing without regeneration:
- **Drag & Drop:** Reposition legend and diagram elements
- **Feature Click:** Change individual feature colors with popup picker
- **Legend Editor:** Reorder, rename, delete legend entries
- **Canvas Padding:** Adjust whitespace around diagram
- **Real-time Color Sync:** Palette changes update SVG instantly

### 4. Color Management
- **Default Colors (-d):** Base palette with 16+ presets
- **Specific Rules (-t):** Regex-based conditional coloring
- **Feature Overrides:** Per-feature color customization
- Automatic legend entry creation for custom colors

### 5. Export Options
| Format | Method | DPI Options |
|--------|--------|-------------|
| SVG | Direct download | N/A |
| PNG | Canvas rendering | 72-600 DPI |
| PDF | jsPDF + svg2pdf.js | Selected DPI |

## State Management

### Reactive State (Vue refs)
```javascript
// System state
pyodideReady        // Pyodide initialization status
processing          // Diagram generation in progress
loadingStatus       // Current loading step message
errorLog            // Error messages

// Results
results             // Array of {name, content} SVG outputs
selectedResultIndex // Currently viewed result
svgContent          // Computed sanitized SVG

// UI state
mode                // 'circular' | 'linear'
zoom                // Preview zoom level (0.1 - 2.0+)
showFeaturePanel    // Feature color editor visibility
showLegendPanel     // Legend editor visibility
showCanvasControls  // Canvas padding panel visibility
```

### Form Data
```javascript
form = {
    prefix, species, strain,      // Metadata
    track_type,                   // 'tuckin' | 'middle' | 'spreadout'
    legend,                       // Position: 'right' | 'left' | 'top' | 'bottom' | 'none'
    show_labels, separate_strands,
    suppress_gc, suppress_skew,   // Circular options
    show_gc, show_skew,           // Linear options
    align_center, normalize_length
}

adv = {
    features,                     // Feature types to draw
    window_size, step_size, nt,   // GC calculation
    def_font_size, label_font_size,
    block_stroke_*, line_stroke_*, axis_stroke_*,  // Styling
    legend_box_size, legend_font_size,
    // Linear-specific
    feature_height, gc_height, comparison_height,
    min_bitscore, evalue, identity
}
```

## Security Considerations

### Content Security Policy
Strict CSP headers limit resource loading:
- Scripts: Self + specific CDNs only
- Styles: Self + Google Fonts + CDNs
- Images: Self + data: + blob:
- Workers: Self + blob:
- Frames: None (frame-ancestors 'none')

### SVG Sanitization
All generated SVG passes through DOMPurify:
```javascript
DOMPurify.sanitize(rawSvg, {
    USE_PROFILES: { svg: true },
    FORBID_TAGS: ['style', 'script', 'foreignObject', 'iframe', ...],
    FORBID_ATTR: ['onload', 'onclick', 'onerror', ...]
});
```

### Regex Validation
User-provided regex patterns are validated:
- Length check (warn if >50 chars)
- ReDoS pattern detection
- Try-catch around `new RegExp()`

## Python Integration

### Pyodide Initialization Flow
1. Load Pyodide runtime
2. Install micropip
3. Load gbdraw wheel from same origin
4. Install Python dependencies (biopython, svgwrite, pandas, fonttools, bcbio-gff)
5. Initialize gbdraw module
6. Mark `pyodideReady = true`

### Python Execution Pattern
```javascript
// Run Python code
pyodide.runPython(`
    from gbdraw.api.diagram import assemble_circular_diagram_from_record
    result = assemble_circular_diagram_from_record(record, config)
    svg_string = result.tostring()
`);

// Access Python results
const svgOutput = pyodide.globals.get('svg_string');
```

## Custom Vue Components

### HelpTip
Tooltip component with smart positioning:
```html
<help-tip text="Explanation text here"></help-tip>
```

### FileUploader
File input with drag-drop support:
```html
<file-uploader
    label="GenBank File"
    accept=".gb,.gbk"
    v-model="files.c_gb"
    :small="true">
</file-uploader>
```

## Development Notes

### Modifying the UI
- CSS: Uses Tailwind utility classes inline
- Custom classes defined in `<style>` block: `.card`, `.btn-*`, `.form-input`, etc.
- Colors follow Slate palette with Blue/Indigo accents

### Adding New Settings
1. Add reactive state to `setup()` function
2. Add form element in appropriate card
3. Include in config object passed to Python
4. Update Python-side handling

### Debugging
- Browser console shows Pyodide output
- `console.log` statements throughout for key operations
- Vue DevTools compatible

## File Dependencies

When deploying:
- `index.html` (this file)
- `gbdraw-X.X.X-py3-none-any.whl` (Python wheel, served from same origin)
- CDN dependencies (Vue, Pyodide, TailwindCSS, icons, jsPDF)

## Known Limitations

1. **Memory:** Large genomes may exhaust browser memory
2. **Performance:** Initial Pyodide load takes 5-15 seconds
3. **Fonts:** Custom fonts require paths accessible to Pyodide
4. **CairoSVG:** Not available in browser (PNG/PDF use canvas instead)

## Related Files

- `gbdraw/cli.py` - CLI entry point with `gui` command
- `app.py` - Flask wrapper for development server
- `gbdraw/data/config.toml` - Default configuration values
- `gbdraw/data/color_palettes.toml` - Color palette definitions
