# CLAUDE.md - gbdraw Web App (index.html)

This document provides guidance for working with the gbdraw serverless web application.

## Overview

`index.html` is a **self-contained single-page application (SPA)** that runs gbdraw entirely in the browser using WebAssembly. No server is required - all genome data processing happens client-side via Pyodide (Python compiled to WebAssembly).

- **File size:** ~4500+ lines (HTML + embedded JavaScript + CSS)
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
│   │   │   │   ├── Circular: GenBank or GFF3+FASTA
│   │   │   │   └── Linear: Multi-sequence with BLAST tracks
│   │   │   ├── Basic Settings card
│   │   │   │   ├── Output prefix, Legend position
│   │   │   │   ├── Track layout (circular only)
│   │   │   │   ├── Scale style (linear only)
│   │   │   │   └── Checkboxes (labels, strands, GC, etc.)
│   │   │   ├── Colors & Filters (collapsible)
│   │   │   │   ├── Default Colors (-d) with palette picker
│   │   │   │   ├── Specific Rules (-t) with regex support
│   │   │   │   ├── Label Filtering (Blacklist/Whitelist)
│   │   │   │   └── Qualifier Priority
│   │   │   ├── Advanced Options (collapsible)
│   │   │   │   ├── GC window/step settings
│   │   │   │   ├── Font sizes, stroke styles
│   │   │   │   ├── Legend styling
│   │   │   │   ├── Linear-specific settings (heights, BLAST filters, scale)
│   │   │   │   ├── Circular-specific offsets
│   │   │   │   └── Feature type selection
│   │   │   ├── About & Citation (collapsible)
│   │   │   └── Generate button
│   │   └── Right Panel (8 columns)
│   │       ├── Error display
│   │       ├── Result preview (SVG container with zoom)
│   │       ├── Floating controls
│   │       │   ├── Reset positions button
│   │       │   ├── Canvas padding toggle
│   │       │   └── Zoom controls
│   │       ├── Canvas padding controls panel
│   │       ├── Feature Color Picker Popup (on click)
│   │       ├── Legend Editor (slide-out drawer)
│   │       └── Feature Color Editor (slide-out drawer)
│   ├── Vue component templates
│   │   ├── HelpTip (tooltip component with smart positioning)
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
- **Circular mode:** Single genome visualization with track type options (tuckin/middle/spreadout)
- **Linear mode:** Multi-genome comparison with BLAST tracks and scale style options (bar/ruler)

### 3. Interactive SVG Editing
Post-generation interactive editing without regeneration:
- **Drag & Drop:** Reposition legend and diagram elements
- **Feature Click:** Change individual feature colors with popup picker (includes legend name option)
- **Legend Editor:** Reorder (up/down), rename, delete legend entries
- **Canvas Padding:** Adjust whitespace around diagram (top/right/bottom/left)
- **Real-time Color Sync:** Palette changes update SVG instantly
- **Position Reset:** Reset all element positions to original

### 4. Color Management
- **Default Colors (-d):** Base palette with 16+ presets, dynamically editable
- **Specific Rules (-t):** Regex-based conditional coloring with optional legend captions
- **Feature Overrides:** Per-feature color customization with automatic legend entry creation
- **Dual Legend Support:** Both horizontal and vertical legends synchronized in linear mode

### 5. Label Filtering
- **Blacklist Mode:** Exclude labels containing specified keywords
- **Whitelist Mode:** Only show labels matching feature/qualifier/keyword rules
- **Qualifier Priority:** Define label priority order per feature type (e.g., CDS → product,gene,locus_tag)

### 6. Export Options
| Format | Method | DPI Options |
|--------|--------|-------------|
| SVG | Direct download | N/A |
| PNG | Canvas rendering | 72, 96, 150, 300, 600 DPI |
| PDF | jsPDF + svg2pdf.js | Selected DPI |

## State Management

### Reactive State (Vue refs)
```javascript
// System state
pyodideReady            // Pyodide initialization status
processing              // Diagram generation in progress
loadingStatus           // Current loading step message
errorLog                // Error messages

// Results
results                 // Array of {name, content} SVG outputs
selectedResultIndex     // Currently viewed result
svgContent              // Computed sanitized SVG
pairwiseMatchFactors    // {pathId: factor} for match re-interpolation

// UI state
mode                    // 'circular' | 'linear'
zoom                    // Preview zoom level (0.1 - 2.0+)
showFeaturePanel        // Feature color editor visibility
showLegendPanel         // Legend editor visibility
showCanvasControls      // Canvas padding panel visibility
clickedFeature          // {id, svg_id, label, color, legendName} for popup
clickedFeaturePos       // {x, y} popup position

// Legend state
legendEntries           // Extracted legend entries [{caption, color, yPos}]
circularLegendPosition  // Saved legend position for circular mode
linearLegendPosition    // Saved legend position for linear mode
addedLegendCaptions     // Set of manually added legend captions

// Feature editor state
extractedFeatures       // Features from last generation
featureRecordIds        // Record IDs for multi-record files
selectedFeatureRecordIdx // Currently selected record index
featureSearch           // Search filter text
featureColorOverrides   // {featureKey: {color, caption}}

// Canvas state
canvasPadding           // {top, right, bottom, left}
diagramElements         // Draggable diagram elements
diagramElementBaseTransforms // Base positions for drag offset
skipCaptureBaseConfig   // Flag to prevent re-capture on internal updates
```

### Form Data
```javascript
form = {
    prefix, species, strain,      // Metadata
    track_type,                   // 'tuckin' | 'middle' | 'spreadout' (circular)
    legend,                       // Position: 'right' | 'left' | 'top' | 'bottom' | 'upper_left' | 'upper_right' | 'none'
    scale_style,                  // 'bar' | 'ruler' (linear)
    show_labels,                  // boolean (circular)
    show_labels_linear,           // 'none' | 'all' | 'first' (linear)
    separate_strands,
    allow_inner_labels,           // Circular option for inner labels
    suppress_gc, suppress_skew,   // Circular options
    show_gc, show_skew,           // Linear options
    align_center, normalize_length // Linear options
}

adv = {
    features,                     // Array of feature types to draw
    window_size, step_size, nt,   // GC calculation
    def_font_size, label_font_size,
    block_stroke_color, block_stroke_width,
    line_stroke_color, line_stroke_width,
    axis_stroke_color, axis_stroke_width,
    legend_box_size, legend_font_size,
    // Linear-specific
    resolve_overlaps,             // Experimental overlap resolution
    feature_height, gc_height, comparison_height,
    min_bitscore, evalue, identity,
    scale_interval, scale_font_size, scale_stroke_width, scale_stroke_color,
    // Circular-specific
    outer_label_x_offset, outer_label_y_offset,
    inner_label_x_offset, inner_label_y_offset
}

// Manual rules
manualSpecificRules     // [{feat, qual, val, color, cap}]
manualPriorityRules     // [{feat, order}]
manualBlacklist         // Comma-separated string
manualWhitelist         // [{feat, qual, key}]
filterMode              // 'None' | 'Blacklist' | 'Whitelist'
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
    ADD_TAGS: ['use', 'g', 'defs', 'linearGradient', 'radialGradient', 'stop', ...],
    ADD_ATTR: ['xlink:href', 'href', 'transform', 'viewBox', ...],
    FORBID_TAGS: ['style', 'script', 'foreignObject', 'iframe', 'animate', ...],
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
5. Initialize gbdraw module with helper functions
6. Mark `pyodideReady = true`

### Python Helper Functions
```python
get_palettes_json()        # Load color palettes from TOML
run_gbdraw_wrapper()       # Execute circular or linear mode
generate_legend_entry_svg() # Generate SVG for dynamic legend entry
extract_features_from_genbank() # Extract features for UI color editor
```

### Python Execution Pattern
```javascript
// Run Python code
pyodide.runPython(`
    from gbdraw.circular import circular_main
    result = circular_main(args)
`);

// Access Python results via JSON
const resultJson = pyodide.runPython("run_gbdraw_wrapper('circular', args)");
const results = JSON.parse(resultJson);
```

## Custom Vue Components

### HelpTip
Tooltip component with smart positioning (avoids viewport overflow):
```html
<help-tip text="Explanation text here"></help-tip>
```

### FileUploader
File input with drag-drop support:
```html
<file-uploader
    label="GenBank File (.gb)"
    accept=".gb,.gbk,.txt"
    v-model="files.c_gb"
    :small="true">
</file-uploader>
```

## Key JavaScript Functions

### Legend Management
- `getAllFeatureLegendGroups(svg)` - Get all legend groups (handles dual legends)
- `getVisibleFeatureLegendGroup(svg)` - Get visible legend for extraction
- `addLegendEntry(caption, color)` - Add entry to all legend groups
- `removeLegendEntry(caption)` - Remove entry from all legend groups
- `updateLegendEntryColorByCaption(caption, color)` - Update color without removing
- `extractLegendEntries()` - Extract entries for Legend Editor panel
- `moveLegendEntryUp/Down(idx)` - Reorder entries
- `swapLegendEntries(idx1, idx2)` - Swap Y positions in SVG

### Drag & Drop
- `startLegendDrag(e)` / `onLegendDrag(e)` / `endLegendDrag()`
- `startDiagramDrag(e)` / `onDiagramDrag(e)` / `endDiagramDrag()`
- `parseTransform(str)` - Extract x,y from transform attribute
- `resetAllPositions()` / `resetLegendPosition()` - Reset to original

### Feature Color Editing
- `setFeatureColor(feat, color, customCaption)` - Set color and update rules
- `applyInstantPreview(feat, color, caption)` - Update SVG instantly
- `getFeatureColor(feat)` / `canEditFeatureColor(feat)` - Color lookup
- `updateClickedFeatureColor(color)` - Handle popup color change

### Canvas & Export
- `updateCanvasPadding()` - Apply padding to SVG viewBox
- `downloadSVG()` / `downloadPNG()` / `downloadPDF()` - Export functions
- `setDpiInPng(blob, dpi)` - Inject pHYs chunk for DPI metadata

## Development Notes

### Modifying the UI
- CSS: Uses Tailwind utility classes inline
- Custom classes defined in `<style>` block: `.card`, `.btn-*`, `.form-input`, etc.
- Colors follow Slate palette with Blue/Indigo accents
- Collapsible sections use `<details>` element with custom styling

### Adding New Settings
1. Add reactive state to `setup()` function
2. Add form element in appropriate card
3. Include in config object / args array passed to Python
4. Update Python-side handling if needed

### Debugging
- Browser console shows Pyodide output
- `console.log` statements throughout for key operations
- Vue DevTools compatible
- `[DEBUG]` prefixed logs for incremental edit tracking

## File Dependencies

When deploying:
- `index.html` (this file)
- `gbdraw-0.8.3-py3-none-any.whl` (Python wheel, served from same origin)
- CDN dependencies (Vue, Pyodide, TailwindCSS, icons, jsPDF, DOMPurify)

## Known Limitations

1. **Memory:** Large genomes may exhaust browser memory
2. **Performance:** Initial Pyodide load takes 5-15 seconds
3. **Fonts:** Custom fonts require paths accessible to Pyodide
4. **CairoSVG:** Not available in browser (PNG/PDF use canvas instead)
5. **Incremental Updates:** Some edits require full regeneration

## Related Files

- `gbdraw/cli.py` - CLI entry point with `gui` command
- `app.py` - Flask wrapper for development server
- `gbdraw/data/config.toml` - Default configuration values
- `gbdraw/data/color_palettes.toml` - Color palette definitions
