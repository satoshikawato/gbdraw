[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Gallery](./GALLERY.md) | [Recipes](./RECIPES.md) | [FAQ](./FAQ.md) | [ABOUT](./ABOUT.md)

# Web App Guide

The gbdraw web app runs entirely in your browser using WebAssembly (Pyodide). Your genome files are processed locally and are not uploaded to a server.

## Access

### Hosted Web App
Open the hosted app at:
https://gbdraw.app/

### Local Web App (CLI)
Run the local web UI from an installed gbdraw package:
```bash
gbdraw gui
```
This serves the app at `http://localhost:<port>` and opens your browser. An internet connection is still required to load CDN dependencies.

## Browser Requirements

- Modern browser with WebAssembly support (Chrome, Edge, Firefox, Safari).
- Sufficient memory for your genome size (large genomes may exceed browser limits).
- File upload enabled (drag and drop or file picker).

## Workflow

1. Choose diagram mode (circular or linear).
2. Upload inputs (GenBank/DDBJ or GFF3 + FASTA).
3. Optional: add BLAST files for linear comparisons (outfmt 6 or 7).
4. Adjust styling, tracks, labels, and legend settings.
5. Click Generate Diagram and preview the results.
6. Export results (SVG, PNG, PDF).

## Tutorials

### Tutorial 1: First circular diagram (GenBank)
1. Click the Circular mode toggle in the header.
2. In the Input section, upload a `.gb`/`.gbk` file.
3. Choose a Legend position and Track Layout if needed.
4. Enable Show Labels if you want annotations on the map.
5. Click Generate Diagram.
6. Use the preview zoom controls and drag to reposition if needed.
7. Download SVG, PNG, or PDF from the result toolbar.

### Tutorial 2: Linear comparison with BLAST
1. Switch to Linear mode.
2. Upload 2 or more genomes (GenBank or GFF3 + FASTA).
3. Upload BLAST files in outfmt 6 or 7.
4. Ensure the BLAST file order matches the genome order (A vs B, B vs C).
5. If ribbons are missing, adjust BLAST Filters in Advanced Options (bitscore, e-value, identity).
6. Click Generate Diagram and select a result if multiple outputs appear.

### Tutorial 3: Colors and feature rules
1. Open the Colors & Filters panel.
2. Pick a palette from the Default Colors dropdown.
3. Click a color chip to edit the base color for a feature type.
4. Add a custom feature color with Add Feature and a color picker.
5. Add Specific Rules (-t) to color features by qualifier or regex.
6. Optional: upload -d or -t TSV files if you already have them.

### Tutorial 4: Label control
1. Turn on Show Labels.
2. In Label Filtering, choose Blacklist or Whitelist.
3. For Blacklist, add common terms like "hypothetical".
4. For Whitelist, add rules (feature, qualifier, keyword).
5. Use Qualifier Priority to prefer gene names over product names.
6. Adjust Label Font Size in Advanced Options if labels are crowded.

### Tutorial 5: Legend and layout edits
1. Click the Edit legend button to open the Legend Editor.
2. Rename, reorder, recolor, or delete legend entries.
3. Drag the legend or diagram to reposition on the canvas.
4. Use the Canvas Padding control to adjust margins.
5. Click Reset all element positions if needed.

### Tutorial 6: Save and reuse settings
1. Click Save Config to download a JSON settings file.
2. Use Load Config to restore settings later.
3. Re-upload genome and BLAST files (configs do not include files).
4. Regenerate the diagram after loading.

## Screenshots and Release Capture

### Capture setup (recommended)
1. Use Chrome or Edge with 100% zoom and OS scale.
2. Set the window size to 1440x900 (or keep it consistent across all shots).
3. Use a small dataset from `examples/` for fast rendering.
4. Wait until the app finishes loading (Pyodide ready) before capturing.
5. Capture each required screen below in a clean, uncluttered state.

### Required shots checklist
- Header with mode toggle and Save/Load Config buttons
- Input files selected (GenBank or GFF3 + FASTA)
- Colors & Filters panel (palette + specific rules)
- Advanced Options (GC/labels and mode-specific settings)
- Processing state (Generate Diagram -> Processing)
- Result preview with DPI selector and export buttons
- Legend Editor open (rename/reorder example)
- Drag and reposition example (before/after if possible)
- Canvas padding panel open
- Feature color edit (clicked feature with updated color)

### Annotation template
See `docs/WEBAPP_SCREENSHOTS.md` for file naming and caption library.

```text
File: webapp_01_overview.png
Title: Web App Overview
Purpose: Show mode toggle and config actions
Steps: 1) Toggle mode 2) Open Save/Load config
Callouts:
  [1] Mode toggle
  [2] Save Config
  [3] Load Config
Notes: Browser zoom 100%, window 1440x900
Alt: gbdraw web app header with mode toggle and config buttons
```

## Inputs

- **GenBank/DDBJ**: Upload one or more `.gb`/`.gbk` files.
- **GFF3 + FASTA**: Upload matching annotation and sequence files. The `##FASTA` section inside a GFF3 file is not supported; a separate FASTA is required.
- **BLAST (linear mode only)**: Provide outfmt 6 or 7 files. The order of BLAST files must match the order of genome files.

## Editing and Customization

- Click features in the SVG preview to change colors.
- Use the legend editor to rename, reorder, or remove entries.
- Drag and drop the legend and diagram elements to reposition them.
- Adjust canvas padding and other layout settings for spacing.
- Save and load configuration presets for reuse.

## Export

The web app supports:
- **SVG** (recommended for publication and external editing)
- **PNG**
- **PDF**

If you need EPS/PS, export SVG and convert with an external tool such as Inkscape. Mixed-format text may not convert reliably; SVG is the safest format for complex labels.

## Limitations

- Large genomes can hit browser memory limits.
- Initial load can take 5-15 seconds while Pyodide initializes.
- Non-SVG exports are generated without CairoSVG.

## Troubleshooting

- **No similarity ribbons**: confirm BLAST outfmt 6/7 and correct file order.
- **GFF3 input fails**: ensure the FASTA is provided as a separate file.
- **Slow or unresponsive UI**: try a smaller dataset or use the CLI.

[Home](./DOCS.md) | [Installation](./INSTALL.md) | [Quickstart](./QUICKSTART.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Gallery](./GALLERY.md) | [Recipes](./RECIPES.md) | [FAQ](./FAQ.md) | [ABOUT](./ABOUT.md)
