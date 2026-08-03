[Documentation home](../DOCS.md) | [Web export how-to](../HOW_TO/GUI/export-publication-and-interactive-figures.md) | [CLI export how-to](../HOW_TO/CLI/export-static-and-interactive-outputs.md) | [Publication guidance](../EXPLANATION/prepare-publication-and-reproducible-figures.md)

# Output format and export reference

| Format | Web app | CLI | Python | Main use |
|---|---:|---:|---:|---|
| SVG | Yes | Yes | Yes | Editable vector master and web display |
| Interactive SVG | Yes | Yes | Yes | Browser search and popups |
| PNG | Yes | Yes | Yes | Fixed-pixel raster output |
| PDF | Yes | Yes | Yes | Vector print and submission |
| EPS | No | Yes | Yes | Requested legacy vector workflow |
| PS | No | Yes | Yes | Requested legacy PostScript workflow |

SVG is the base render. CLI and Python conversion to PNG, PDF, EPS, or PS requires CairoSVG and its runtime dependencies. The browser uses its packaged export path and validates the downloaded file independently.

Static SVG uses `<prefix>.svg`. Interactive SVG uses `<prefix>.interactive.svg`. A session uses `.gbdraw-session.json` or `.gbdraw-session.json.gz`. Output prefixes are filename components, not paths. Existing files are replaced only when overwrite is explicitly enabled; directories, special files, unsafe parents, and invalid prefix components remain errors.

DPI affects a raster only together with physical dimensions. Compute pixels from final inches multiplied by requested DPI. Check output dimensions and file signatures rather than trusting an extension.

Interactive SVG contains controls and metadata that static SVG does not. It must remain free of unsafe input-derived scripts or event-handler attributes. Downstream SVG optimizers may strip required IDs or data attributes. See [Interactive SVG and semantic hooks](interactive-svg-and-semantic-hooks.md).
