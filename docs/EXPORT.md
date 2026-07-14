[Home](./DOCS.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Workflow guide](./WORKFLOW_GUIDE.md) | **Export** | [CLI Reference](./CLI_Reference.md)

# Export figures for publication

Choose the final format from the journal or production workflow, not from the screen preview.

## Format chooser

| Format | Use it for | Notes |
|---|---|---|
| SVG | Editing, web display, archival vector source | Preserves paths, groups, text, and metadata. |
| PDF | Submission and print workflows that accept vector figures | Check font substitution in the final PDF viewer. |
| PNG | Raster-only systems, slides, or a fixed pixel deliverable | Set dimensions from final print size and required DPI. |
| EPS/PS | Legacy publication workflows | Use only when requested; test text and transparency conversion. |
| Interactive SVG | Feature/match inspection in a browser | Not a substitute for a static submission figure. Keep the session too. |

SVG is always available. PDF, PNG, EPS, and PS require CairoSVG; install the export extra when needed:

```bash
python -m pip install "gbdraw[export]"
gbdraw circular --gbk genome.gb -o figure -f svg,pdf,png
```

Use `interactive_svg` to create `figure.interactive.svg` in addition to the base SVG.

## PNG size and DPI

DPI only has meaning with a physical size. Compute pixels as:

```text
pixels = inches × DPI
```

A 180 mm-wide figure is about 7.09 inches. At 300 DPI it needs about 2,127 pixels; at 600 DPI it needs about 4,254 pixels. Do not label a small raster as high-DPI without increasing its pixel dimensions.

## Fonts and text

- Keep the fonts used during generation available on the system that converts or edits the vector file.
- Check italic markup such as `<i>Genus species</i>` after conversion.
- Embed fonts where the journal workflow permits it, or convert text to paths only in a final derivative. Keep an editable SVG with real text.
- Reopen PDF/EPS output on another machine to catch fallback fonts and changed line breaks.

## Editing vector output

Edit a copy of the SVG. Preserve element IDs, data attributes, groups, and embedded metadata when downstream interactivity or provenance depends on them. Flattening, optimizing, or round-tripping through a vector editor can remove popup metadata even when the static figure still looks correct.

## Accessibility and print checks

- Do not rely on hue alone; combine color with labels, order, stroke, or orientation where possible.
- Check common color-vision deficiencies and a grayscale preview.
- Verify text size at final placement, not at browser zoom.
- Print or rasterize a proof to catch thin strokes and low-contrast legend entries.
- Record the palette and all manual color overrides with the figure.

For reproducibility, retain the source input identifiers/checksums, complete command or session, gbdraw version, editable SVG, and the exact submitted derivative.

[Home](./DOCS.md) | [Tutorials](./TUTORIALS/TUTORIALS.md) | [Workflow guide](./WORKFLOW_GUIDE.md) | **Export**
