[Home](../README.md) | [Installation](./INSTALL.md) | [Usage](./USAGE.md) | [Gallery](./GALLERY.md) | [RECIPES](./RECIPES.md) | [Advanced](./ADVANCED.md) | [FAQ](./FAQ.md)
# Roadmap

## Planned features
### Short-term plans
- Automatic adjustment of GC content/skew windows and steps (10kb window-1kb step for genomes larger than 1Mb; window size makes a huge difference in PNG/PDF processing time!)
- Tick label configuration (e.g. per Mb)
- `--separate_strands` by default (change in args)
- Legend placement adjustment
- Feature-wise SVG element grouping (for the ease of manual editing)
### Long-term plans
- PAF/MAF support (pairwise matches) 
- Feature overlap resolution (overlapping genes, transcript isoforms etc.;unsure how to impelment mutiple tracks with label overlap resolution)
- Custom tracks (read depth, motifs, etc.)
## Known issues
- **Trans-introns** are not currently visualized.
- **Mixed-format text** (e.g., combining italic and block elements like `<i>Ca.</i> Tyloplasma litorale`) cannot be reliably converted from SVG to PDF/PNG/EPS/PS.  
  → As a workaround, export to **SVG format** and convert to other formats using external tools like [**Inkscape**](https://inkscape.org/).

[Home](../README.md) | [Installation](./INSTALL.md) | [Usage](./USAGE.md) | [Gallery](./GALLERY.md) | [RECIPES](./RECIPES.md) | [Advanced](./ADVANCED.md) | [FAQ](./FAQ.md)