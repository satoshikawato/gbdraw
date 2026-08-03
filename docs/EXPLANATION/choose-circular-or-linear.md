[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [How-to guides](../HOW_TO/README.md) | [Reference](../REFERENCE/README.md)

# Choose Circular or Linear

The layout should follow the biological question. Circular diagrams work best when a complete circular replicon is the subject. Linear diagrams make order, orientation, cropped regions, and links between records easier to compare.

| Question | Start with | Reason |
|---|---|---|
| What features and quantitative tracks occur around one complete circular genome? | Circular | The coordinate axis closes at the origin and uses radial space efficiently. |
| How do several chromosomes or plasmids relate as members of one assembly? | Circular multi-record | A shared canvas keeps each replicon circular while preserving record identity. |
| How do homologous spans line up between records? | Linear | Ribbons and curves have explicit record endpoints and direction. |
| Do I need cropped regions, reverse-complement display, rows, or a common bp-per-pixel scale? | Linear | These operations have a direct left-to-right interpretation. |

Circular similarity rings compare one reference against one or more evidence sets. They summarize where evidence lies around the reference, but do not show the same pairwise topology as Linear links. Use [Choose a genome-comparison method](choose-a-genome-comparison-method.md) before selecting a comparison layout.

The layout can be changed later, but placement settings are mode-specific. A Circular grid position does not become a Linear row, and a Linear crop or reverse-complement choice has no Circular equivalent.

Start with the [first Circular tutorial](../TUTORIALS/GUI/first-circular-genome-diagram.md) or [first Linear tutorial](../TUTORIALS/GUI/first-linear-genome-diagram.md). Exact controls and options are listed in the [Web app reference](../REFERENCE/web-app.md) and [Command-line reference](../REFERENCE/command-line.md).
