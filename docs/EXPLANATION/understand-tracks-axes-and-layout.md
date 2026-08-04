[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [How-to guides](../HOW_TO/README.md) | [Reference](../REFERENCE/README.md)

# Understand tracks, axes, and layout ownership

A track slot says what is drawn, where it sits relative to the record axis, and how much space it may use. The axis remains the coordinate reference; moving a track does not change genomic coordinates.

Simple layout controls build a default stack from enabled features, GC content, skew, and depth tracks. Custom slots replace that stack with an explicit ordered model. When a custom stack is active, the saved slot order and axis index own placement, so simple placement controls become inactive. Turning the custom stack off preserves it; Reset rebuilds it.

Circular slots use radial sides and geometry. Slots may sit outside, on the feature axis, or inside. Linear slots use an ordered list split by an axis index; slots before and after that boundary occupy opposite sides, while overlay slots attach to a drawable anchor.

Feature height and reserved slot height are separate. Feature height changes glyph thickness. Slot height reserves space for a renderer and may expand when labels or feature lanes require more room. Spacing is clearance from the adjacent track, not part of the biological coordinate range.

## Outer composition follows visible content

Circular and Linear figures finish record-local layout before placing the legend and plot title. The renderer reports the visible bounds of the primary diagram, including records, tracks, labels, definitions, scales, annotations, and comparisons. One outer composition step then places decorations against that box. Circular grids pack records directly from the same authoritative bounds rather than recovering geometry from serialized SVG.

Docked positions (`left`, `right`, `top`, and `bottom`) use fixed renderer spacing: 24 px from primary content to the legend, 20 px from primary content to the title, 20 px between a legend and title on the same side, and 16 px around the final figure. These are internal defaults, not user settings. They keep the default canvas close to painted content. Titles center on the primary diagram rather than the width after legend docking.

Changing a docked legend or title position changes only the outer transforms and final viewBox. It does not alter Circular radii, record axes, feature labels, track geometry, or comparisons. Corner legend positions remain overlays and use 8 px of obstacle clearance when collision handling shifts them. In the Web editor, current SVGs carry composition metadata schema 1; the editor reapplies automatic placement before a saved drag offset. Supported older SVG results pass through one explicit legacy adapter when loaded.

Numeric tracks own their axes and ticks. GC content may be drawn as deviation from the mean or as percentage; skew depends on the selected dinucleotide. Depth series keep a logical series index even when a record has no values for that series. Empty cells do not reserve painted geometry.

For worked procedures, see [How to add depth, GC content, and skew tracks](../HOW_TO/GUI/add-depth-gc-and-skew-tracks.md) and [How to use input tables](../HOW_TO/CLI/use-input-tables.md). Renderer tokens and constraints are in [Palettes, feature rules, labels, shapes, and track renderers](../REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md).
