[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [How-to guides](../HOW_TO/README.md) | [Reference](../REFERENCE/README.md)

# Understand tracks, axes, and layout ownership

A track slot says what is drawn, where it sits relative to the record axis, and how much space it may use. The axis remains the coordinate reference; moving a track does not change genomic coordinates.

Simple layout controls build a default stack from enabled features, GC content, skew, and depth tracks. Custom slots replace that stack with an explicit ordered model. When a custom stack is active, the saved slot order and axis index own placement, so simple placement controls become inactive. Turning the custom stack off preserves it; Reset rebuilds it.

Circular slots use radial sides and geometry. Slots may sit outside, on the feature axis, or inside. Linear slots use an ordered list split by an axis index; slots before and after that boundary occupy opposite sides, while overlay slots attach to a drawable anchor.

Feature height and reserved slot height are separate. Feature height changes glyph thickness. Slot height reserves space for a renderer and may expand when labels or feature lanes require more room. Spacing is clearance from the adjacent track, not part of the biological coordinate range.

Numeric tracks own their axes and ticks. GC content may be drawn as deviation from the mean or as percentage; skew depends on the selected dinucleotide. Depth series keep a logical series index even when a record has no values for that series. Empty cells do not reserve painted geometry.

For worked procedures, see [How to add depth, GC content, and skew tracks](../HOW_TO/GUI/add-depth-gc-and-skew-tracks.md) and [How to use input tables](../HOW_TO/CLI/use-input-tables.md). Renderer tokens and constraints are in [Palettes, feature rules, labels, shapes, and track renderers](../REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md).
