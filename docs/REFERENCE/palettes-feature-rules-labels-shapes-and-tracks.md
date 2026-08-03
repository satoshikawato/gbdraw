[Documentation home](../DOCS.md) | [Styling how-to](../HOW_TO/GUI/style-features-labels-titles-and-legends.md) | [Palette Explorer](../PALETTE_EXPLORER.md)

# Palettes, feature rules, labels, shapes, and track renderers

## Color precedence

Specific feature-color rules override the selected default-color table for matching features. A palette selects a complete default table. Direct preview edits apply to the current generated feature state and should be saved in a session or exported handoff if they must be reproduced.

Feature rule tables match a feature type and, where configured, qualifier/value conditions. Qualifier-priority tables control which available qualifier supplies a label; they do not change biological feature identity. Label whitelist, blacklist, and override rules are separate filters with explicit precedence.

## Feature presentation

Supported presentation choices include directional arrows, rectangles, and underlays for applicable feature types. Block stroke, connector-line stroke, width, fill, visibility, and overlap resolution remain separate controls. Hiding a feature removes its glyph; it does not alter the underlying input record or comparison evidence.

## Track renderers

| Renderer | Content | Important constraint |
|---|---|---|
| `features` | Genomic feature glyphs and labels | Owns the feature axis/overlay position. |
| `dinucleotide_content` | GC percentage or configured dinucleotide content | May own percentage axis and ticks. |
| `dinucleotide_skew` | Configured dinucleotide skew | Depends on record sequence. |
| `depth` | One logical depth series | Requires a mapped depth source for the record. |
| `annotations` | Named annotation set | Requires `set_id`; overlay annotations require a valid drawable anchor. |
| `ticks` | Coordinate ticks and labels | Placement follows mode and scale settings. |

Circular slots use side, radius/width, gaps, z-order, and renderer parameters. Linear slots use order, side/overlay ownership, reserved height, spacing, and an axis boundary. Invalid annotation anchor chains, missing logical depth indices, or contradictory overlay z-order are rejected.

The palette names and exact swatches are displayed in the [Palette Explorer](../PALETTE_EXPLORER.md). Use the focused [web-app styling guide](../HOW_TO/GUI/style-features-labels-titles-and-legends.md) or [command-line styling guide](../HOW_TO/CLI/set-colors-labels-visibility-shapes-and-strokes.md) for runnable procedures.
