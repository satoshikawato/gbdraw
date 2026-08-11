[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [Technical documentation](README.md) | [Palette Explorer](../PALETTE_EXPLORER.md) | [FAQ](../FAQ.md)

# Palettes, feature rules, labels, shapes, tracks, and layout

## Diagram layout

The record axis is the coordinate reference. Moving a record, track, title, or
legend changes presentation, not genomic coordinates.

| Layout | Record behavior |
|---|---|
| Circular single record | One complete circular replicon and its radial tracks |
| Circular multi-record | Complete circular records placed in rows on one canvas; size can be automatic, length-scaled, or equal |
| Linear single row | Ordered records or selected regions on one shared bp-per-pixel scale |
| Linear multiple rows | Explicit row placement, with input order preserved within each row |

Circular layout accepts row positions such as `#1@1`; record IDs can replace
indices when order may change. Linear layout uses the same selector-and-row
form. Regions are 1-based and inclusive. Reverse-complement display changes
the shown orientation and endpoint mapping but leaves the source record
unchanged.

Simple controls assemble a track stack from enabled features, coordinate
ticks, GC content, skew, and depth. An explicit slot list is authoritative.
Turning a custom stack off preserves it for later use; resetting settings
rebuilds the simple stack.

Circular slots can sit outside, on, or inside the feature axis. Linear slots
are ordered around an axis boundary, with above, below, and overlay placement
where the renderer supports it. Slot order controls placement. It does not
change which depth series, annotation set, or comparison evidence a renderer
uses.

Coordinate rulers belong to the genome axis. Numeric tracks own their own axes
and ticks. Hiding the coordinate scale therefore leaves the record axis and
quantitative axes available. Plot titles, record definitions, and legends are
placed against visible diagram bounds after record-local layout; moving a
docked decoration does not alter features, tracks, or comparisons.

Outer composition uses the visible bounds of the primary diagram after
record-local layout. Docked placement leaves 24 px between primary content and
a legend, 20 px between primary content and a title, 20 px between a legend and
title on the same side, and 16 px around the final figure. These distances are
renderer defaults, not user settings. Titles center on the primary diagram.
Corner legends remain overlays and use 8 px of obstacle clearance. Current Web
SVGs store composition metadata as schema 1; the editor reapplies automatic
placement before a saved drag offset. Supported older SVGs enter the editor
through one explicit legacy adapter.

The exact surface controls are listed in the [Web app](web-app.md),
[command-line](command-line.md), and [Python API](python-api.md) pages.

## Tracks, axes, and annotations

| Renderer | Content | Main constraint |
|---|---|---|
| `features` | Genomic feature glyphs and labels | Owns the feature axis or overlay anchor. |
| `ticks` | Genome-coordinate ticks and labels | Follows mode and coordinate-scale settings. |
| `dinucleotide_content` | GC percentage, deviation from the mean, or the configured base-pair content | Can own a percentage axis and ticks. |
| `dinucleotide_skew` | Signed skew for the configured base pair | Requires record sequence. |
| `depth` | One logical numeric depth series | Requires a mapped source for each record that has data. |
| `annotations` | One named annotation set | Requires `set_id`; an overlay also requires a drawable anchor. |

Circular slots use side, radius or width, gaps, z-order, and renderer
parameters. Linear slots use order, side or overlay ownership, reserved height,
spacing, and an axis boundary. A feature slot's reserved band and the feature
glyph height are different settings. Invalid anchor chains, missing logical
depth indices, and contradictory overlay order are rejected.

GC-content mode can show deviation from the record mean or an absolute
percentage. Window and step control smoothing and sampling: larger values
produce a smoother, less local trace. Reversing a configured base pair leaves
its content trace unchanged and reverses the sign of its skew. GC-content and
depth axes and ticks are independent of the genome-coordinate scale.

Each repeated depth input is one logical series. It can use one source for all
records or one source, `DataFrame`, or `None` per displayed record. A missing
entry is not converted to zero and does not borrow another record's data. In a
Linear layout it reserves no painted track space for that record, while the
series identity remains stable where data exists. A logarithmic depth axis
requires a positive minimum.

An annotation renderer draws one named `set_id`, not every row in an annotation
table. Coordinate targets and feature targets are mutually exclusive.
Coordinate positions are 1-based and inclusive; only Circular annotations can
cross the origin. Explicit lanes are zero-based. Separate slots allow sets to
use different placement or dimensions, and overlay annotations require a
compatible enabled anchor. See [Input formats and TSV
schemas](input-formats-and-tsv-schemas.md#annotation-table-fields) for the
table fields and mark types.

The default `overflow=error` rejects a track that needs more than its reserved
extent. For line, bracket, and band slots, it also rejects explicit-lane
collisions. A highlight-only slot ignores explicit lanes, flattens its marks
into lane 0, and retains overlaps during lane assignment. `overflow=compress`
permits explicit-lane collisions and, only for an overfull multi-lane track,
reduces the gaps between lanes; it still fails when the fixed lane bands and
padding cannot fit. `overflow=clip` permits explicit-lane collisions and clips
an overfull track. Increase the slot extent or reduce padding and label offset
when the content must remain fully visible.

## Feature presentation

Presentation rules do not edit the input annotation. Their order is:

1. The selected feature types establish the baseline set.
2. A matching specific-color rule or visibility `show` rule can reveal another
   feature; visibility `off` removes it.
3. A matching specific-color rule overrides the selected palette or default
   feature-type color.
4. Qualifier priority chooses candidate label text. Whitelist or blacklist
   filtering runs next, then ordinary label overrides replace selected text.
5. Shape, stroke, overlap, title, definition, and legend settings affect only
   the drawing.

A partial default-color table changes only the listed feature types. Omitted
feature types retain their values from the selected palette.

Whitelist and blacklist label filters are mutually exclusive. An ordinary
label override cannot restore a label removed by the active filter. Exact-match
regular expressions should be anchored when a broader qualifier match would be
unsafe.

Feature visibility actions are `show`, `off`, and `exclude_matching`.
`exclude_matching` keeps the feature visible but removes it from protein-search
input; `off` removes the glyph and the search input. The first matching
visibility rule wins.

Supported feature renderings include directional `arrow`, nondirectional
`rectangle`, and full-band `underlay`. Underlays require an enabled feature
anchor and render behind foreground glyphs. Block strokes, connector strokes,
fill, width, visibility, and overlap resolution are separate controls.
Overlap resolution may add drawing lanes but does not crop a feature or change
comparison evidence. With Circular simple controls, separate-strand placement
owns the lanes and disables overlap resolution.

Trans-introns are not currently visualized.

Web preview edits become exact qualifier rules when one safe rule represents
the chosen scope; otherwise they remain feature-specific rules. The full
editor behavior is in [Preview, search, and
editor](web-app.md#preview-search-and-editor).

The [Palette Explorer](../PALETTE_EXPLORER.md) shows current palette names and
swatches. The schemas for default-color, specific-color, qualifier-priority,
label-filter, label-override, and visibility tables are in [Input formats and
TSV schemas](input-formats-and-tsv-schemas.md#styling-tables).
