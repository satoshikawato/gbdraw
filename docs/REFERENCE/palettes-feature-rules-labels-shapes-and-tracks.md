[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [Technical documentation](README.md) | [Palette Explorer](https://gbdraw.app/gallery/palettes/) | [FAQ](../FAQ.md)

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
ticks, GC content, skew, and depth. An explicit slot list overrides these simple controls.
Turning a custom stack off preserves it for later use; resetting settings
rebuilds the simple stack.

Circular slots can sit outside, on, or inside the feature axis. Linear slots
are ordered around an axis boundary, with above, below, and overlay placement
where the renderer supports it. Slot order controls placement. It does not
change which depth series, annotation set, or comparison results a renderer
uses.

Coordinate rulers belong to the genome axis. Numeric tracks draw their own axes
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
placement before a saved drag offset. The editor converts supported older SVGs when it opens them.

The exact surface controls are listed in the [Web app](web-app.md),
[command-line](command-line.md), and [Python API](python-api.md) pages.

## Tracks, axes, and annotations

| Renderer | Content | Main constraint |
|---|---|---|
| `features` | Genomic feature glyphs and labels | Provides the feature axis or overlay anchor. |
| `ticks` | Genome-coordinate ticks and labels | Follows mode and coordinate-scale settings. |
| `dinucleotide_content` | GC percentage, deviation from the mean, or the configured base-pair content | Can draw its own percentage axis and ticks. |
| `dinucleotide_skew` | Signed skew for the configured base pair | Requires record sequence. |
| `depth` | One logical numeric depth series | Requires a mapped source for each record that has data. |
| `annotations` | One named annotation set | Requires `set_id`; an overlay also requires a drawable anchor. |

Circular slots use side, radius or width, gaps, z-order, and renderer
parameters. In the Web app, Width and Radius use numeric fields with px/×R
selectors; R is the base circle radius. A unit change preserves the number,
and blank means Auto. An inside `gc_content`, `gc_skew`, or `depth` row with a
blank Radius is placed as in the default stack of the selected preset: at the
preset's radius for that row when it fits, otherwise packed below the row
above it, so a custom stack that repeats the default stack draws the same
figure. A slot with a radius is drawn at that radius, and the
other rows on its side keep the stack order around it: rows before it are
placed farther from the center, rows after it closer. When the inside rows do
not fit, numeric rows with an Auto width are compressed together as in an Auto
stack, a row with a radius included (centred on it); an explicit width is kept. Two adjacent
rows on one side keep the larger of their facing gaps, the inner row's outer
gap (`outer_gap_px`) or the outer row's inner gap (`inner_gap_px`); a gap that
is not set may shrink with compressed rows, to no less than 1 px. See [Circular track Width and Radius](web-app.md#circular-track-width-and-radius)
for input, History, Session and Generate behavior. Linear slots use order,
side or overlay placement, reserved height,
spacing, and an axis boundary. A feature slot's reserved band and the feature
glyph height are different settings. Invalid anchor chains, missing logical
depth indices, and contradictory overlay order are rejected.

Circular inside tracks stack outward from the center definition. If they do not
fit only because of the definition, the species line is wrapped at word
boundaries into the fewest lines that let them fit, and the tracks are placed
again; the text is unchanged. The wrap is not applied to a figure that already
fits or to an explicit **Center Reserved Radius** (`--center_reserved_radius`).
`definition_font_size` 18 is the default, so 18 counts as not explicit and may
wrap; any other value is explicit and never wraps. If the tracks still do
not fit, the error names the slot that could not be placed and the reserved
definition radius. Shorten the species or strain text, reduce the definition
font size, set a smaller center reserved radius, or place tracks outside. When
an explicit center reserved radius is the limit, the error names that radius
instead; set a smaller radius or place tracks outside. When GC content or GC
skew tracks lie inside the slot that could not be placed, the error also names
the room kept for them and offers turning them off or narrowing them
(`--no-gc`, `--no-skew`, `--gc_content_width`, `--gc_skew_width`; in the web
app, **Hide GC Content**, **Hide GC Skew**, or their widths).

GC-content mode can show deviation from the record mean or an absolute
percentage. Window and step control smoothing and sampling: larger values
produce a smoother, less local trace. Reversing a configured base pair leaves
its content trace unchanged and reverses the sign of its skew. GC-content and
depth axes and ticks are independent of the genome-coordinate scale.

Each repeated depth input is one logical series. It can use one source for all
records or one source, `DataFrame`, or `None` per displayed record. A missing
entry is not converted to zero and does not borrow another record's data. In a
Linear layout it reserves no painted track space for that record, while the
series identity remains stable where data exists. On a Circular multi-record
canvas the record keeps an empty depth slot; without any depth input no depth
slot is reserved. A logarithmic depth axis requires a positive minimum.

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
into lane 0, and keeps overlaps during lane assignment. `overflow=compress`
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

When several specific-color rules match one feature, the first match wins.
Rules are tried by qualifier key, in the order the keys appear on the feature,
and then by row order within one key; the row order of rules with different
qualifier keys has no effect. The Web app's **Move rule up** and **Move rule
down** therefore change the result only between rules with the same qualifier
key, and the CLI orders the rules the same way.

A partial default-color table changes only the listed feature types. Omitted
feature types keep their values from the selected palette.

The generated feature legend reflects which specific color rules actually match
features in the diagram. Default-colored CDS are captioned `CDS` when no specific
CDS rule is used. When a specific CDS rule is used and default-colored CDS remain,
the legend contains the specific caption and `other proteins` for the remainder.
An unused rule does not change `CDS`; removing the used rule and regenerating
restores it. Other feature types use `other <feature type>s` for the same remainder
(for example, `other genes`). The same rule applies to Web regeneration after
a **This feature only** color edit creates a specific rule: the default caption
changes because its category now contains only the remaining members.

Legend rows that have one name follow one rule in the CLI, the Python API, and
the Web app:

- Rows with the same name and the same swatch (fill and line) are one legend
  row.
- When two or more rows with the same name and different swatches are drawn,
  every one of them is numbered in legend order, for example `Packaging (1)`
  and `Packaging (2)`. A name drawn once has no number.

The rule covers every row: specific rules that give one caption several
colors, a used specific rule whose caption equals a generated row (a feature
type such as `tRNA`, an `other <feature type>s` row, or a track row), and
track and region-annotation rows with one `legend_label`. Two depth tracks
labelled `Coverage` draw one row when their colors are equal, and
`Coverage (1)` and `Coverage (2)` when they differ. Release 0.13.0 drew
`Coverage` and `Coverage (2)` in both cases. Two GC content tracks always have
the palette's GC content color, so they draw one `GC content` row.

The numbers belong to the drawn legend only. They follow the legend order, so
they change when a row is added, removed, or moved, and a name left with one
row loses its number. gbdraw does not rewrite a rule's caption: the color
table, the Web app's rule list and popup, the downloaded color-table TSV, and
a saved Session keep it as written, and no legend name carries a color ID.
One compatibility reading applies when a color table file or a saved Session
is read: a caption that ends in its own rule's color ID, such as
`Packaging [#ffe700]` on a rule colored `#ffe700`, loses that ID, also when
you typed that name yourself. The ID is the form older builds saved: a space,
then the rule's color as lower-case `#rrggbb` (or `none`) in square brackets.
Another spelling of the color (`[#FFE700]`, a color name), the ID of another
color, and a caption that is only the ID stay as written. A caption typed in
the Web app is drawn as typed until its Session is saved and loaded, and a
color table passed to the Python API as a `DataFrame` is used as written.

Whitelist and blacklist label filters are mutually exclusive. An ordinary
label override cannot restore a label removed by the active filter. Exact-match
regular expressions should be anchored when a broader qualifier match would be
unsafe.

A per-feature override is a label-override row whose qualifier is `hash`. It
decides one feature's label before the label filters and the label display
scope (`--show_labels` or `--labels`, `labels.linear.scope`,
`labels.circular.scope`). Non-empty text shows the label, also on a record or
with a scope that otherwise draws no labels. Empty text hides it. No label
is drawn in three cases: the feature is hidden, the feature is drawn as
`underlay`, or, with `labels.rendering = "embedded_only"`, the label does not
fit inside its feature. Such a row shows a label only once one of these
changes. The Web
feature popup writes these rows for **Label visibility** **On** and **Off**. A
label text edit uses `locus_tag`, `gene`, or `record_location` when one of them
identifies the feature, and falls back to `hash` only when none does.

Feature visibility actions are `show`, `off`, and `exclude_matching`.
`exclude_matching` keeps the feature visible but removes it from protein-search
input; `off` removes the glyph, its label, and the search input. The first matching
visibility rule wins.

A [feature override table](input-formats-and-tsv-schemas.md#feature-override-table)
row (`--feature_override_table`, or typed `feature_overrides`) names one feature
by its source identity and decides its Feature visibility, Label visibility,
and label text before the visibility and label tables.

Supported feature renderings include directional `arrow`, nondirectional
`rectangle`, and full-band `underlay`. Underlays require an enabled feature
anchor and render behind foreground glyphs. Block strokes, connector strokes,
fill, width, visibility, and overlap resolution are separate controls.
Overlap resolution may add drawing lanes but does not crop a feature or change
comparison results. Separate-strand placement and overlap resolution can be enabled together in
both modes. Undefined-strand features share the nominal negative pool with
negative-strand features during automatic allocation.

Trans-introns are not currently visualized.

Web preview edits become exact qualifier rules when one safe rule represents
the chosen scope; otherwise they remain feature-specific rules. The full
editor behavior is in [Preview, search, and
editor](web-app.md#preview-search-and-editor).

The [Palette Explorer](https://gbdraw.app/gallery/palettes/) applies any
built-in palette to one fixed Circular diagram, so only the colors change. The
palette tokens and semantic colors come from
[`gbdraw/data/color_palettes.toml`](../../gbdraw/data/color_palettes.toml), and
the generated [color palette reference](../../examples/color_palette_examples.md)
lists every palette with representative Circular and Linear figures.

Palette names describe appearance, not a guaranteed accessibility grade. Check
contrast against the final background and at final print size. Do not rely on
color alone for a biological distinction; combine it with labels, order, track
position, feature shape, stroke, or another redundant cue.

The schemas for default-color, specific-color, qualifier-priority,
label-filter, label-override, and visibility tables are in [Input formats and
TSV schemas](input-formats-and-tsv-schemas.md#styling-tables).

## Labels that do not fit

After a figure is laid out, gbdraw checks every drawn label: Circular outer,
inner, and embedded labels in either placement, and Linear labels. Two
horizontal Circular labels overlap when their boxes come within 1 px
horizontally and 0.1 px vertically; radial, embedded, and Linear labels
overlap when their text areas intersect. Two leader lines cross when they
intersect; lines that only touch do not cross. A Circular label with a leader
covers the definition text in the centre when its box overlaps a line of that
text by more than 0.5 px in both directions.

Before the check, gbdraw repairs small defects of horizontal Circular labels.
When a record has 1 to 100 overlapping label pairs, crossing leader pairs, and
labels over the definition text in total, counted over its outer, inner, and
embedded labels, gbdraw moves outer and inner labels by up to 24 px, bends
their leader lines by up to 12 px, and may move the point where a leader meets
its label; embedded labels do not move. One move can shift up to 24 following
labels of the same side, or a run of up to 12 neighbouring labels that reaches
4 labels beyond the defects. Labels keep their side and order. The repaired
layout is kept only if it leaves no defect; otherwise the record keeps its
layout and the rules below apply. Records without a defect, and records with
radial labels, are not changed.

After a successful repair, the labels around the repaired defects make room:
the labels of the defects and up to 12 neighbouring labels on each side of
them, on the same side of the circle, may move within the same limits toward
shorter, more radial leader lines, and labels next to a crowded run may move
away from it along the circle. Each leader in the window meets its label at the
point its direction of arrival selects. This step never adds a defect or a
leader line through another label's text; if it would, the repaired layout is
kept. It can enlarge the canvas by up to about 24 px and adds about 1 to 8 s to
the drawing of a repaired map. Records the repair does not change are not
affected.

- If any labels overlap, or radial labels need more room than the circular
  tracks leave, nothing is drawn. The command line prints `ERROR:` with the
  message and exits with status 1. Python raises `LabelsDoNotFitError`, a
  `GbdrawError`.
- If leader lines cross, labels cover the definition text, or Linear region
  (annotation) labels overlap other labels or lines, but no labels overlap,
  the figure is drawn as before. The command line prints the same
  message after `WARNING:` and exits with status 0. Python warns with
  `LabelLeadersCrossWarning`, a `UserWarning`.

The message counts the overlapping label pairs, crossing leader pairs,
labels over the definition text, and region label overlaps of each record and
side, and lists settings
that make the labels fit. Each setting changes one thing and is checked by
drawing the whole figure again. The rows depend on which labels do not fit:

| Row | Change | Offered for |
|---|---|---|
| outward, both directions | `--outer_label_x_radius_offset` and `--outer_label_y_radius_offset` together, in steps of 0.025, up to 2 | outer labels |
| outward, vertical only | `--outer_label_y_radius_offset`, in steps of 0.025, up to 2 | outer labels |
| outward, horizontal only | `--outer_label_x_radius_offset`, in steps of 0.025, up to 2 | outer labels |
| smaller font | `--label_font_size`, in steps of 0.5 pt, down to 6 pt | all labels |
| tighter rows | `--circular_label_spacing`: 2.5, 2, 1.5, 1, 0.5, 0.25, or 0.1 px, below the current value | outer labels |
| combination | the font reduced by half the smaller-font change (at least 0.5 pt), with the smallest Y offset that then fits | outer labels |
| all labels outside | `--labels out` | records with inner labels, in either placement |
| radial labels | `--label_placement radial` | inner and embedded labels |
| separate overlapping features | `--resolve_overlaps` | labels drawn on features that are drawn on top of each other, Circular and Linear |

Each row gives the mildest value that fits, or says that none fits within its
limit. When nothing was drawn and a row has no value that fits, it gives the
mildest value that leaves no overlapping labels, if there is one, and says how
many leader pairs still cross or labels still cover the definition text. Such
a row comes after the others; running again with it draws the figure with the
warning. With several records, one value must fit every record. The outer-label
rows do not move inner, embedded, or radial labels. For inner, embedded, or radial
labels, the switches are tried first; the smaller font is then the first
font size, going down from the current one, at which the whole figure fits, so
a smaller size than the one offered may not fit.

When the figure is drawn anyway, gbdraw searches for 5 seconds and names the
rows it did not reach. When nothing is drawn, every outer-label row is
searched; for inner, embedded, or radial labels the switches and then the font
search share 30 seconds, and a setting is tried only if drawing the figure once
more ends within that time. Linear feature labels get only the separate
overlapping features row, because the Linear layout keeps its other labels
apart. Linear region labels get no row.

A Linear region label moves along its region to a place clear of the feature
labels, leader lines, record names and definition text, other records of its
row, and the region labels placed before it. When no place along its region
is clear, gbdraw tries, in this order:

1. the empty area inside the region's own filled box (a `highlight` or `band`
   with a solid fill), clear of feature glyphs, the axis, and other tracks by
   the label spacing;
2. at the label's place along its region, each feature label whose letters
   it meets moves to another row on its side: the nearest row first, then
   one new row beyond the others, for which the figure is drawn again;
3. a label that then meets only the leader lines of its record's feature
   labels is drawn in front of them on a white rounded box at 70% opacity,
   and a leader line behind the box does not count.

Steps 2 and 3 may combine. gbdraw uses a step only when it clears every
conflict of the label and makes no new overlap, checked by drawing the
figure; otherwise the label stays on one line along its region. The font
size and the number of lines never change. A figure in which every region
label has a clear place along its region is drawn as before. When a conflict
remains, each thing the label still overlaps counts as one region label
overlap (a feature label and its leader line count once). This only warns:
the Web app shows "Region labels overlap other labels or lines (N)." over
the figure, with no **Choose a fix** button. SVG, PNG, and PDF keep the box's
opacity and drawing order; EPS and PS keep them by rasterizing the area of
the box.

In the Web app, the rows appear in the **Labels do not fit** dialog. Generate
searches no row of a figure drawn anyway; the dialog computes each row that
Generate did not reach, one after another, for at most 30 seconds each. A row
that finds no value in that time shows "No value found within 30 seconds."
The next time the dialog opens, that row is searched again for up to 30
seconds, continuing from where it stopped, so it can still find a value.
While a row is being computed, Apply and Cancel first wait for it to finish,
at most 30 seconds, and the dialog shows "Finishing the choice being
computed…".

With 50 or more outer labels, gbdraw packs the labels in rows with
`--circular_label_spacing` between them. If those rows do not fit, the labels
are packed 0.1 px apart, as in 0.13, so a map that fit in 0.13 still draws.

In Python, `LabelsDoNotFitError.report` and `LabelLeadersCrossWarning.report`
hold a `LabelFitReport`: the records that do not fit, by label kind
(`circular_outer`, `circular_inner`, `circular_embedded`, `linear`, or
`linear_region` for Linear region labels), with
their label counts, font sizes, and counts, and `refused` for radial labels
without room; and the rows. A row's `values` maps config paths to the values
that fit, and its `status` is `found`, `found_with_notice` (the figure then
draws with the warning; `remaining` holds what is left), `not_found`, or
`not_computed`.
`RequestRenderResult.label_fit_report` returns the report of a figure drawn
anyway.

## Manual feature placement

A placement belongs to a whole biological feature within its existing feature
slot. Every exon, connector and fragment created by rotation uses one assignment.
Glyphs, labels, leaders, reserved space and canvas bounds follow that assignment.
Feature-associated comparison ribbons attach to the corresponding endpoint
feature's paint with 4 SVG units of clearance; a record's other lanes do not
replace the endpoint identity.

| Resolved feature layout | Available requests |
|---|---|
| Circular split, combined or separate strands | Auto, Main, Outward lane 1, Inward lane 1 |
| Linear overlay, combined strands | Auto, Main, Above lane 1, Below lane 1 |
| One-sided, either mode; Linear separate strands | Auto, Main |

Main means the nominal lane in the feature's strand pool. In Circular split
layout, Outward/Inward lane 1 sits one lane beyond the nominal outer/inner lane
when strands are separate, or one lane from the shared Main lane when combined.
Moving one feature leaves the other features' nominal lanes in place.
Custom slots determine
availability from their final geometry; a preset name cannot override that
geometry. An unsupported direction fails with a reason that names the feature
and the lanes that accept it; Auto and Main are available in every slot.

| Request | Resolver off | Resolver on |
|---|---|---|
| Auto (no override) | Nominal lane; overlaps allowed | Allocated deterministically after fixed features |
| Main or supported lane 1 | Fixed placement | Fixed placement reserved first |

Two fixed features conflict when their envelope intersection on the same physical
lane exceeds `feature_overlap_tolerance_bp`. The tolerance is a non-negative
integer number of bases, default 0; touching intervals do not conflict. Circular
wrap intersections are summed. An intron remains inside its multipart parent's
collision envelope. Conflicts identify the features and lane; they do not silently
move either fixed feature.

In the Web Feature Editor choose **Feature placement**, or use **Selected feature
placements** for a selection. Auto deletes the override; Undo/Redo restores the
requested state. Apply rotation, placement or tolerance changes with **Generate
Diagram**. They are drafts until generation succeeds; the previous Result remains
available. A source replacement breaks its link to the old source. Known invisible
features can be dormant; an unknown or stale feature cannot. Placements stay
with the mode in which they were made: a Generate in the other mode ignores
them, and they apply again after you switch back.

Lane 2 and higher, arbitrary pixel dragging, per-exon placement, annotation-lane
placement and a general label-collision redesign are outside this feature.
