[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [Technical documentation](README.md) | [Output formats](output-formats-and-export.md) | [FAQ](../FAQ.md)

# Interactive SVG and semantic hooks

Every gbdraw SVG marks its records, tracks, features, comparison matches,
annotations, and legends with `data-gbdraw-*` attributes. Use these attributes
to select elements in a script, a stylesheet, or a post-processing step. Do not
select by SVG `id`, path order, or path geometry: those can change between
releases even when the figure shows the same features.

## Interactive SVG

Interactive SVG embeds controls, searchable feature metadata, feature and match
popups, group inspection, and the supported sequence downloads. Static SVG keeps
the record, track, feature, match, and annotation attributes below; it has no
`data-gbdraw-interactive-*` markers and no embedded application.

Search in an interactive SVG uses the same fields as
[Feature Search in the Web app](web-app.md#preview-search-and-editor): **All**
does not search nucleotide or amino-acid sequences or `/translation` values,
and **Location** matches the displayed 1-based INSDC location. A feature popup
shows each part of a split or origin-spanning location and their summed length.
A file keeps the interactive runtime it was exported with; export it again to
get the current search and popup behavior.

gbdraw escapes text that comes from your input files. Input text never becomes
an executable `<script>` element or an `on*` event-handler attribute. Species
text accepts `<i>` for italics, as in `<i>Escherichia coli</i>`; it is never
run as HTML.

## SVG IDs

For the same gbdraw version, ordered input, and settings, every SVG `id` is the
same from run to run. IDs are valid XML identifiers and unique within the file,
including Circular multi-record canvases and copied definitions such as clips
and hatch patterns. Every local `href`, `xlink:href`, and `url(#...)` reference
points to one emitted ID.

The exact spelling of an ID can change in any release. Use IDs only for
references inside one file, and select elements by the attributes below.

## Records, tracks, and legends

Record indexes are zero-based and follow the displayed record order. Record IDs
come from the source record ID (`SeqRecord.id`) and need not be unique, so use
the index when two records can share an accession. In Circular multi-record
output, the outer group of each complete record carries both attributes.

| Element | Attributes | Meaning |
|---|---|---|
| Record group | `data-gbdraw-record-id`, `data-gbdraw-record-index` | Source record ID and displayed position |
| Linear record group of a cropped or reverse-complemented record | `data-gbdraw-record-source-start`, `data-gbdraw-record-source-end`, `data-gbdraw-record-source-step` | Input-file span shown by the record (1-based, inclusive) and its direction (`1` or `-1`); record-local base `x` is source `start + x - 1` (step `1`) or `end - x + 1` (step `-1`). Absent when the record shows its input file one to one |
| Record definition | `data-gbdraw-role="record-definition"` or `"record-definition-row"`, `data-gbdraw-definition-part`, record ID/index | Main record text or a row-level definition. A row has `record-definition-row` only when some text describes the whole row, so do not assume one per row |
| Plot title | `data-gbdraw-role="plot-title"` | Shared Circular title |
| Comparison legend | `data-gbdraw-role="comparison-legend"`, `data-gbdraw-orientation` | Identity legend; orientation is `h`, `v`, or `circular` |
| Track group | `data-gbdraw-slot-id`, `data-gbdraw-slot-renderer` | Slot ID (yours in a custom stack, gbdraw's default ID otherwise) and the renderer drawn in that slot |

Typical selectors:

```css
g[data-gbdraw-role="record-definition"][data-gbdraw-record-index="0"]
g[data-gbdraw-slot-renderer="depth"][data-gbdraw-slot-id="coverage"]
[data-gbdraw-role="comparison-legend"][data-gbdraw-orientation="h"]
```

Common renderer values are `features`, `ticks`, `dinucleotide_content`,
`dinucleotide_skew`, `depth`, `annotations`, and `sequence_conservation`. Other
documented track renderers use their renderer name.

## Features, comparison matches, and annotations

| Element | Attributes | Meaning |
|---|---|---|
| Drawn feature part | `data-gbdraw-feature-id`, `data-gbdraw-stable-feature-id`, `data-gbdraw-feature-part`, record ID/index | This drawn element, the biological feature, the kind of part, and the record it belongs to |
| Interactive feature | `data-gbdraw-interactive-feature="true"` | The element has feature metadata in an interactive SVG |
| Comparison match | `data-gbdraw-match-id`; Linear files may also carry `data-gbdraw-pairwise-match-id` | Match ID, unique within the file |
| Interactive match | `data-gbdraw-interactive-match="true"` | The element has match metadata in an interactive SVG |
| Annotation mark | `data-gbdraw-annotation-id`, `data-gbdraw-annotation-set-id`, `data-gbdraw-annotation-track-id`, record index | Annotation, annotation set, track slot, and the record it belongs to |

A feature with a split (joined) location, such as exons separated by introns,
can be drawn as several parts. Use `data-gbdraw-stable-feature-id` for the
biological feature, and `data-gbdraw-feature-id` with `data-gbdraw-feature-part`
for one drawn part. A comparison match carries both its query and subject ends;
an annotation mark carries its set, track, and displayed record. Do not infer a
feature or match from an element's `id` or path geometry.

Other `data-gbdraw-*` attributes support layout, editor state, or the embedded
interactive runtime. They can change in any release.

## Match a feature between renders

To find the same feature in two renders, use the pair `(recordKey,
biologicalFeatureId)`. When a source feature carries it explicitly,
`(recordIndex, stableFeatureId)` or `(recordIndex, sourceFeatureIndex)` also
identifies it. A rendered feature ID names one drawn element: the feature's
source-record hash, also for a cropped or reverse-complemented record, plus a
suffix for the record when several records are drawn, and for identical
features. It changes when the records or their order change. Protein handles
created during a comparison are valid only inside the saved result and protein
identity manifest that created them. `protein_id` and `sourceProteinId` are for
display and export; do not join on them.

Stop rather than guess when an identity is missing, ambiguous, or inconsistent.
Every identity field given for one feature, including a rendered ID, stable
feature ID, and source feature index, must name the same feature in the same
record. Do not fall back to matching a protein label.

[Documentation home](../DOCS.md) | [Technical documentation](README.md) | [Output formats and export](output-formats-and-export.md) | [FAQ](../FAQ.md)
