[Documentation](./DOCS.md) | [Interactive SVG tutorial](./TUTORIALS/8_Interactive_SVG_Sessions.md) | [Export](./EXPORT.md)

# SVG IDs and semantic hooks

gbdraw exposes semantic `data-gbdraw-*` attributes for first-party interaction
and targeted post-processing. Use these attributes instead of matching internal
SVG `id` prefixes.

## ID guarantee

For the same gbdraw version, ordered input, and resolved configuration, every
emitted SVG `id` is deterministic. IDs are valid XML/SVG identifiers and are
unique within the document, including Circular multi-record canvases and copied
definitions such as clips and hatch patterns. Every local `href`, `xlink:href`,
and `url(#...)` reference resolves to one emitted ID.

Exact ID spelling is not a cross-version API. A release may change an internal
ID while preserving the semantic attributes below. Static SVG consumers should
therefore rely on IDs only for document-local references.

## Supported structural hooks

Record indexes are zero-based and refer to displayed-record order. Record IDs
contain the source `SeqRecord.id` and need not be unique, so use the index when
duplicate accessions are possible.
In Circular multi-record output, the outer wrapper for each complete record
carries both attributes.

| Element | Supported attributes | Meaning |
|---|---|---|
| Record group | `data-gbdraw-record-id`, `data-gbdraw-record-index` | Source record identity and displayed instance |
| Record definition | `data-gbdraw-role="record-definition"` or `"record-definition-row"`, `data-gbdraw-definition-part`, record ID/index | Main record text or a row-level definition |
| Plot title | `data-gbdraw-role="plot-title"` | Shared Circular title |
| Comparison legend | `data-gbdraw-role="comparison-legend"`, `data-gbdraw-orientation` | Identity legend; orientation is `h`, `v`, or `circular` |
| Track group | `data-gbdraw-slot-id`, `data-gbdraw-slot-renderer` | Raw logical slot ID and the renderer selected for that slot |

Typical selectors are:

```css
g[data-gbdraw-role="record-definition"][data-gbdraw-record-index="0"]
g[data-gbdraw-slot-renderer="depth"][data-gbdraw-slot-id="coverage"]
[data-gbdraw-role="comparison-legend"][data-gbdraw-orientation="h"]
```

The common renderer values include `features`, `ticks`,
`dinucleotide_content`, `dinucleotide_skew`, `depth`, `annotations`, and
`sequence_conservation`. Additional documented track renderers use their
canonical renderer name.

## Supported feature, match, and annotation hooks

| Element | Supported attributes | Meaning |
|---|---|---|
| Rendered feature part | `data-gbdraw-feature-id`, `data-gbdraw-stable-feature-id`, `data-gbdraw-feature-part`, record ID/index | Rendered instance, biological feature identity, part kind, and owning record |
| Interactive feature | `data-gbdraw-interactive-feature="true"` | The element has an interactive metadata entry |
| Comparison match | `data-gbdraw-match-id`; Linear files may also carry `data-gbdraw-pairwise-match-id` | Stable document-local match identity |
| Interactive match | `data-gbdraw-interactive-match="true"` | The element has an interactive match payload |
| Annotation mark | `data-gbdraw-annotation-id`, `data-gbdraw-annotation-set-id`, `data-gbdraw-annotation-track-id`, record index | Annotation, set, slot, and owning record |

Compound features can produce several rendered parts. Use
`data-gbdraw-stable-feature-id` for biological identity and
`data-gbdraw-feature-id` plus `data-gbdraw-feature-part` for a particular
rendered instance. Do not infer feature or match identity from the element's
`id` or path geometry.

Other `data-gbdraw-*` attributes may support layout bookkeeping, transient
editor state, or the embedded standalone runtime. They are not part of this
semantic selector contract unless listed above.

[Documentation](./DOCS.md) | [Interactive SVG tutorial](./TUTORIALS/8_Interactive_SVG_Sessions.md) | [Export](./EXPORT.md)
