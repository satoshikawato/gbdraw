[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [Technical documentation](README.md) | [Output formats](output-formats-and-export.md) | [SVG hook inventory](../SVG_SEMANTIC_HOOKS.md)

# Interactive SVG and semantic-hook reference

gbdraw SVG groups expose stable semantic identifiers and `data-gbdraw-*` attributes for records, definitions, feature tracks, quantitative tracks, features, comparisons, annotations, ticks, legends, and other rendered roles. Consumers should select by the documented semantic hook, not path order or generated geometry.

Feature identity binds the biological record and feature, including split locations and reverse-complement display. Comparison hooks bind query and subject endpoints. Annotation hooks bind set, annotation, track, and displayed record. Stable hooks persist across rendering changes when the biological identity is unchanged; layout coordinates are not stable API values.

Identity fields have separate roles. Analysis runtime handles are valid only within the
artifact and protein-identity manifest that created them. The canonical cross-render key is
`(recordKey, biologicalFeatureId)`. During normalization, a source feature may instead be
identified by `(recordIndex, stableFeatureId)` or `(recordIndex, sourceFeatureIndex)` when that
source identity is carried explicitly. A rendered feature ID identifies one DOM instance and
may change after reverse-complement display or another rendering transformation. Public
`protein_id` and `sourceProteinId` values are display and export metadata, not join keys.

Consumers must fail closed when a biological identity is missing, ambiguous, or inconsistent.
Every identity field supplied for one endpoint, including a rendered ID, stable feature ID, and
source feature index, must resolve to the same feature in the same record. Do not recover a failed
identity lookup by matching a public protein label.

Interactive SVG adds embedded controls, searchable metadata, feature and match popups, group inspection, and supported sequence downloads. Static SVG may retain semantic groups without embedding the interactive application.

Generated interactive output sanitizes input-derived text and markup. Untrusted input must not become executable `<script>` content or `on*` event-handler attributes. Species text accepts the documented limited markup for display, not arbitrary HTML execution.

The [SVG ID and attribute inventory](../SVG_SEMANTIC_HOOKS.md) lists the exact
integration tokens. This page documents their current meaning; path order and
exact path bytes are not stable hooks.
