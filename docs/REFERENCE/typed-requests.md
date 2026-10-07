[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [Technical documentation](README.md) | [FAQ](../FAQ.md)

# Typed request reference

The `gbdraw.api` namespace exposes explicit input, planning, rendering, output, analysis-artifact, table, track, and session types and functions for pipelines and integrations. Use the package-root [Python API](python-api.md) for ordinary in-memory drawing.

## Request models

| Type | Represents |
|---|---|
| `CircularDiagramRequest` | one Circular result in `single` or `grid` grouping |
| `CircularBatchRequest` | one resolved Circular output per selected record |
| `LinearDiagramRequest` | one ordered Linear layout and its comparison plan |
| `RecordInput` | source, cardinality, selector, region, orientation, stable key, and presentation for one declaration |
| `RenderOutputRequest` | output directory, one-component prefix, formats, and overwrite policy |
| `CircularBatchOutputPolicy` | collision-free output naming after a batch expands |

`DiagramRequest` is the union of the three request forms.

## Input sources and cardinality

`GenBankInputSource` holds one GenBank path. `GffFastaInputSource` holds a matched GFF3 and FASTA pair. `InMemoryRecordSource` holds one or more Biopython `SeqRecord` values.

Every `RecordInput` declares a `RecordCardinality`:

| Value | Behavior |
|---|---|
| `EXACTLY_ONE` | default; zero or multiple selected records is an error |
| `FIRST` | deliberately selects the first source candidate |
| `ALL` | expands all candidates in source order |

Selectors and selector-qualified regions identify one record and therefore require `EXACTLY_ONE`. The planner loads each unique source once, then applies selection, reverse complement, and record-local region transforms. Collection-level ordering and layout are applied after expansion.

For Linear Similarity Group alignment, `similarity_alignment` stores the exact
reference and each record's anchor or Skip decision. It has no orientation
field. A `SimilarityAlignmentReference` input is replaced by that plan when the
planner resolves it. Each `RecordInput` sets its orientation through
`RecordPresentation.reverse_complement` or its region setting. The planner
projects anchor centers after resolving those record transforms.

## Planning and rendering lifecycle

| Function | Result and side effects |
|---|---|
| `resolve_request(request)` | materialized request; no drawing or file output |
| `plan_request(request)` | `CircularRequestPlan`, `CircularBatchRequestPlan`, or `LinearRequestPlan`; no diagram assembly |
| `build_request_plan_diagram(plan)` | builds an already resolved plan without loading or planning again |
| `build_request_diagram(request)` | validates and builds a prepared diagram without writing output |
| `render_request(request)` | resolves, plans, builds, and writes requested formats |
| `render_prepared_request(prepared, output)` | writes a previously prepared request |

Each plan exposes `preflight_outputs()`. A single Circular or Linear plan validates its materialized output set; a Circular batch validates all resolved targets together before diagram construction. Format generation is sequential, not transactional: a later conversion failure does not remove formats already written.

`resolve_request()` returns one in-memory, `EXACTLY_ONE` input per displayed record, resolved table-backed options and outputs, and no collection-level record transforms. `render_request()` accepts unresolved or already materialized requests through the same planning boundary. A request with a `SimilarityAlignmentReference` is the one exception to "no drawing": planning it runs the requested orthogroup analysis once, in memory, and the planned request carries the resolved plan and that analysis as precomputed comparisons.

## Output rules

`RenderOutputRequest.output_prefix` is one filename component, not a path. It rejects POSIX and Windows separators, `.` and `..`, ASCII control characters, and Windows-reserved device, stream, and wildcard names. Put directories in `output_directory`. Dots inside a valid prefix are preserved. It is at most 200 bytes in UTF-8. A rejected prefix raises `INPUT_INVALID` with `field: output_prefix` and reason `REQUIRED`, `FILENAME` or `FILENAME_LENGTH`. The web app sends the typed **Output Prefix** as typed; a prefix it derives from a record ID replaces the characters a browser replaces in a download name (`"~*/:<>?\|` and control characters) with `_`, drops leading and trailing dots and spaces, and puts `_` before a Windows device name, so a Result is saved under its own name.

A render always writes the base `.svg` plus requested additional formats. Existing regular files are replaceable only with `overwrite=True`; directories, special files, dangling symlinks, invalid parents, and target collisions are errors. Circular batches require unique resolved targets.

## Mode-specific options and metadata

Circular requests use `CircularDiagramOptions`, `CircularRequestTrackOptions`, `CircularMultiRecordOptions`, and `CircularOutputOptions`. Linear requests use `LinearDiagramOptions`, `LinearRequestTrackOptions`, `LinearMultiRecordOptions`, and `LinearOutputOptions`.

The shorter `gbdraw.api.CircularTrackOptions` and `gbdraw.api.LinearTrackOptions` names are compatibility aliases for the typed request track classes. They are not the beginner-facing `gbdraw.CircularTrackOptions` and `gbdraw.LinearTrackOptions` classes.

`PreparedDiagramRequest.linear_metadata` and `RequestRenderResult.linear_metadata` use `LinearDiagramMetadata`. It contains computed comparisons, the compatibility `orthogroups` field for Similarity-group membership, and the collinearity result.

### Ortholog path representation

Protein comparison inference returns `OrthogroupGraphResult`. Its
`path_indexes_by_orthogroup_id` maps group IDs to `OrthologPathCollection` values.
Each collection provides an exact integer `count`, one-based `path_at(rank)`,
`path_by_id(id)`, `rank_of(protein_ids)`, and `containing_count(protein_id)`.
Normal analysis, metadata, and saved resources do not expand all paths.

`materialize_ortholog_paths(result)` explicitly returns the original
`OrthogroupResult` with tuple-valued `ortholog_paths_by_orthogroup_id`. The protein
comparison producers also accept `path_representation="exhaustive"`.
`collection.iter_paths()` explicitly traverses every path. All three operations
take time proportional to their output; retaining the complete output can require
exponential memory. Collections deliberately have no implicit iterator or `len`.

The legacy `OrthogroupResult` constructor and supplied path tuples remain
supported. Saving or reading a legacy corpus preserves its exact contents,
order, IDs, and shared information; it does not infer additional paths from its
edges. Graph payloads store decimal count strings to avoid JavaScript integer
rounding. Typed analysis resource schema 3 is required to read newly saved data.

## Resolved Similarity Group alignment

A Linear request may carry `similarity_alignment=SimilarityAlignmentPlan(...)`
and a complete `LinearMultiRecordOptions.record_translations` sequence keyed by
stable `recordKey`, or `similarity_alignment=SimilarityAlignmentReference(...)`,
which the planner resolves to such a plan after the orthogroup analysis. The schema-2 plan records the exact reference, one
anchor or Skip decision per record, and rationale. Record presentation or region
state sets the orientation; the plan does not store direction settings. It must be
fully resolved; a group-ID string, partial record coverage, or schema-1 plan
is invalid in a current request. See the [executable typed Python
example](python-api.md#typed-linear-similarity-group-alignment) and the
[Session and request compatibility](session-and-request-compatibility.md#similarity-alignment-in-requests-and-sessions).

## Depth tracks

One `DepthTrackInput` represents one logical series. `source` accepts a path or `DataFrame` shared by all displayed records, or one path, `DataFrame`, or `None` per record. Linear entries may set `height`.

The legacy `depth_table`, `depth_file`, `depth_tables`, `depth_files`, and `depth_track_*` inputs remain compatibility inputs. Do not combine them with `depth_tracks`; normalization rejects mixed forms. Current request and session writers serialize an accepted form as `depthTracks`.

## Current analysis artifacts

Pass `CurrentRequestArtifacts` when an integration already has current raw comparison results, derived grouping/block results, or the protein identity manifest. The fresh render boundary accepts only the current artifact schemas, and every derived runtime handle must resolve through the supplied manifest.

`CurrentRequestArtifacts` does not accept an arbitrary session JSON mapping.
Use `render_session()` for a supported saved session and
`CurrentRequestArtifacts` only for current analysis-artifact models.

## Session lifecycle

| Function | Purpose |
|---|---|
| `build_session_document()` | resolve one request, or several drawings, and embed their resources in a session document |
| `save_session_document()` | build and write the document |
| `load_session_document()` | parse and validate a saved document |
| `upgrade_session_document()` | return a Session 31–44 in the current version without rendering it, with a warning for each dropped Result |
| `materialize_session()` | expose embedded resources as temporary paths |
| `session_to_request()` | convert one drawing of a materialized session to a typed request |
| `with_request_output()` | replace output settings without mutating the request |
| `render_session()` | migrate supported persisted state and replay one drawing's request plus saved analysis artifacts |
| `render_session_drawings()` | render several drawings together, with distinct output names |
| `derive_region_drawing()` | derive a new drawing of selected regions from a materialized session |

`build_session_document()` and `save_session_document()` accept optional
`title` and `created_at` values. If `created_at` is omitted, the writer records
the current UTC save time. A fixed timestamp can make a test fixture
byte-reproducible; it does not make a replay universally reproducible.

Materialized paths expire when the materialization context closes. `session_to_request()` followed by `render_request()` renders only the decoded typed request; use `render_session()` when saved comparison artifacts must also be replayed.

`derive_region_drawing(materialized, regions)` returns a `RegionDrawing`. Its
source is the session's only drawing, or the drawing that `drawing=` names by
ID or name, as in `session_to_request()`. Each
`RegionSelection(record_key, start, end)` names a record of the source drawing
and a region in source coordinates (1-based, inclusive); give at most one
region per record. `drawing` in the result is a `SessionDrawingSpec` for
`build_session_document(drawings=...)`, and `drawing.request` draws the same
files cut to those regions, in source order. The current Session version names
each drawing by its mode, so it cannot write a `name` other than `Linear` or
`Circular`. It keeps the drawing's settings, colors and rules,
label tables, and the per-feature edits, placements, annotations, and Depth
inside the regions:

- `margin` adds bases on each side, clamped to the record ends.
- The mode is Linear, or the source's mode for one whole record. `mode` can
  choose Circular for one record; a Circular region is drawn as a closed
  circle numbered from 1.
- Each record keeps its orientation unless `reverse_complement` is given.
- A region whose start is after its end (one that crosses the origin of a
  circular record) is refused.
- A size you set (label font, stroke and axis widths, feature height, track
  widths, windows, tick interval, tick fonts) is kept only when its Auto value
  is the same at the new length and mode. Otherwise it returns to Auto and is
  listed in `adaptation.reset` with both Auto values. Pass `adapt_sizes=False`
  to keep every size.
- Comparison tables and ring tables are not carried, because they use the
  coordinates of the whole records. A LOSAT search runs again on the new
  records. `dropped` lists every setting or item that was not carried.

The request reads the materialized files, so build or render it before the
materialization context closes.

Session conversion rejects values from the wrong mode. For example, a Circular
request containing Linear track values raises `SessionConversionError`.

### Drawings

A Session (project) holds one or more drawings. `SessionDocument.drawings`
lists them in document order as `SessionDrawing(id, name, mode,
has_canonical_request)`; a `SessionDrawing` is a drawing of the project, not an
SVG drawing such as `RequestRenderResult.drawing`. A drawing without a
committed render, such as the one drawing of a settings-only Session, has
`has_canonical_request=False`.

`SessionDocument.drawing(selector, mode=...)` selects a drawing by its exact ID,
else by a unique exact name; without a selector it returns the only drawing, or
the only drawing of `mode`. `SessionDocument.active_drawing_id` names the
drawing the Web app opens. `SessionDocument.mode` and
`SessionDocument.has_canonical_request` describe the only drawing. With several
drawings, they and an unselected or unknown selection raise
`SessionDrawingSelectionError`, whose message lists each drawing's ID, mode and
name.

Pass `drawing=` (an ID or a name) to `session_to_request()` and
`render_session()` to choose the drawing. `render_session_drawings()` renders
every drawing that has a committed render, or the drawings named in
`drawings=`; it skips the others with a logged notice, and naming one of them
raises `SessionDrawingSelectionError`. One drawing keeps its output names.
Several drawings write `<base>_<id>`, where `<base>` is `output_prefix` or the
drawing's own prefix, and a Circular batch inside still appends `_<n>`. Every
output path of every selected drawing is checked before the first file is
written, and each embedded resource is parsed once for all drawings. The
result maps drawing IDs to render results in document order.

The current Session version, 46, holds at most one drawing of each mode: a
Circular drawing with ID `circular` and name `Circular`, and a Linear drawing
with ID `linear` and name `Linear`. The Web app writes the second drawing's
Result set in `otherModeResult`. Each drawing keeps its settings and Legend
edits in its mode's slice of `modes`; both drawings share the LOSAT caches.
`build_session_document(drawings=[...])` takes typed requests or
`SessionDrawingSpec(request, mode=, id=, name=, state=)` values, where `state`
holds Web-owned drawing fields such as `results`, `editorState`, or the
drawing's slice in `modes`; `active_drawing=` names the drawing to open. It
rejects what version 46 cannot hold: two drawings of one mode, a second drawing
without its Results, other IDs or names, or a second drawing whose shared
fields differ from the first's.

`upgrade_session_document()` returns a `SessionUpgrade`: the current
`document` and its `warnings`. A current document is returned unchanged. A
Session 31–44 gets the migrations of a CLI re-save without a render: its
request is decoded, adapted to current typed state, and encoded again with the
same resource IDs, and its Web-owned fields are migrated. Results with a
feature catalog (Sessions 40–44) are kept; Sessions 31–39 saved none, so their
Results are dropped until the next render. A drawing whose Results are dropped
gets one line in `warnings` that names each dropped Result, and the line is
logged as a warning. Render the drawing and save the Session to write new
Results. Sessions 27–30 have no canonical request and raise
`SessionVersionError`.

Request schemas 6 and 7 record each input's runtime cardinality. This lets a
selectorless source retain `RecordCardinality.ALL` until record planning expands
it. Deferred table paths, record-derived output naming, and collection-level
transforms still require `resolve_request()` before encoding. Session writers
perform that resolution automatically.

## Other exported names

`gbdraw.api` also exports table readers and row models, record and region selectors, annotation and track models, request plans and render results, output-byte helpers, current web-runtime capability constants, and the session exception hierarchy. The full list of exports is `gbdraw.api.__all__`; names outside that list are not part of this public namespace.

## Related

- [Tutorials](../TUTORIALS/README.md)
- [Python API reference](python-api.md)
- [Session and request compatibility](session-and-request-compatibility.md)
- [Input formats and TSV schemas](input-formats-and-tsv-schemas.md)
- [Output format and export reference](output-formats-and-export.md)

## Feature identity overrides

Per-feature edits name one original-source feature by its request record key
and biological feature ID. The ID is the original-coordinate hash, with
`~<n>` added for identical features of one record, so an edit keeps naming the
same feature after crop, reverse complement, reordering, and record
duplication. A record input loaded twice has two record keys, so each copy has
its own edits.

`CircularDiagramOptions` and `LinearDiagramOptions` accept exactly one of
`feature_placements`, `feature_placement_table`, and
`feature_placement_table_file`. Exact overrides use the public
`FeaturePlacementOverride(record_key, biological_feature_id, target)` and
`FeaturePlacementTarget(kind, side=None, level=None)` types. Main uses
`kind="main"`; directional lane 1 uses `kind="lane"` with a mode-compatible
side and `level=1`. Final feature-slot geometry determines direction support.

`feature_overrides` holds `FeatureOverride(record_key, biological_feature_id,
feature_visibility=None, label_visibility=None, label_text=None)` rows
(`diagramOptions.featureOverrides`, request schema 9). `None` keeps the
rule-based result, a row must set at least one value, and an identity may appear
once. Instead of rows, pass a `feature_override_table` DataFrame or a
`feature_override_table_file` path with the
[feature override table](input-formats-and-tsv-schemas.md#feature-override-table)
columns; the three inputs are mutually exclusive. A row decides before the
feature visibility and label tables:

| Field | Values and effect |
|---|---|
| `feature_visibility` | `on` draws the feature, `off` hides it, and `exclude_matching` ignores the visibility table and keeps the feature-type selection. `off` and `exclude_matching` also leave LOSATP protein extraction. |
| `label_visibility` | `off` hides the label. `on` shows it regardless of label scope, whitelist, and blacklist, using `label_text`, else the qualifier text with the non-`hash` label rules applied, else `<type> <start>..<end>`. |
| `label_text` | One line of text. With `label_visibility` unset, it replaces the text of a label that is shown anyway and never shows a label. |

`FeatureIdentitySpan(record_key, biological_feature_id, envelope, circular_path)`
is an annotation target for one such feature. A feature that is not drawn skips
the annotation with the `feature_selector_unmatched` warning.

The shared planner turns both tables into exact rows when the records load,
before rendering or encoding, and resolves every identity once. A
table row must name a feature of the source. A record key outside the request is
an error.
An edit whose feature a cropped record does not have (`crop_excluded`: outside
the crop, or removed while loading), that an uncropped record does not have
(`absent`, for example removed by a GFF3 type filter), or that is not in the
source (`unresolved`) does not fail the render: the edit stays dormant and is
reported in `feature_identity_notices` on the render result and on `Diagram`, in
the Web metadata, and as one CLI log line per notice. Placement applies the same
rule. A GFF3 input also loads the type of each feature whose Feature visibility a
row sets, so a feature a row turns `on` is drawn and one it turns `off` stays in
the Web feature catalog.
Requested placement can be combined with each `RecordInput.display` and the
`canvas.feature_overlap_tolerance_bp` config override. It does not change source
coordinates, sequences, or feature identities.

## Record display intent

`RecordInput.display` is `RecordDisplayOptions(is_circular=None,
start_coordinate=None)`. Set it on each selected biological record. The planner
validates topology, length, and crop constraints after source resolution, then
passes the same display transform to both diagram modes. An unset start retains
the existing reverse-complement display; explicit 1 anchors source base 1.
Resolved plans expose where each record came from and its transforms without rotating the
source `SeqRecord`. Reordering records does not change their exact feature
placement identities.

## Combined typed request example

In a new directory, save the complete [combined Python example](python-api.md#combined-rotation-and-placement-example)
as `rotated_placed_chloroplast.py` and obtain its four inputs. The following
program runs that recipe, then expresses the same presentation as a typed
request. It writes `typed_rotated_placed_chloroplast.svg` and checks whole-SVG
parity with the package-root API. Save it as `typed_chloroplast.py` and run
`python typed_chloroplast.py`.

<!-- executable:joint-typed:start -->
```python
import runpy
from pathlib import Path
from xml.etree import ElementTree
from gbdraw.api import (
    CircularDiagramRequest, CircularDiagramOptions, CircularOutputOptions,
    CircularRequestTrackOptions, ColorOptions, InMemoryRecordSource,
    RecordInput, RecordDisplayOptions, RenderOutputRequest, render_request,
)

recipe = runpy.run_path("rotated_placed_chloroplast.py")
options = recipe["options"]
request = CircularDiagramRequest(
    records=(RecordInput(InMemoryRecordSource(recipe["record"]),
        display=RecordDisplayOptions(start_coordinate=5500)),),
    options=CircularDiagramOptions(
        config_overrides=options.config_overrides,
        colors=ColorOptions(color_table_file=str(options.features.color_table)),
        selected_features_set=options.features.types,
        feature_placement_table=options.features.placements,
        qualifier_priority_file=str(options.labels.qualifier_priority),
        annotations=options.annotations,
        tracks=CircularRequestTrackOptions(circular_track_slots=recipe["track_slots"]),
        output=CircularOutputOptions(legend=options.legend),
        species=options.species,
    ),
    output=RenderOutputRequest(output_prefix="typed_rotated_placed_chloroplast"),
)
result = render_request(request)
assert ElementTree.tostring(ElementTree.parse("typed_rotated_placed_chloroplast.svg").getroot()) == ElementTree.tostring(ElementTree.fromstring(recipe["chloroplast_bytes"]))
print("Typed and package-root SVG trees are identical")
```
<!-- executable:joint-typed:end -->
