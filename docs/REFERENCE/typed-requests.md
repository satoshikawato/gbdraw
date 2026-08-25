[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [Technical documentation](README.md) | [FAQ](../FAQ.md)

# Typed request reference

The `gbdraw.api` namespace exposes explicit input, planning, rendering, output, analysis-artifact, table, track, and session contracts for pipelines and integrations. Use the package-root [Python API](python-api.md) for ordinary in-memory drawing.

## Request models

| Type | Ownership |
|---|---|
| `CircularDiagramRequest` | one Circular result in `single` or `grid` grouping |
| `CircularBatchRequest` | one resolved Circular output per selected record |
| `LinearDiagramRequest` | one ordered Linear layout and its comparison plan |
| `RecordInput` | source, cardinality, selector, region, orientation, stable key, and presentation for one declaration |
| `RenderOutputRequest` | output directory, one-component prefix, formats, and overwrite policy |
| `CircularBatchOutputPolicy` | collision-free output naming after a batch expands |

`DiagramRequest` is the union of the three request forms.

## Input sources and cardinality

`GenBankInputSource` owns one GenBank path. `GffFastaInputSource` owns a matched GFF3 and FASTA pair. `InMemoryRecordSource` owns one or more Biopython `SeqRecord` values.

Every `RecordInput` declares a `RecordCardinality`:

| Value | Behavior |
|---|---|
| `EXACTLY_ONE` | default; zero or multiple selected records is an error |
| `FIRST` | deliberately selects the first source candidate |
| `ALL` | expands all candidates in source order |

Selectors and selector-qualified regions identify one record and therefore require `EXACTLY_ONE`. The planner loads each unique source once, then applies selection, reverse complement, and record-local region transforms. Collection-level ordering and layout are applied after expansion.

## Planning and rendering lifecycle

| Function | Result and side effects |
|---|---|
| `resolve_request(request)` | canonical materialized request; no drawing or file output |
| `plan_request(request)` | `CircularRequestPlan`, `CircularBatchRequestPlan`, or `LinearRequestPlan`; no diagram assembly |
| `build_request_plan_diagram(plan)` | builds an already resolved plan without loading or planning again |
| `build_request_diagram(request)` | validates and builds a prepared diagram without writing output |
| `render_request(request)` | resolves, plans, builds, and writes requested formats |
| `render_prepared_request(prepared, output)` | writes a previously prepared request |

Each plan exposes `preflight_outputs()`. A single Circular or Linear plan validates its materialized output set; a Circular batch validates all resolved targets together before diagram construction. Format generation is sequential, not transactional: a later conversion failure does not remove formats already written.

`resolve_request()` returns one in-memory, `EXACTLY_ONE` input per displayed record, resolved table-backed options and outputs, and no collection-level record transforms. `render_request()` accepts unresolved or already materialized requests through the same planning boundary.

## Output rules

`RenderOutputRequest.output_prefix` is one filename component, not a path. It rejects POSIX and Windows separators, `.` and `..`, ASCII control characters, and Windows-reserved device, stream, and wildcard names. Put directories in `output_directory`. Dots inside a valid prefix are preserved.

A render always writes the base `.svg` plus requested additional formats. Existing regular files are replaceable only with `overwrite=True`; directories, special files, dangling symlinks, invalid parents, and target collisions are errors. Circular batches require unique resolved targets.

## Mode-specific options and metadata

Circular requests use `CircularDiagramOptions`, `CircularRequestTrackOptions`, `CircularMultiRecordOptions`, and `CircularOutputOptions`. Linear requests use `LinearDiagramOptions`, `LinearRequestTrackOptions`, `LinearMultiRecordOptions`, and `LinearOutputOptions`.

The shorter `gbdraw.api.CircularTrackOptions` and `gbdraw.api.LinearTrackOptions` names are compatibility aliases for the typed request track classes. They are not the beginner-facing `gbdraw.CircularTrackOptions` and `gbdraw.LinearTrackOptions` classes.

`PreparedDiagramRequest.linear_metadata` and `RequestRenderResult.linear_metadata` use `LinearDiagramMetadata`. It contains computed comparisons, the compatibility `orthogroups` field for Similarity-group membership, and the collinearity result.

## Depth tracks

One `DepthTrackInput` represents one logical series. `source` accepts a path or `DataFrame` shared by all displayed records, or one path, `DataFrame`, or `None` per record. Linear entries may set `height`.

The legacy `depth_table`, `depth_file`, `depth_tables`, `depth_files`, and `depth_track_*` inputs remain compatibility inputs. Do not combine them with `depth_tracks`; normalization rejects mixed forms. Current request and session writers serialize an accepted form as canonical `depthTracks`.

## Current analysis artifacts

Pass `CurrentRequestArtifacts` when an integration already owns current raw comparison results, derived grouping/block results, or the protein identity manifest. The fresh render boundary accepts only the current artifact schemas, and every derived runtime handle must resolve through the supplied manifest.

`CurrentRequestArtifacts` does not accept an arbitrary session JSON mapping.
Use `render_session()` for a supported saved session and
`CurrentRequestArtifacts` only for current analysis-artifact models.

## Session lifecycle

| Function | Purpose |
|---|---|
| `build_session_document()` | resolve a request and embed its resources in a session document |
| `save_session_document()` | build and write the document |
| `load_session_document()` | parse and validate a saved document |
| `materialize_session()` | expose embedded resources as temporary paths |
| `session_to_request()` | convert a canonical materialized session to a typed request |
| `with_request_output()` | replace output settings without mutating the request |
| `render_session()` | migrate supported persisted state and replay the request plus saved analysis artifacts |

`build_session_document()` and `save_session_document()` accept optional
`title` and `created_at` values. If `created_at` is omitted, the writer records
the current UTC save time. A fixed timestamp can make a test fixture
byte-reproducible; it does not make a replay universally reproducible.

Materialized paths expire when the materialization context closes. `session_to_request()` followed by `render_request()` renders only the decoded typed request; use `render_session()` when saved comparison artifacts must also be replayed.

Session conversion rejects values from the wrong mode. For example, a Circular
request containing Linear track values raises `SessionConversionError`.

Canonical request schema 6 records each input's runtime cardinality. This lets a
selectorless source retain `RecordCardinality.ALL` until record planning expands
it. Deferred table paths, record-derived output naming, and collection-level
transforms still require `resolve_request()` before encoding. Session writers
perform that resolution automatically.

## Exported supporting contracts

`gbdraw.api` also exports table readers and row models, record and region selectors, annotation and track models, request plans and render results, output-byte helpers, current web-runtime capability constants, and the session exception hierarchy. The authoritative export inventory is `gbdraw.api.__all__`; names outside that list are not part of this public namespace contract.

## Related

- [Python Tutorials](../TUTORIALS/PYTHON/README.md)
- [Python API reference](python-api.md)
- [Session and request compatibility](session-and-request-compatibility.md)
- [Input formats and TSV schemas](input-formats-and-tsv-schemas.md)
- [Output format and export reference](output-formats-and-export.md)
