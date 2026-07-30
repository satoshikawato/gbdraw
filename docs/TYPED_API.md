[Home](./DOCS.md) | [Python API](./PYTHON_API.md) | **Typed API** | [Session compatibility](./SESSION_COMPATIBILITY.md)

# Build typed render requests

Use the package-root API for ordinary in-memory drawing. The `gbdraw.api`
namespace is for integrations that need explicit input sources, output policy,
request planning, or session round trips.

## Build and render a request

`RecordInput` keeps the source and per-record presentation together.
`RenderOutputRequest` owns the output directory, prefix, formats, and overwrite
policy.

Every `RecordInput` also states how many records may be taken from its source:

- `RecordCardinality.EXACTLY_ONE` is the default. After an optional selector is
  applied, resolving zero or multiple records is an error.
- `RecordCardinality.FIRST` deliberately takes the first candidate in source
  order.
- `RecordCardinality.ALL` expands every candidate in source order.

Selectors and selector-qualified regions identify one record and therefore use
`EXACTLY_ONE`. The planner loads each unique source once per request, then applies
selection, reverse complement, and per-input region transforms. Collection-level
regions, labels, and subtitles are applied after expansion.

```python
from pathlib import Path

from gbdraw.api import (
    CircularDiagramOptions,
    CircularDiagramRequest,
    GenBankInputSource,
    RecordCardinality,
    RecordInput,
    RenderOutputRequest,
    render_request,
)

request = CircularDiagramRequest(
    records=(
        RecordInput(
            source=GenBankInputSource(Path("examples/MjeNMV.gb")),
            cardinality=RecordCardinality.EXACTLY_ONE,
            record_key="reference",
        ),
    ),
    options=CircularDiagramOptions(),
    output=RenderOutputRequest(
        output_prefix="typed-circular",
        output_directory=Path("out"),
        formats=("svg", "interactive_svg"),
        overwrite=False,
    ),
)
result = render_request(request)
print(result.output_paths)
```

Use `GffFastaInputSource` for paired GFF3 and FASTA input, or
`InMemoryRecordSource` for an existing Biopython `SeqRecord`.

## Mode-specific request types

Typed Circular requests use `CircularDiagramOptions`,
`CircularRequestTrackOptions`, and `CircularOutputOptions`. Typed Linear
requests use `LinearDiagramOptions`, `LinearRequestTrackOptions`, and
`LinearOutputOptions`.

`CircularDiagramRequest` represents one `single` or `grid` render.
`CircularBatchRequest` represents separate diagrams. A materialized batch has one
`RenderOutputRequest` per record; an unresolved `ALL` batch instead uses
`CircularBatchOutputPolicy`, which assigns collision-free names after expansion.
`LinearDiagramRequest` represents one Linear layout. Comparison and Circular
conservation tables are likewise resolved only after their record collections have
expanded, so row selectors refer to displayed records.

Call `build_request_diagram()` when an integration needs a validated prepared
render object without writing output. Call `render_request()` to build and save
the requested formats. `CircularRequestPlan`,
`CircularBatchRequestPlan`, and `LinearRequestPlan` expose normalized records,
layout, and builder selection for integrations that need to inspect the plan.
`plan_request()` resolves that plan without assembling a drawing.
`build_request_plan_diagram()` then builds the same plan without loading or
planning it again. This split lets an integration preflight materialized output
paths before any rendering work begins. Each plan exposes
`preflight_outputs()`. Single Circular and Linear plans validate their
materialized output set; a Circular batch validates every resolved target
as one complete set before diagram construction starts. Output generation is
sequential rather than transactional: if a later conversion fails, formats
already written by that item remain.

For Linear requests, `PreparedDiagramRequest.linear_metadata` and
`RequestRenderResult.linear_metadata` use the exported
`LinearDiagramMetadata` type. It carries the computed comparisons,
the compatibility `orthogroups` field containing similarity-group memberships,
and the collinearity result alongside the drawing.

Call `resolve_request()` when an integration needs the planner's canonical,
materialized projection without rendering:

```python
from gbdraw.api import resolve_request

materialized_request = resolve_request(request)
```

The returned request contains one in-memory, `EXACTLY_ONE` input for each displayed
record, resolved table-backed options and outputs, and no collection-level record
transforms. `render_request()` accepts either unresolved or materialized requests and
uses the same planning boundary.

When an integration already owns current saved-analysis artifacts, pass them as
`CurrentRequestArtifacts`. The render boundary accepts only the current raw,
derived, and identity-manifest schemas. Derived runtime handles must resolve through
the supplied manifest. This boundary does not accept a session JSON mapping.
Use `render_session()` for persisted-session replay so supported older artifacts
are migrated by the session compatibility adapter before they enter rendering.

## Typed depth tracks

Add one `DepthTrackInput` for each logical series. `source` accepts one path or
DataFrame shared by all displayed records, or one path, DataFrame, or `None` per
record. Linear entries also accept `height`.

```python
from gbdraw.api import DepthTrackInput, LinearDiagramOptions

typed_options = LinearDiagramOptions(
    depth_tracks=(
        DepthTrackInput(
            source=(
                Path("record-a.depth.tsv"),
                Path("record-b.depth.tsv"),
            ),
            label="Coverage",
            color="#2563eb",
            height=24,
        ),
    ),
)
```

The older `depth_table`, `depth_file`, `depth_tables`, `depth_files`, and
`depth_track_*` fields remain compatibility inputs. Do not combine them with
`depth_tracks`; request normalization rejects mixed forms. Current request and
session writers serialize either accepted input style as canonical
`depthTracks`.

The shorter `gbdraw.api.CircularTrackOptions` and
`gbdraw.api.LinearTrackOptions` names are compatibility aliases for the typed
request classes. They are distinct from the beginner-facing
`gbdraw.CircularTrackOptions` and `gbdraw.LinearTrackOptions`.

## Comparisons and grouping

Optional Circular comparison FASTA paths belong to
`CircularDiagramOptions.conservation_fasta_files`. This keeps interactive
matched-span sequence metadata in the same request as the corresponding BLAST
input.

`RenderOutputRequest.output_prefix` is one filename component. It rejects POSIX
and Windows directory separators, `.` and `..`, ASCII control characters, and
Windows-reserved device, stream, or wildcard names; put the directory in
`output_directory`. Dots inside a valid prefix are preserved: `sample.v1`
produces `sample.v1.svg`. A render always writes the base
`.svg` plus any additional requested formats. Existing regular files are
replaceable only when `overwrite=True`; directories, special files, dangling
symlinks, and invalid parents are rejected even then.

Circular batches require unique resolved output targets. When names are derived
from record IDs, the same one-component restriction applies; use an explicit
batch output prefix for a path-like ID. A batch result returns one item per
record.

For the accepted persisted versions, retired inputs, and migration boundaries,
see [Session and request compatibility](./SESSION_COMPATIBILITY.md).

## Session round trips

The typed session lifecycle is:

1. `build_session_document()` or `save_session_document()` resolves the request and
   embeds its resources.
2. `load_session_document()` validates a saved document.
3. `materialize_session()` provides temporary resource paths.
4. `session_to_request()` converts a canonical session to a typed request, or
   `render_session()` renders it directly.

```python
from gbdraw.api import (
    load_session_document,
    materialize_session,
    session_to_request,
    with_request_output,
)

document = load_session_document("figure.gbdraw-session.json.gz")
with materialize_session(document, output_directory="out") as materialized:
    replay = session_to_request(materialized)
    replay = with_request_output(
        replay,
        output_prefix="replayed",
        output_directory=Path("out"),
        overwrite=False,
    )
    replay_result = render_request(replay)
```

Materialized paths expire when the context closes. The compatibility reference
documents which session versions have a public typed conversion. Calling
`session_to_request()` followed by `render_request()` renders only the decoded
typed request; call `render_session()` when saved analysis artifacts must also be
replayed.

Canonical request schema 5 remains a materialized wire format: it does not encode
runtime cardinality, deferred table paths, record-derived output naming, or
collection-level record transforms. The low-level canonical encoder rejects such an
unresolved request instead of silently dropping fields; call `resolve_request()`
first. Session-document writers do this resolution automatically.

Pin a gbdraw version in reproducible pipelines and compare representative output
after upgrading. SVG geometry may change intentionally while the Python call
remains valid.

[Home](./DOCS.md) | [Python API](./PYTHON_API.md) | **Typed API** | [Session compatibility](./SESSION_COMPATIBILITY.md)
