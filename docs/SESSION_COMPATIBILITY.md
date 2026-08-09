[Home](./DOCS.md) | [Current compatibility reference](./REFERENCE/session-and-request-compatibility.md) | [CLI inventory](./CLI_Reference.md) | [Python API](./PYTHON_API.md) | [Typed API](./TYPED_API.md) | **Compatibility history**

# Session and request compatibility history

This page retains the detailed version-by-version migration history for gbdraw
session files, canonical render requests, and saved LOSAT results. The concise
[session and request compatibility reference](./REFERENCE/session-and-request-compatibility.md)
is the canonical owner of current support. Tutorials and the FAQ describe what
a user should do; release notes record when a format changed.

## Supported versions

Current writers emit one session and request format:

| Format | Current writer | Accepted by current readers |
|---|---:|---|
| gbdraw session | 40 | 27–33 and 39–40 |
| Canonical `renderRequest` | 5 | 1, 2, and 5 |

Session versions 34–38 and canonical request schemas 3–4 were development-only
formats. They were never released on the supported history and are rejected.

The public typed-session bridge can convert session versions 31–33 and 39–40 to
a typed request. Versions 27–30 remain supported only as CLI replay inputs
because they do not contain a canonical `renderRequest`. Use the same
`circular` or `linear` subcommand that created the session.

`render_session()` is the compatibility boundary for canonical session replay.
It migrates supported persisted artifacts into `CurrentRequestArtifacts`, then
calls the same current typed renderer used by fresh requests. Fresh
`render_request()` calls never parse old session envelopes or rewrite retired
protein identifiers.

Both `.gbdraw-session.json` and lossless gzip-compressed
`.gbdraw-session.json.gz` files are accepted. The Web app writes the compressed
form by default; the CLI writes uncompressed JSON unless the output name ends in
`.gz`.

## Current request ownership

Canonical request schema 5 records Circular grouping explicitly as `single`,
`grid`, or `batch`. A single diagram or grid has one output object. A Circular
batch has one resolved output object per record. `renderRequest.output.prefix`
is the output-prefix owner.

Current sessions keep mode-specific layout values under
`ui.layoutPreferences`. Supported older parallel title, legend, and grouping
fields are migrated when read and are not written again.

The Web writer stores file bytes once under `resources`. `webFiles` binds those
resources to active and inactive input controls, so files shared by the
committed request and an editable draft are not copied into a second payload.
For Linear comparisons, `webFiles.bindings.linearComparisons` contains only a
stable comparison-edge ID and its file binding. Endpoint, source, inclusion,
and filename metadata are not duplicated there. Version 39 sessions with
legacy embedded `files` remain readable.

`renderRequest` owns the last committed render. Web config retains inactive
Custom Track stacks, disabled rows, draft Axis positions, and per-mode
comparison profiles. `editorState.featureCatalog` holds the schema-3 feature
catalog used by the saved preview and editor. Current sessions store one base
SVG Result per logical diagram; readers collapse paired base and interactive
Results from supported older sessions.

Linear comparison intent is an independent editable draft under
`config.linearComparisonPlan`. Its mode is `none`, `adjacent`, or `selected`;
it also stores the default adjacent source and stable edge metadata. Placement
alone remains under `config.linearRecordLayout`. Current writers do not store a
global `blastSource`, nested layout comparisons, per-record BLAST files, or
per-record LOSAT filenames. The editable plan can therefore opt out without
changing the last committed comparison artifacts in `renderRequest`.

Supported pre-40 Web sessions migrate directly to this plan. A disabled or
absent legacy layout retains the old adjacent LOSAT/upload behavior, an enabled
explicit list becomes `selected`, and an authoritative empty explicit list
becomes `none`. Legacy per-record uploads and custom filenames are attached to
their original positional gap by stable record UID. CLI-only replay sessions
do not gain a synthetic Web comparison draft. The accepted session versions
remain 27–33 and 39–40.

## Retired inputs

Fresh CLI and Python requests reject these retired names or values. Supported
older sessions and canonical request schemas 1–2 migrate them before replay.

| Retired input | Current input |
|---|---|
| Circular `--multi_record_size_mode sqrt` | `--multi_record_size_mode auto` |
| Linear `--label_placement on_feature` | `--label_placement above_feature` |
| Linear `--track_layout spreadout` / `tuckin` | `--track_layout above` / `below` |
| `--depth_tick_interval` | `--depth_large_tick_interval` |
| `--feature_table` | `--feature_visibility_table` |
| `--collinear_max_gene_gap` | `--collinear_max_unit_gap` |
| Circular slot `spacing` | `inner_gap_px` and `outer_gap_px` |
| Circular slot `strict`, `compress`, or `reserve` | No direct replacement; geometry and reservation are derived from `side` |

Current multiword long options use underscore spelling except for the documented
active aliases. `--annotation-table` remains an alias for
`--annotation_table`, and `--gc_content_tick_interval` remains an alias for
`--gc_content_large_tick_interval`.

The private `__gbdraw_legacy_spacing` key is read only from canonical request
schemas 1–2 and is never written by schema 5. Pixel spacing migrates to
`inner_gap_px` and `outer_gap_px`. Factor-based spacing can be replayed but
cannot be saved losslessly in the current format; replace it with explicit pixel
gaps before saving a migrated session.

## Saved protein-comparison results

A protein-search cache hit requires the same amino-acid sequences, selected
protein set, record-instance and feature bindings, query/subject direction,
program, and meaningful search arguments. Upload filenames, modification times,
resource names, and display-only aliases are not part of that identity.

Current sessions use these independent payload schemas:

| Saved payload | Schema |
|---|---:|
| Protein raw search result | 4 |
| Derived protein comparison | 3 |
| Protein identity manifest | 2 |
| Nucleotide raw search result | 2 |

Generated protein FASTA, raw QUERY/SUBJECT fields, protein maps, and derived
references use deterministic session-internal handles of the form
`h_[a-z2-7]{26}`. Those handles bind a record instance to a complete CDS
identity. They are not user-facing protein names.

Sessions 27–33 may contain schema-2 protein candidates and schema-1 derived
evidence. Current readers keep those candidates separate from current cache
entries. Generation verifies the full FASTA content, program and arguments,
direction, and one-to-one feature mapping. A verified candidate is promoted
without rerunning LOSAT; an unverifiable candidate is a cache miss only for its
record pair.

**Save Raw LOSAT TSV** resolves generated protein handles through the identity
manifest immediately before download. It writes readable, percent-encoded
`protein_id`, `locus_tag`, GFF `ID`, or location-based aliases, while preserving
comments, row order, columns 3–12, numeric spelling, and line endings. The whole
download fails if any handle cannot be resolved. User-uploaded comparison TSV is
never rewritten.

## Typed-session resource lifetime

`materialize_session()` decodes embedded resources to temporary paths. Those
paths are valid only inside the active materialization context:

```python
from gbdraw.api import (
    load_session_document,
    materialize_session,
    render_session,
)

document = load_session_document("figure.gbdraw-session.json")
with materialize_session(document, output_directory="out") as materialized:
    result = render_session(materialized)
```

Do not retain a decoded resource path after the `with` block ends. Use
`with_request_output()` inside the same context when replay needs a different
prefix, output directory, format, or overwrite policy.

[Home](./DOCS.md) | [Current compatibility reference](./REFERENCE/session-and-request-compatibility.md) | [CLI inventory](./CLI_Reference.md) | [Python API](./PYTHON_API.md) | [Typed API](./TYPED_API.md) | **Compatibility history**
