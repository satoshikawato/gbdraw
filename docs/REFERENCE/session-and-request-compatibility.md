[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [Technical documentation](README.md) | [Compatibility history](../SESSION_COMPATIBILITY.md) | [FAQ](../FAQ.md)

# Session and request compatibility

Current writers emit session version 44 and canonical `renderRequest` schema 8.

| Persisted format | Current writer | Accepted by current readers |
|---|---:|---|
| gbdraw session | 44 | 27–33, 39–42, and 44 |
| Canonical `renderRequest` | 8 | 1, 2, 5, 6, 7, and 8 |
| Web file bindings | 2 | 1; 2 in sessions 41–42 and 44 |

Session 44 documents written with canonical request schema 7 remain readable.
Current saves promote that request to schema 8 while retaining feature
catalog schema 4 and the saved preview until the next Generate.

Session versions 34–38 and request schemas 3–4 were development-only and are
rejected. Do not change a version number, resource hash, or runtime binding by
hand; changing metadata does not migrate its content.

## What a session preserves

A current session embeds its input resources, normalized settings, last
committed render request, generated result, and supported editor and comparison
state. Replay does not depend on the original file path remaining valid.
Treat the session as sensitive when its embedded source data is sensitive.

Both `.gbdraw-session.json` and lossless `.gbdraw-session.json.gz` are
accepted. The web app writes compressed sessions by default. The command line
writes `<output>.gbdraw-session.json` by default. `--session_output` selects an
explicit path and implies saving; a `.gz` suffix selects compression.

The web app's **Save Session** writes the current committed result and editable
state. **Load Session** restores saved values instead of applying fresh browser
defaults. Generate when you want the Result to reflect changed controls; saving before
Generate deliberately preserves the newer draft alongside the earlier Result. Loading and regenerating should preserve biological identities, labels,
record placement, comparison artifacts, and supported editor state; SVG bytes
or text metrics can still differ across gbdraw versions.

Fresh CLI sessions omit `config` because they have no independent Web draft.
Web initializes their settings from `renderRequest` and restores original input
files from their bindings. A present `config` must contain valid `form` and `adv`
objects; a partial draft is rejected. CLI replay preserves a supplied Web draft.

The Web file inventory uses binding schema 2. An explicit composite `c_gb`
restores one editable GenBank File from ordered component resource bindings.
Each component retains its own filename, MIME type, modification time and exact
bytes; duplicate payloads can share storage while retaining every occurrence.
Combining appends LF to each nonempty component that lacks a final LF. Composition
is independent of the last committed `renderRequest`, which remains the replay
authority. Python validates and preserves the draft binding without combining
its components for rendering.

Schema-1 ordinary bindings and File arrays remain supported. Existing sessions
without explicit bindings retain their request-derived source initialization;
original components cannot be recovered if their membership was never saved.
Schema 2 is accepted with sessions 41–42 and 44. Unknown or malformed bindings reject
before import replaces the current work. Older schema-1 readers reject new
schema-2 documents; changing the schema number does not convert them.

## Linear file-level defaults

A Linear source card carries a default **Organism / strain** and **Subtitle /
title** that apply to every record read from that file. A record that leaves its
own field empty inherits the file default; a record that fills it overrides.

`renderRequest.records[].presentation` always carries the resolved text that is
drawn. `webFiles.linearRecordMetadata[]` records the inheritance separately, so
loading restores the same editable state:

| Field | Meaning |
|---|---|
| `fileDefinition`, `fileSubtitle` | the source card's file-level defaults |
| `recordDefinition`, `recordSubtitle` | the record's own value, empty when it inherits |

Both pairs are written only when the source has a file-level default. Sessions
written before these fields existed carry neither; those readers compare the
resolved text against the file default instead, which reads an override that
happens to repeat the default as inheritance. Request schema 8 adds no field for
these Web-only defaults,
because these fields describe Web editing state rather than the
render request: a reader that ignores them still renders identical output.

Session 44 replaces the two editable Linear visibility booleans with independent
selected modes: `linear_accession_visibility` and `linear_length_visibility`,
each set to `auto`, `show`, or `hide`. Auto resolves from the effective rendered
rows and is projected to the existing request-schema-7 booleans; the request and
Python render model did not gain a new field. A version-42 `true` becomes Show,
`false` becomes Hide, and a missing boolean becomes historical Show. If a
selected-mode field is already present, it takes precedence. Current writers do
not write the retired booleans.

When several records share a Linear row, a label or subtitle that no record of
the row contradicts describes the whole row and is drawn once beside it. An
empty value does not contradict anything, so a records table may name a row on
its leading record alone. A value that differs from the one the leading record
carries is record-local, and every record of the row then draws that line above
itself, including the record that leads the row. A subtitle follows its label,
so it is never left beside the row on its own.

Older Linear sessions may contain automatically inferred replicon or organelle
names saved as subtitles. These values remain visible when **Show Replicon** is
off: loading does not guess their origin or remove matching text. Clear an
unwanted record subtitle to restore its file default; clear that default too if
no subtitle is wanted. New uploads leave automatic names to **Show Replicon**.
Loading preserves the saved preview. **Generate** applies the current settings
and definition alignment, so its placement can differ from an older preview.

## Saving settings before loading a source

**Save Session** also works before any biological source is loaded. It preserves
the full editable configuration, Circular and Linear mode profiles, supported
editor preferences, and auxiliary files such as colors, filters and qualifier
priorities. **Load Session** replaces the current Session with those settings,
clearing any previous sources and Result. A rejected Load restores the previous
work. Load a real source and Generate to apply the saved settings to a diagram.

This settings-only variant was introduced in session 42 with explicit
`renderRequest: null`,
empty `results`, and a null feature catalog. It has no committed render. Missing
requests, dangling resources, or biological inputs in either mode cannot select
this variant. Auxiliary files retain their ordinary resource bindings and bytes.

Python can load and materialize a settings-only Session. CLI replay,
`session_to_request()` and `render_session()` report that it has no biological
render request. Existing supported full Sessions remain readable. Current
settings-only writers emit session 44, while current readers also accept the
released session-42 form. Readers whose maximum version is 42
reject newly written session-44 files.

Older settings JSON without a `format` field (containing `form` or `adv`) still
uses the legacy configuration import. It does not need a render request. This
is distinct from a malformed `format: "gbdraw-session", version: 41` envelope.

## Replay boundaries

On the command line, replay a session with the same `circular` or `linear`
subcommand that created it. Output prefix, format, session-output, and overwrite
options may replace their saved counterparts. Other diagram options are
rejected because they would combine persisted and new settings ambiguously.

With `--session_output`, canonical CLI replay writes the regenerated Result and
preserves the editable draft's component bytes, order and File metadata. Resource
IDs may change when the output table is rebuilt. Explicit Web bindings, including
null and empty lists, take precedence over historical direct-source lists.

In Python, `render_session()` is the persisted-session entry point.
`load_session_document()` validates a document, and `materialize_session()`
exposes embedded resources only while its context is active.
`session_to_request()` converts a current materialized session to a typed
request. Rendering that request alone does not replay saved comparison
artifacts; use `render_session()` when those artifacts belong in the result.

`render_request()` accepts current typed requests, not historical session
envelopes. Public typed session conversion accepts full versions 31–33, 39–42, and 44;
versions 27–30 are CLI replay inputs only. Canonical schema 8 retains schema 7's
display values and schema 6's input cardinality, including selectorless `all`
inputs. Resolve a typed request
before encoding when it still contains deferred paths or collection-level transforms.

## Similarity alignment request ownership

For Linear requests, `renderRequest` schema 8 stores `recordTranslations` and
`similarityAlignment` inside `renderRequest.layout`. The nested alignment plan
is schema 2; this does not change the Session version or request schema. Every
translation has one stable `recordKey` and finite `x` and `y` values. An active
plan covers those displayed keys and stores the exact reference, each target's
Select or Skip outcome, requested `preserve` or `match_reference` policy, and
effective source-relative reverse-complement result. Current readers reject
partial, malformed, mismatched, or unsupported plans. Schema 1 of the nested
plan was never released and has no reader. There is no Circular form or generic
transform matrix.

Released request schemas 1, 2, 5, 6, and 7 remain readable. Their
`alignOrthogroupFeature` protein-setting string is confined to a reader-only
legacy path. Current writers emit neither that field nor the old Session-only
`orthogroupState.selectedOrthogroupAlignmentFeature` copy. Released legacy
Sessions with saved feature-catalog and orthogroup identity metadata materialize
that state to the current schema-2 plan before saving. Malformed or unmappable
legacy values produce an actionable error without replacing the last successful
Result. Current requests reject group-only input. Loading a saved preview does
not initialize the diagram Worker or start LOSATP.

A current Session round trip retains exact feature and record identities,
Select/Skip rationale, requested and effective orientation, base translations,
and the immediate pre-align Reset baseline. Ordinary Generate renders the saved
plan through the canonical typed path without guessing a new anchor. A stable
reorder resolves by `recordKey` and biological feature identity. Source
replacement, crop, selector, manual orientation, and record drag clear the plan
with a visible reason. A stale reference requires Reselect/Clear and a stale
target requires Select/Skip; pending or failed repair keeps the last successful
Result. Undo/Redo restores the complete artifact. Preview-only guides,
candidate markers, and recommendation badges are never saved.

In Web **Run Info**, **Source recipe** uses the original input filenames and
public CLI settings. Keep those original files and download any listed generated
helpers with **Download reproducibility files**. That bundle also includes the
canonical session referenced by **Exact replay**, with its embedded resources
and saved analysis artifacts. An unavailable Source recipe reports why it cannot
express the committed semantics losslessly.

Both commands target the successful generation represented by Run Info, even
after Undo/Redo or edits to the current controls. Exact replay reconstructs that
committed generation; it does not apply later preview-only editor changes or an
ungenerated draft. **Save Session** is the route for preserving editable work
alongside the earlier Result. Exact replay does not promise byte-identical SVG
across gbdraw versions or font environments.

## Saved comparison results and cache reuse

Sessions can retain raw comparison rows, derived Similarity groups or
Collinear blocks, and the protein identity information needed to bind them to
features. Every derived runtime handle must resolve to the identity information
stored with that artifact. A missing or inconsistent identity fails instead of
falling back to a display label.

Comparison cache reuse requires the same sequence content, selected proteins,
record and feature bindings, query/subject direction, program, and meaningful
search arguments. Filenames and display labels do not define cache identity.
Only affected record pairs rerun when one valid cache key changes.

Pairwise hit limits, Similarity-group member limits, and Collinear block
settings are derived options. Changing one recomputes the affected derived
result while retaining eligible raw search rows. Current derived artifacts
record their upstream raw keys and requested settings. Collinear provenance
also records the effective `cds` or `locus` unit kinds produced by `auto` and
whether orthogroup inference was enabled. Sessions preserve each LOSATP mode's
hit-limit drafts. A saved Collinear pipeline without `collinearInferOrthogroups`
uses its historical inference behavior (ON); new Web state defaults to OFF.

**Save Raw LOSAT TSV** resolves generated protein handles to stable readable
aliases. Uploaded comparison TSV is never rewritten. Export raw results for a
durable evidence record; the cache exists to avoid repeated work.

For the version-by-version record of retired input names and saved-result
formats, see the [compatibility history](../SESSION_COMPATIBILITY.md). Release
notes record when support changed; this page documents current support.

## Record rotation and feature placement

Session 41 and request schema 7 add requested record display and feature placement
intent together. Schema 6 keeps its original cardinality and row-inheritance
meaning; session 40 retains its committed-request and editable-config authority.

Each schema-7 record has `display: {isCircular, startCoordinate}`. Both values
are nullable. A null start and an explicit source coordinate 1 remain distinct.
`diagramOptions.featurePlacements` contains sorted exact record/biological-feature
identities with a Main or mode-compatible lane-1 target. Auto removes an override.
Tolerance belongs to `canvas.feature_overlap_tolerance_bp`, a non-negative integer
with default 0. Resolved lanes, pixel coordinates and display fragments are not
requested persistence fields.

Supported older requests have unset display, empty placements and tolerance 0.
Saving a schema-5/6 Web session promotes it without requiring Generate. Versions
27–30 still support CLI replay only, and unknown or development-only versions
remain rejected.

Editable Web rotation and placement drafts are saved in config, separately from
the last successful request and Result. An inactive rotation start stays in the
draft but is omitted from the effective request. Load restores both states;
Generate applies the draft. Source replacement invalidates the replaced source's
bindings, even when its filename or input-card UID is unchanged.
