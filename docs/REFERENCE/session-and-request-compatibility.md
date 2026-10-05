[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [Technical documentation](README.md) | [Compatibility history](../SESSION_COMPATIBILITY.md) | [FAQ](../FAQ.md)

# Session and request compatibility

Current writers emit session version 45 and canonical `renderRequest` schema 9.

| Persisted format | Current writer | Accepted by current readers |
|---|---:|---|
| gbdraw session | 45 | 27–33, 39–42, 44, and 45 |
| Canonical `renderRequest` | 9 | 1, 2, 5, 6, 7, 8, and 9 |
| Web file bindings | 2 | 1; 2 in sessions 41–42, 44, and 45 |

Session 44 documents written with canonical request schema 7 or 8 remain
readable. Current saves promote that request to schema 9 and the Session to 45,
and keep the saved preview until the next Generate. Schema 9 adds
`diagramOptions.featureOverrides` (see
[Feature identity overrides](typed-requests.md#feature-identity-overrides));
a schema-8 request reads as one without overrides.

Session 45 stores the Web app's per-feature edits (Feature visibility, Label
visibility, and label text) in `features.featureOverrides`, one row per
original-source feature named by `recordKey` and `biologicalFeatureId`, as the
request does, and by the mode the edit belongs to (`scope`); Feature placement
drafts name their mode the same way. Both modes can use the same record key, so
a request carries only the rows of its own mode. Its feature catalog is schema 5, which records the selector
values each drawn feature had (`drawnSelector`). Loading a Session 44, or a
Web Session 31–33, moves its edits keyed by rendered feature ID onto the
feature they name: a rendered ID in the saved feature catalog names its
feature; otherwise a copy or record suffix is removed when exactly one feature
remains. A Session before 40 is matched through its GenBank sources read again
with its crops and orientations (the drawn hash and record position of each
rendered ID), or, without readable sources, through its saved feature metadata
for records drawn without a crop or reverse complement. An edit that names no
feature is dropped, and the Web app reports how many. An older Feature
visibility edit that hid every feature with its hash (each copy of a duplicated
record) now applies to the edited feature only; the Web app reports how many.
Each moved edit takes the mode of the Session's diagram. A Feature placement
draft of a Session 41–44 takes the mode of its lane side; a Main placement,
which reached requests of both modes, is kept for both. A catalog of schema 3 or
4 reads as schema 5 without selector values until the next Generate; a feature
whose rendered ID carries its source hash was drawn with its source
coordinates, so its source values serve until then.

Session versions 34–38 and request schemas 3–4 were development-only and are
rejected. Do not change a version number, resource hash, or runtime binding by
hand; changing metadata does not migrate its content.

## What a session preserves

A current session embeds its input resources, normalized settings, last
committed render request, generated result, and supported editor and comparison
state. Replay does not depend on the original file path remaining valid.
Treat the session as sensitive when its embedded source data is sensitive.

Replay reads embedded BLAST outfmt 6/7 resources with the current [comparison
table rules](input-formats-and-tsv-schemas.md#comparison-and-numeric-tables).
A saved table with more than 12 columns replays from its first 12 columns. A
malformed saved table, or one whose IDs name the wrong Linear records, stops
replay with an error.

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

A Web Session also preserves supported manual legend, plot title, and Linear
scale positions in its saved Result. After loading, **Generate Diagram** carries
each matched item's offset into the new automatic layout for the same diagram;
it does not pin an absolute page position. If the source, selected region,
record set, or layout no longer matches a moved item, Generate keeps the saved
Result and asks for that item's position reset, **Reset Layout**, or matching
settings. Preview search placement and Layout edit hints are not saved. Keep
the Session when browser positioning must be reproduced; a raw Python render
request alone does not contain those manual positions.

Fresh CLI sessions omit `config` because they have no independent Web draft.
Web initializes their settings from `renderRequest` and restores original input
files from their bindings. When the request draws every record of an input file,
its record source names that file's resource, so feature-popup record rotation
works right after Load. A Linear request, or a single Circular request, that
draws only some records of a multi-record file (for example with `--record_id`
or a records table) also reads that file and selects each drawn record by
record ID, or by `#n` when another record of the file has the same ID; the
file stays one File that draws only those records. A Linear record, or the
record of a single Circular request, that is drawn cropped or
reverse-complemented stores its crop in source coordinates as `region` and its
orientation as `region.reverseComplement` or `presentation.reverseComplement`,
as Web Save does; a crop taken after
`--reverse_complement` is written as the same span in source coordinates. Each
Linear row takes its crop and orientation from the request, so **Generate
Diagram** draws the record the CLI drew. A cropped or reversed record of a
Circular batch or grid, the records of a file that a Circular batch or grid
draws only partly, and every record of an older CLI session keep their own
drawn copy; rotating such a record needs **Generate Diagram** first. A present
`config` must contain valid `form` and `adv` objects; a partial draft is
rejected. CLI replay preserves a supplied Web draft.
A CLI binding uid such as `cli-seq-1` is only an initial value: each Linear
file takes the record identity of `renderRequest.records[].recordKey`. A
multi-record file, whose records are `record-1:1`, `record-1:2`, and so on,
becomes one Linear row per record with that recordKey and the request's record
selector (the `#n` of its recordKey when an older CLI session stored no
selector), so **Inherit saved comparison** finds every record it names. A CLI
comparison stays read-only in Web; **Inherit saved comparison** reuses it after
promoting an older committed request, such as a version 42 sidecar, to the
current schema. A Linear Session without a stored comparison plan loads with
**No comparison**. A CLI Session written with `-b` therefore offers
**Replace with current controls** only after a comparison is set up, and starts
no LOSAT run before that. A CLI Linear Session written without `-b` or
`--losat` commits only a disabled protein pipeline (mode `none`,
no pairs). It loads with **No comparison** and without selecting LOSATP, so
Generate draws no ribbons and starts no LOSAT run. A CLI Session written with
`--losat losatp` and an editable protein pipeline loads with the
adjacent LOSATP comparison it drew. Without saved
`ui.layoutPreferences`, the legend and plot-title positions come from the
committed `diagramOptions.output` for its mode and grouping, so a CLI
`--legend` survives loading and the first Generate.

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
Schema 2 is accepted with sessions 41–42, 44, and 45. Unknown or malformed bindings reject
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

`webFiles.bindings.linearSeqs[].inferred_definition` stores the definition
inferred from the record that row selects. It is written when a GenBank File is
uploaded or a row selects another record, and a record's definition is its own
value, then the file default, then this value. Sessions saved before the field
existed, and Sessions written by the CLI, carry none and are not re-inferred when
loaded, so their records keep the definitions they were drawn with.

Session 44 replaces the two editable Linear visibility booleans with independent
selected modes: `linear_accession_visibility` and `linear_length_visibility`,
each set to `auto`, `show`, or `hide`. Auto resolves from the effective rendered
rows and is projected to the existing request-schema-7 booleans; the request and
Python render model did not gain a new field. A version-42 `true` becomes Show,
`false` becomes Hide, and a missing boolean becomes historical Show. If a
selected-mode field is already present, it takes precedence. Current writers do
not write the retired booleans.

Plot-title text (`plot_title`), plot-title font size (`plot_title_font_size`),
and the record-label default font size (`def_font_size`) are kept per mode in
`config.modeProfiles`; the active mode's values also remain in the
flat `form` and `adv` fields. When a Session has no per-mode entry, its saved
flat value belongs to the active mode only, and the other mode starts with an
empty title and automatic font sizes. Other explicit values saved for the
inactive mode are kept. A missing `config.linearRecordLayout` means **Arrange in
rows** is on. A missing `ui.linearTypographyLinked` means the scale and ruler
label font sizes are linked while they are equal. A Circular request supplies no
Linear display values, so a Session without a saved draft keeps Accession and
Length at Auto and Replicon off.

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
settings-only writers emit session 45, while current readers also accept the
session-42 and session-44 forms. Readers whose maximum version is 42 or 44
reject newly written session-45 files.

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
envelopes. Public typed session conversion accepts full versions 31–33, 39–42, 44, and 45;
versions 27–30 are CLI replay inputs only. Canonical schema 9 retains schema 7's
display values and schema 6's input cardinality, including selectorless `all`
inputs. Resolve a typed request
before encoding when it still contains deferred paths or collection-level transforms.

## Similarity alignment request ownership

For Linear requests, `renderRequest` schemas 8 and 9 store `recordTranslations` and
`similarityAlignment` inside `renderRequest.layout`. The nested alignment plan
is schema 2; this does not change the Session version or request schema. Every
translation has one stable `recordKey` and finite `x` and `y` values. An active
plan covers those displayed keys and stores the exact reference, each target's
Select or Skip outcome and rationale. Record presentation or region state
owns the orientation used to project each anchor center. Current readers reject
partial, malformed, mismatched, or unsupported plans. Schema 1 of the nested
plan was never released and has no reader. There is no Circular form or generic
transform matrix. A Python `SimilarityAlignmentReference` is not a persisted
form: writers store the plan it resolves to.

Web **Align…** and **Review alignment options…** take the direct ortholog
evidence for the selected Similarity Group from the orthogroup result in the
committed `renderRequest`, not from the feature catalog. An edge counts only
when both endpoints bind to current members of their own source records;
missing or unbound evidence leaves the usual Review rules in place instead of
a guess. A Web LOSATP Generate commits its typed orthogroup or Collinear result
with the request, so Align, record rotation, and Save and Load use the same
evidence without running LOSATP again.

Released request schemas 1, 2, 5, 6, and 7 remain readable. Their
`alignOrthogroupFeature` protein-setting string is confined to a reader-only
legacy path. Current writers emit neither that field nor the old Session-only
`orthogroupState.selectedOrthogroupAlignmentFeature` copy. Released legacy
Sessions with saved feature-catalog and orthogroup identity metadata materialize
that state to the current schema-2 plan before saving. Malformed or unmappable
legacy values produce an actionable error without replacing the last successful
Result. Current requests reject group-only input. Loading a saved preview does
not start LOSATP; [Save and Load Sessions](web-app.md#save-and-load-sessions)
describes when Load starts the diagram engine.

A current Session round trip retains exact feature and record identities,
Select/Skip rationale, record orientation, base translations, and the
immediate pre-align Reset baseline and trustworthy restoration evidence.
Ordinary Generate renders the saved plan through the canonical typed path without guessing a new anchor. A stable
reorder resolves by `recordKey` and biological feature identity. Source
replacement, crop, selector, and record drag clear the plan with a visible reason.
Manual Reverse keeps the plan and aligns the same anchors on the next Generate.
**Reset alignment…** offers positions only (preserving current directions) or
positions plus the latest Align's actual direction changes. Combined Reset
restores absolute pre-Align directions only for those records, including a
changed reference, and replaces later manual direction edits on those targets.
Either scope consumes plan and evidence; Undo restores them before another scope
can be chosen. See [Web Reset](web-app.md#reset-positions-or-directions).
A stale reference requires Reselect/Clear and a stale target requires Select/Skip; pending or failed repair keeps the last
successful Result. Undo/Redo restores the complete artifact. Preview-only guides,
candidate markers, and recommendation badges are never saved.

Restoration evidence is stored in `editorState.alignmentResetReceipt`, separately
from the direction-independent plan. It binds source bytes, record/selector/crop
identity and exact anchors, and stores only actual absolute direction deltas
plus any reference x correction. It never stores a review mode or Custom policy.
Save and fresh Load retain this binding; style regeneration, stable reorder and
ordinary Reverse do not rewrite the pre-Align directions. Source/selector/crop
invalidation and record drag clear evidence with the plan.

Supported older Sessions without historical direction evidence retain their
position baseline, with combined Reset disabled for missing evidence. A current
receipt with an empty direction list instead reports that the latest Align
changed no directions. Malformed or stale current evidence is rejected before
replacing the current artifact; it is not silently dropped. Session version 44,
request schema 8 and plan schema 2 remain unchanged by this restoration data.

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
Only affected record pairs rerun when one valid cache key changes. A Session
stores each raw result once: rows that share a cache key, such as two ring
files with one sequence, share one entry named by the first row.

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
