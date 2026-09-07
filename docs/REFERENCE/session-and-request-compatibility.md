[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [Technical documentation](README.md) | [Compatibility history](../SESSION_COMPATIBILITY.md) | [FAQ](../FAQ.md)

# Session and request compatibility

Current writers emit session version 41 and canonical `renderRequest` schema 7.

| Persisted format | Current writer | Accepted by current readers |
|---|---:|---|
| gbdraw session | 41 | 27–33 and 39–41 |
| Canonical `renderRequest` | 7 | 1, 2, 5, 6, and 7 |

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

## Replay boundaries

On the command line, replay a session with the same `circular` or `linear`
subcommand that created it. Output prefix, format, session-output, and overwrite
options may replace their saved counterparts. Other diagram options are
rejected because they would combine persisted and new settings ambiguously.

In Python, `render_session()` is the persisted-session entry point.
`load_session_document()` validates a document, and `materialize_session()`
exposes embedded resources only while its context is active.
`session_to_request()` converts a current materialized session to a typed
request. Rendering that request alone does not replay saved comparison
artifacts; use `render_session()` when those artifacts belong in the result.

`render_request()` accepts current typed requests, not historical session
envelopes. Canonical schema 6 records each input's cardinality, including
selectorless `all` inputs. Resolve a typed request before encoding when it still
contains deferred paths or collection-level transforms.

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
also records the effective `cds` or `locus` unit kinds produced by `auto`.

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
