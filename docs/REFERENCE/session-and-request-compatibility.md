[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [Technical documentation](README.md) | [Compatibility history](../SESSION_COMPATIBILITY.md) | [FAQ](../FAQ.md)

# Session and request compatibility

Current writers emit session version 40 and canonical `renderRequest` schema 5.

| Persisted format | Current writer | Accepted by current readers |
|---|---:|---|
| gbdraw session | 40 | 27–33 and 39–40 |
| Canonical `renderRequest` | 5 | 1, 2, and 5 |

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
defaults. Generate again before saving when the result or feature catalog is
stale. Loading and regenerating should preserve biological identities, labels,
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
envelopes. The current canonical request is a materialized wire format. Resolve
a typed request before encoding when it still contains deferred paths,
collection-level transforms, or unresolved cardinality.

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

**Save Raw LOSAT TSV** resolves generated protein handles to stable readable
aliases. Uploaded comparison TSV is never rewritten. Export raw results for a
durable evidence record; the cache exists to avoid repeated work.

For the version-by-version record of retired input names and saved-result
formats, see the [compatibility history](../SESSION_COMPATIBILITY.md). Release
notes record when support changed; this page documents current support.
