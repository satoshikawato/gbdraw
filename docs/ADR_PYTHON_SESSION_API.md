# ADR: Version 31 session documents in the public Python API

- Status: accepted; version 31 bridge implemented
- Date: 2026-07-14
- Gate review: 2026-07-15、version 31 `renderRequest` round-trip により公開 gate を開放
- Related plan: `PYTHON_API_IMPROVEMENT_PLAN.md`, Phase 4

## Context

gbdraw currently has two session shapes: GUI-authored documents and CLI sidecar
documents. Both are validated and materialized by `gbdraw.session_io`, but versions
27–30 are regenerated through `session_to_cli_args()`. For sidecars,
`cliInvocation.args` and `cliInvocation.fileBindings` supply that argument list.
GUI sessions are also converted into CLI arguments.

The public Python API has typed diagram options, but it does not yet have typed
input request models that cover every supported source: records tables, GFF/FASTA
pairs, selectors, regions, per-record depth tracks, comparison inputs, embedded
LOSAT cache data, and output/session metadata. Re-exporting the existing conversion
function would therefore make argparse option names and positional argument indices
part of the library contract.

## Decision

Session loading, materialization, regeneration, and saving remain internal in this
API revision. `load_session`, `validate_session`, `session_to_cli_args`, and
related types are not re-exported from `gbdraw.api`.

Phase 4 stops at this design gate. A public session bridge may be added only after
both circular and linear workflows have CLI-independent typed input request models
and a pure conversion function can produce those models without importing CLI
modules or returning argument strings.

The follow-up [70-field audit](DIAGRAM_OPTIONS_AUDIT.md) found no dead
`DiagramOptions` field, but it did not add typed models for record inputs, output
requests, or session metadata. The prerequisites above remain unmet.

The E0 [compatibility matrix](PYTHON_SESSION_COMPATIBILITY_MATRIX.md) separated
file acceptance, CLI regeneration, and CLI-independent typed conversion.
Versions 27 through 30 are accepted envelopes, but that version range does not
guarantee a lossless typed conversion:

- CLI sidecars can be regenerated only when `cliInvocation` schema 1 and all indexed
  file/table bindings are complete. Their conversion still depends on a CLI
  argument array.
- GUI-only sessions are converted from `config` and `files` by a best-effort
  mapping. Several table inputs, editor states, custom track slots, output policy,
  and cache/result distinctions are not reproduced losslessly.
- Version 30 acquired `config.cliOptions` and multi-record-position recovery after
  its version bump. The version number alone therefore does not identify one
  document shape.

## E0 amendment: `renderRequest` boundary

The gate remains closed, with the following decisions for the next implementation
series.

1. CLI-independent request models belong to `gbdraw.api.requests`. They are not
   exported from `gbdraw.api` until the complete session gate passes.
2. Versions 27 through 30 retain the CLI paths that rebuild an argument list. They
   are not promised as lossless inputs to a future public `session_to_request()`
   function.
3. A public bridge requires a versioned typed render-request payload in the session
   document. The planned top-level field is `renderRequest`, with its own schema
   number and mode-specific Circular/Linear body.
4. Adding `renderRequest` requires session version 31. Python and Web constants,
   writers, readers, migrations, and tests must move in the same change. No further
   unversioned additions may be required for typed conversion.
5. For the future public conversion, settings come from `renderRequest`. GUI state
   is an editor projection, saved SVG `results` are artifacts, and `cliInvocation`
   supplies CLI arguments for versions 27–30; none overrides `renderRequest`.
6. A materialization context owns embedded temporary files. Request paths are valid
   only inside that context, and cleanup must cover conversion and render failures.
7. Public symbols may be added only in the release that proves version 31
   save → load → materialize → typed request → render round-trips for both modes,
   without importing CLI parsers or exposing option names/argument indices.

The next authorized phase is therefore the session-independent request model and
version 31 schema design, not a partial re-export of the current session helpers.

## Consequences

- Python callers can use all newly documented diagram, table, interactive SVG, and
  export APIs without depending on session internals.
- Existing CLI and GUI session behavior remains unchanged.
- Session files are not advertised as a stable Python library interchange format.
- Version 27～30 files remain loadable through the Circular or Linear CLI
  `--session` option; public typed conversion starts with version 31.
- A future proposal must define temporary-file lifetime, typed version/embedded-file
  errors, LOSAT cache round-tripping, and conversion tests for both session shapes
  before adding public symbols.

This is the explicit non-public outcome allowed by item 7 of the improvement plan's
definition of done; exposing a CLI-coupled partial bridge would create a less stable
contract than leaving the boundary internal.

## Version 31 implementation amendment

The prerequisites above are now satisfied for version 31 documents. The public
bridge consists of typed Circular/Linear requests, session document
load/build/save functions, a context-owned resource materializer, typed conversion,
and request/session render functions. Python and Web rendering read version 31
settings from `renderRequest`; editor state and prior SVG results remain adjunct
artifacts.

Versions 27 through 30 deliberately remain outside public typed conversion. The
Circular and Linear CLI `--session` paths can regenerate them from
`cliInvocation` or GUI state. A caller attempting `session_to_request()` on those
versions receives `SessionVersionError` rather than a guessed request.

Every path embedded in a decoded request belongs to the active
`materialize_session(...)` context. Callers must decode and render inside the
`with` block; using the materialized session after exit raises
`SessionResourceError`.

## Web restoration amendment for versions 31–33

The Web reader treats `renderRequest` and its referenced resources as the sole
authority for render semantics. Stored `config`, `ui`, `features`,
`editorState`, and `results` are projected only through explicit Web-metadata or
artifact allowlists and cannot overwrite canonical settings. Canonical resource
tables are decoded and parsed during side-effect-free preflight; hydration does
not replay user-upload watchers.

For Linear comparison height, versions 31–33 map an invalid historical value to
Auto without a user warning. New Web, CLI, and Python config inputs continue to
reject non-finite and non-positive values.
