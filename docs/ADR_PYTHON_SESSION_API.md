# ADR: Keep session orchestration outside the public Python API

- Status: accepted, amended by E0 compatibility review
- Date: 2026-07-14
- Gate review: 2026-07-15、version 27～30 compatibility matrix 後も判断を維持
- Related plan: `PYTHON_API_IMPROVEMENT_PLAN.md`, Phase 4

## Context

gbdraw currently has two session shapes: GUI-authored documents and CLI sidecar
documents. Both are validated and materialized by `gbdraw.session_io`, but replay is
normalized through `session_to_cli_args()`. For sidecars, `cliInvocation.args` and
`cliInvocation.fileBindings` are the authoritative replay representation. GUI
sessions are also converted into the same CLI argument list.

The public Python API has typed diagram options, but it does not yet have typed
input request models that cover every replayable source: records tables, GFF/FASTA
pairs, selectors, regions, per-record depth tracks, comparison inputs, embedded
LOSAT cache data, and output/session metadata. Re-exporting the existing replay
function would therefore make argparse option names and positional argument indices
part of the library contract.

## Decision

Session loading, materialization, replay, and saving remain internal in this API
revision. `load_session`, `validate_session`, `session_to_cli_args`, and related
types are not re-exported from `gbdraw.api`.

Phase 4 stops at this design gate. A public session bridge may be added only after
both circular and linear workflows have CLI-independent typed input request models
and a pure conversion function can produce those models without importing CLI
modules or returning argument strings.

The follow-up [70-field audit](DIAGRAM_OPTIONS_AUDIT.md) found no dead
`DiagramOptions` field, but it did not add typed models for record inputs, output
requests, or session metadata. The prerequisites above remain unmet.

The E0 [compatibility matrix](PYTHON_SESSION_COMPATIBILITY_MATRIX.md) separated
envelope acceptance, internal CLI replay, and CLI-independent typed conversion.
Versions 27 through 30 are accepted envelopes, but that version range does not
guarantee a lossless typed conversion:

- CLI sidecars are replayable only when `cliInvocation` schema 1 and all indexed
  file/table bindings are complete. The authoritative representation is still a
  CLI argument array.
- GUI-only sessions are converted from `config` and `files` by a best-effort
  mapping. Several table inputs, editor states, custom track slots, output policy,
  and cache/result distinctions are not reproduced losslessly.
- Version 30 acquired `config.cliOptions` and multi-record-position recovery after
  its version bump. The version number alone therefore does not identify one
  canonical document shape.

## E0 amendment: canonical request boundary

The gate remains closed, with the following decisions for the next implementation
series.

1. CLI-independent request models belong to `gbdraw.api.requests`. They are not
   exported from `gbdraw.api` until the complete session gate passes.
2. Versions 27 through 30 retain their current internal replay paths. They are not
   promised as lossless inputs to a future public `session_to_request()` function.
3. A public bridge requires a canonical typed render-request payload in the session
   document. The planned top-level field is `renderRequest`, with its own schema
   number and mode-specific Circular/Linear body.
4. Adding `renderRequest` requires session version 31. Python and Web constants,
   writers, readers, migrations, and tests must move in the same change. No further
   unversioned additions may be required for typed replay.
5. For the future public conversion, `renderRequest` is authoritative. GUI state is
   an editor projection, saved SVG `results` are artifacts, and `cliInvocation` is
   a legacy/internal replay representation; none silently overrides the canonical
   request.
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
- Legacy version 27～30 replay remains an internal compatibility promise; public
  typed replay will start at the canonical version 31 boundary.
- A future proposal must define temporary-file lifetime, typed version/embedded-file
  errors, LOSAT cache round-tripping, and conversion tests for both session shapes
  before adding public symbols.

This is the explicit non-public outcome allowed by item 7 of the improvement plan's
definition of done; exposing a CLI-coupled partial bridge would create a less stable
contract than leaving the boundary internal.
