# ADR: Keep session orchestration outside the public Python API

- Status: accepted
- Date: 2026-07-14
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

## Consequences

- Python callers can use all newly documented diagram, table, interactive SVG, and
  export APIs without depending on session internals.
- Existing CLI and GUI session behavior remains unchanged.
- Session files are not advertised as a stable Python library interchange format.
- A future proposal must define temporary-file lifetime, typed version/embedded-file
  errors, LOSAT cache round-tripping, and conversion tests for both session shapes
  before adding public symbols.

This is the explicit non-public outcome allowed by item 7 of the improvement plan's
definition of done; exposing a CLI-coupled partial bridge would create a less stable
contract than leaving the boundary internal.
