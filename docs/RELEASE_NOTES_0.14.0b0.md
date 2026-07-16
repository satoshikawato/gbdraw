[Home](./DOCS.md) | [Python API](./PYTHON_API.md) | [Export](./EXPORT.md)

# gbdraw 0.14.0b0 release notes

`gbdraw.api` grows from 70 to 117 exported symbols in this beta. The new entries
cover table readers, multi-record circular diagrams, interactive SVG metadata,
typed render requests, and canonical session documents.
Two behavior fixes make accepted options and returned paths match the work
performed by the library.

## Canonical Python/Web sessions

- Session version 31 adds a CLI-independent typed `renderRequest` and stable embedded
  resource map shared by Python and the Web app.
- `gbdraw.api` now exposes typed Circular/Linear requests, canonical session
  load/build/save, context-owned materialization, conversion/render results, and
  session-specific error classes.
- Version 31 replay treats `renderRequest` as authoritative. Versions 27–30 retain
  internal CLI compatibility but are intentionally not converted by the public typed
  bridge.

## Behavior corrections

### Collinearity anchor mode is now honored

Previously, `DiagramOptions.collinearity_anchor_mode` accepted `all`,
`one_to_one`, `rbh`, and their supported aliases, but linear assembly always used
`rbh`. The normalized requested mode now reaches the collinearity builder.

The default remains `rbh`, so callers that omit this option keep the previous
output. A caller that explicitly passed a non-default mode but depended on the old
forced-`rbh` behavior should now set `collinearity_anchor_mode="rbh"`.

### Explicit library export is strict

Previously, `save_figure_to` could return PNG, PDF, EPS, or PS paths even when an
optional converter was unavailable or conversion failed. It now returns only paths
that exist when the call completes successfully. A missing optional dependency or
unsupported environment raises `ValidationError`; a failed conversion or missing
output file raises `ExportError`. Both inherit from `GbdrawError`.

Library callers should handle the failure explicitly:

```python
# This example assumes `canvas` and `output_dir` have already been created.
from gbdraw.api import save_figure_to
from gbdraw.exceptions import GbdrawError

try:
    paths = save_figure_to(
        canvas,
        ["svg", "png"],
        output_dir=str(output_dir),
        output_prefix="diagram",
    )
except GbdrawError as exc:
    raise RuntimeError(f"gbdraw export failed: {exc}") from exc
```

PNG, PDF, EPS, and PS still require the `export` extra:

```bash
pip install "gbdraw[export]==0.14.0b0"
```

The CLI retains its warning-and-skip behavior when an optional binary converter is
unavailable. The strict contract applies to the explicit Python library helper
`save_figure_to`.

### High-level builders reject options for the wrong mode

Previously, a non-default Circular-only option passed to `build_linear_diagram`,
or a Linear-only option passed to either Circular high-level builder, could be
silently ignored. The three `build_*` helpers now raise `ValidationError` and name
the incompatible fields. Move the option to the matching builder or leave it at
its default when sharing an option bundle.

`build_circular_diagram` also rejects ambiguous or lossy legacy depth inputs. Pass
one of `depth_table`, `depth_file`, a one-element `depth_tables`, or a one-element
`depth_files`; do not combine singular/plural or table/file forms. The low-level
mode-specific assembler signatures are unchanged.

## Added public Python capabilities

The following are available from `gbdraw.api`:

- Circular multi-record layout: `CircularMultiRecordOptions` and
  `build_circular_multi_diagram`.
- Interactive SVG metadata and enrichment: `InteractiveSvgContext`,
  `build_interactive_svg_context`, and `enrich_svg`.
- Interactive byte rendering through
  `render_to_bytes(..., interactive_context=context)`.
- Records tables: `RecordsTable`, `RecordsTableRow`, `TablePathDependency`, and
  `read_records_table`.
- Conservation tables: `ConservationTable`, `ConservationTableRow`, and
  `read_conservation_table`.
- Circular track tables: `CircularTrackTable` and `read_circular_track_table`.
- Label tables: `read_label_whitelist_table`,
  `read_qualifier_priority_table`, and `read_label_override_table`.
- `DiagramOptions` fields for the DataFrame and file forms of label whitelist,
  qualifier-priority, and label-override inputs.

See the [Python API guide](./PYTHON_API.md) for executable examples and the full
capability matrix.

## Compatibility and migration

- Existing `assemble_circular_diagram_from_record`,
  `assemble_circular_diagram_from_records`, and
  `assemble_linear_diagram_from_records` import paths remain available.
- Existing defaults, including `collinearity_anchor_mode="rbh"`, remain unchanged.
- For each label-table input, pass either a DataFrame or a file path. Passing both
  forms raises `ValidationError` instead of choosing one silently.
- Catch `GbdrawError` for all expected gbdraw library failures, or catch
  `ValidationError` and `ExportError` separately when recovery differs.

## Session API boundary

The public session bridge begins with canonical version 31 documents.
`load_session_document`, `build_session_document`, `materialize_session`,
`session_to_request`, and `render_session` are exported from `gbdraw.api` and use
the typed `renderRequest` payload rather than CLI argument names or positions.

Versions 27 through 30 remain available for internal CLI replay only. Public typed
conversion rejects them with `SessionVersionError` instead of reconstructing a
request from legacy `cliInvocation` or GUI state. The
[session API ADR](./ADR_PYTHON_SESSION_API.md) records the version 31 boundary and
the temporary-resource lifetime contract.

[Home](./DOCS.md) | [Python API](./PYTHON_API.md) | [Export](./EXPORT.md)
