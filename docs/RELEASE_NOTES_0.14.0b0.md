[Home](./DOCS.md) | [Python API](./PYTHON_API.md) | [Export](./EXPORT.md)

# gbdraw 0.14.0b0 release notes

`gbdraw.api` grows from 70 to 87 exported symbols in this beta. The new entries
cover table readers, multi-record circular diagrams, and interactive SVG metadata.
Two behavior fixes make accepted options and returned paths match the work
performed by the library.

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

GUI and CLI session replay remains an internal feature in this release. Session
loading, validation, materialization, and replay are not exported from
`gbdraw.api`, because the current replay representation depends on CLI argument
names and positions.

The [session API ADR](./ADR_PYTHON_SESSION_API.md) records this decision and the
typed-request-model requirements that must be met before a public session bridge
is added. This boundary does not affect the existing CLI or GUI session workflow.

[Home](./DOCS.md) | [Python API](./PYTHON_API.md) | [Export](./EXPORT.md)
