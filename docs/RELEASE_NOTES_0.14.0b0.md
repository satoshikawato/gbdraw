[Home](./DOCS.md) | [Python API](./PYTHON_API.md) | [Export](./EXPORT.md)

# gbdraw 0.14.0b0 release notes

This beta introduces a small top-level Python interface for new library users.
`draw_circular` and `draw_linear` return a first-party `Diagram`, mode-specific
options replace the shared 71-field option bundle, and one Circular function now
handles both single and multi-record input. Lower-level request, session, table,
and rendering components remain available from `gbdraw.api` for integrations.

## New top-level Python interface

- Import the main drawing workflow from `gbdraw`, not `gbdraw.api`.
- Use `read_genbank` or `read_gff`, then call `draw_circular` or `draw_linear`.
- Pass `CircularOptions` or `LinearOptions`; wrong-mode fields are absent instead
  of being accepted and rejected later.
- Pass one `SeqRecord` or a sequence to the same drawing function.
- Use `CircularLayout` only for a multi-record grid.
- Call `Diagram.to_svg()`, `Diagram.to_bytes()`, or `Diagram.save(path)` without
  handling `svgwrite.Drawing` directly.
- `Diagram.save(path)` writes exactly the requested path. It does not create an
  additional base SVG when saving another format.

See the [Python API guide](./PYTHON_API.md) for executable examples.

## Python/Web session version 32

- Session version 32 stores materialized annotation sets, targets, styles, and annotation track bindings in the canonical request. Version 31 canonical sessions remain readable.

- Session version 31 adds a CLI-independent typed `renderRequest` and stable embedded
  resource map shared by Python and the Web app.
- `gbdraw.api` now exposes typed Circular/Linear requests, session
  load/build/save, context-owned materialization, conversion/render results, and
  session-specific error classes.
- Version 31 rendering reads settings from `renderRequest`. Versions 27–30 can be
  regenerated with the Circular or Linear CLI `--session` option, but are not
  converted by the public typed bridge.

## Behavior corrections

### Interactive SVG search remains responsive on large diagrams

Interactive SVG search now prepares reusable field indexes and updates only changed match elements. Applying or clearing a result set uses an SVG-root search state, suppresses bulk feature transitions for two animation frames, and no longer adds a dimmed class to every unmatched feature. Previous and next navigation updates only the old and new active feature parts. The web app preview uses the same difference-based rendering behavior.

New interactive SVG exports use compact metadata schema v2. They omit precomputed FASTA text, deduplicate CDS amino-acid sequence when `qualifiers.translation` already contains it, and store raw match fields instead of expanded popup rows. The embedded runtime derives those views on demand and continues to read schema v1 payloads. Existing standalone v1 files remain unchanged and self-contained.

### Interactive match popups export genomic spans

Circular Homology-ring HSPs, Linear pairwise ribbons, and Linear collinear blocks now open one shared match popup in the web preview and standalone interactive SVG. The popup copies or downloads either genomic span, or both spans as multi-FASTA. Reverse coordinate pairs are reverse-complemented; collinear actions export the complete block envelope, which may include intergenic sequence and non-anchor genes.

Uploaded BLAST rings can supply an optional companion FASTA in the web app, with `--conservation_fasta`, or through the `comparison_fasta` conservation-table column. Without it, the displayed reference span remains available and the popup explains why the comparison span is unavailable. These are ungapped coordinate spans, not reconstructed alignments.

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

## Added lower-level integration capabilities

The following remain available from `gbdraw.api` for CLI, web, session, and custom
integration work:

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
- Region annotation models, TSV loading, coordinate/feature resolution, Circular and Linear annotation track rendering, and interactive SVG annotation metadata.

New drawing code should prefer the top-level interface described above.

## Lower-level compatibility

- Existing `assemble_circular_diagram_from_record`,
  `assemble_circular_diagram_from_records`, and
  `assemble_linear_diagram_from_records` import paths remain available.
- Existing defaults, including `collinearity_anchor_mode="rbh"`, remain unchanged.
- For each label-table input, pass either a DataFrame or a file path. Passing both
  forms raises `ValidationError` instead of choosing one silently.
- Catch `GbdrawError` for all expected gbdraw library failures, or catch
  `ValidationError` and `ExportError` separately when recovery differs.

## Session API boundary

The public session bridge accepts version 31 and 32 canonical documents.
`load_session_document`, `build_session_document`, `materialize_session`,
`session_to_request`, and `render_session` are exported from `gbdraw.api` and use
the typed `renderRequest` payload rather than CLI argument names or positions.

Versions 27 through 30 can be regenerated with `gbdraw circular --session` or
`gbdraw linear --session`. Public typed conversion rejects them with
`SessionVersionError` instead of reconstructing a request from `cliInvocation` or
GUI state. The
[session API ADR](./ADR_PYTHON_SESSION_API.md) records the version 31 boundary and
the temporary-resource lifetime contract.

[Home](./DOCS.md) | [Python API](./PYTHON_API.md) | [Export](./EXPORT.md)
