[Documentation home](./DOCS.md) | [About](./ABOUT.md) | [Reference](./REFERENCE/README.md) | [Session and request compatibility](./REFERENCE/session-and-request-compatibility.md)

# gbdraw 0.14.0b0 release notes

This beta introduces a small top-level Python interface for new library users.
`draw_circular` and `draw_linear` return a first-party `Diagram`, mode-specific
options replace the shared 71-field option bundle, and one Circular function now
handles both single and multi-record input. Typed request/session and table
contracts remain available from `gbdraw.api`; obsolete low-level convenience
re-exports have been removed.

These notes record what changed in this release. For the currently supported
persisted versions and migration boundaries, see
[Session and request compatibility](./SESSION_COMPATIBILITY.md).

## Linear Automatic source-card row layout

- Linear **Arrange in rows** now applies a card's row to every biological record
  loaded from that source. Record order remains card order followed by source
  order, and columns are contiguous within each row.
- A selectorless Linear card uses `RecordCardinality.ALL`; an explicit selector
  remains `EXACTLY_ONE`. An `ALL` card can set `grid_row` but cannot set
  `grid_column`, because columns are assigned after record expansion.
- Canonical `renderRequest` schema 6 persists cardinality and card-level rows.
  Typed readers retain schema 5 as an `EXACTLY_ONE` compatibility input. The
  Web promotes a selectorless Linear schema-5 card to explicit `ALL`; legacy
  multi-record inputs already carry selectors. The session envelope remains
  version 40.

## New top-level Python interface

- Import the main drawing workflow from `gbdraw`, not `gbdraw.api`.
- Use `read_genbank` or `read_gff`, then call `draw_circular` or `draw_linear`.
- Pass `CircularOptions` or `LinearOptions`; wrong-mode fields are absent instead
  of being accepted and rejected later.
- Pass one `SeqRecord` or a sequence to the same drawing function.
- Use `ComparisonRingOptions` and `ComparisonRingTrackOptions` for circular
  BLAST or LOSAT similarity rings through `CircularOptions.comparison_rings`.
  `CircularOptions.conservation` and the older `Conservation*` class names
  remain compatibility aliases.
- Use `CircularLayout` for explicit grid placement; a one-record 1×1 grid is
  valid.
- Call `Diagram.to_svg()`, `Diagram.to_bytes()`, or `Diagram.save(path)` without
  handling `svgwrite.Drawing` directly.
- `Diagram.save(path)` writes exactly the requested path. It does not create an
  additional base SVG when saving another format.

See the [Python API guide](./PYTHON_API.md) for executable examples.

## Architecture/API Phase 0–2 and beta cleanup

- Fresh Python, CLI, and Web requests now resolve one versioned mode profile.
  Circular defaults to comparison thresholds `1e-5` / `50` / `70` / `0` and
  visible GC content/skew; Linear defaults to `1e-2` / `50` / `0` / `0` and
  hidden GC content/skew. Both modes include `misc_RNA`. Linear uses
  `lightgray` for the default axis and `dimgray` when a ruler is on the axis.
- Comparison thresholds are validated eagerly on every current entry point:
  e-value and bitscore must be finite and non-negative, identity must be finite
  and within `0..100`, and alignment length must be a non-negative integer.
  Identity `100` is valid.
- Typed requests now use `CircularDiagramOptions` and `LinearDiagramOptions`.
  Their mode-specific nested contracts are `CircularRequestTrackOptions`,
  `LinearRequestTrackOptions`, `CircularOutputOptions`, and `LinearOutputOptions`.
  The former typed track names remain compatibility aliases and are separate from
  the package-root beginner option classes with the same short names.
  `CircularDiagramRequest` explicitly represents `single` and `grid`;
  `CircularBatchRequest` represents `batch` and carries one output per record.
  `CircularRequestPlan`, `CircularBatchRequestPlan`, and `LinearRequestPlan`
  normalize builder selection. Root API, fresh CLI/Web generation, current
  canonical replay, and legacy internal replay route through those planners.
  A one-record `CircularLayout` produces a valid 1×1 grid.
- Emitted SVG IDs are valid, unique, and deterministic for the same input,
  configuration, and gbdraw version. First-party interactivity uses
  `data-gbdraw-*` roles, record indexes, renderers, and orientations instead of
  depending on undocumented internal ID spelling.
- GC content and GC skew now apply `stroke_width` consistently in Circular and
  Linear rendering. Width `0` produces no visible outline; a positive width uses
  the configured stroke color, and `stroke_color="none"` remains invisible.
- The canonical CLI switches are now `--gc` / `--no-gc`, `--skew` /
  `--no-skew`, and repeatable `--depth_track`. The old executable
  `--show_gc`, `--show_skew`, `--suppress_gc`, `--suppress_skew`, `--depth`,
  and `--show_depth` spellings are removed. Supported saved sessions that
  contain those tokens are migrated by the reader before replay.
- Fresh CLI options and Python request fields also reject the retired
  `depth_tick_interval`, `feature_table`, and `collinear_max_gene_gap` names,
  Circular `multi_record_size_mode=sqrt`, Linear
  `label_placement=on_feature`, and Linear
  `track_layout=spreadout|tuckin`. Use `depth_large_tick_interval`,
  `feature_visibility_table`, `collinear_max_unit_gap`, Circular `auto`,
  Linear `above_feature`, and Linear `above|below`, respectively. CLI
  spellings add the leading `--`.
- Label configuration now uses schema-derived, mode-specific leaves. Circular
  uses `labels.circular.scope` (`none|outer|both`) and
  `labels.circular.placement` (`horizontal|radial`). Linear uses
  `labels.linear.scope` (`none|all|first|orthogroup_top`),
  `labels.linear.placement` (`auto|above_feature`), and
  `labels.linear.rotation`. The shared `labels.rendering`
  (`auto|embedded_only|external_only`) setting remains the policy for embedding
  labels in feature bodies or routing them externally.
- Fresh configuration and override inputs reject the retired flat
  `show_labels` / `allow_inner_labels` aliases, `canvas.show_labels`,
  `canvas.circular.show_labels`, `canvas.linear.show_labels`, and
  `canvas.circular.allow_inner_labels`. Supported persisted-data readers
  migrate them to the canonical `labels.*` leaves; current writers do not emit
  the old forms.
- Fresh Circular slot specifications reject `spacing`, `strict`, `compress`,
  and `reserve`. Use `inner_gap_px` and `outer_gap_px`; reservation and
  compression are derived from current slot geometry and `side`.
- Those retired forms remain readable only through supported persisted-session
  and canonical request schema 1 and 2 migration. The active
  `--annotation-table` alias for `--annotation_table` and
  `--gc_content_tick_interval` alias for
  `--gc_content_large_tick_interval` are unchanged and are not removals.
- The private `__gbdraw_legacy_spacing` transport is confined to canonical
  request schema 1 and 2 readers and is never emitted by the current schema 6
  writer. Pixel spacing migrates to explicit inner and outer pixel gaps.
  Factor-based spacing can still be replayed, but it cannot be re-saved
  losslessly by the current schema 6 writer; replace it with explicit
  `inner_gap_px` and `outer_gap_px` first.
- `gbdraw.api.canvas`, `gbdraw.api.configurators`, and
  `gbdraw.circular_diagram_components` were thin compatibility modules and have
  been removed. `gbdraw.api` no longer re-exports
  `DEFAULT_SELECTED_FEATURES`, the `assemble_*` / `build_*` helpers, or
  `DiagramOptions`, `TrackOptions`, and `OutputOptions`. The obsolete
  `plot_circular_diagram` /
  `plot_linear_diagram` save wrappers and package-level `gbdraw.render`
  aliases were removed; use the root drawing facade, `render_to_bytes`, and
  `save_figure_to`. The CLI-only helpers remain internal to
  `gbdraw.render.export`. Typed request/render/session contracts,
  mode-specific diagram/track/output options, and table APIs remain supported.
- Unused argument/session helpers, the old `suppress_gc_content_and_skew`
  config helper, the CairoSVG availability proxy, the no-op
  `labels.circular_horizontal` adapter, the unreachable `-i` / `--input`
  warning path, inert GC normalization/height config keys, and zero-caller
  layout/label/collinearity helpers were removed. Supported persisted
  config/session readers still tolerate the recognized old fields; fresh
  executable APIs no longer expose them.
- `OutputOptions.output_prefix` was duplicate state and has been removed.
  `RenderOutputRequest.output_prefix` is the sole output-prefix owner for typed
  request and session execution, and the root drawing facade is output-neutral.
  Assembly remains output-neutral.
- Explicit output prefixes preserve dots. One-record Circular batch output uses
  the prefix unchanged; multiple records use `_1`, `_2`, ... suffixes. Without
  a prefix, each record ID is used and duplicate implicit names receive the
  first available `_2`, `_3`, ... suffix. Circular grid output uses the first
  record ID by default and preserves an explicit prefix unchanged. Batch
  sessions preserve the explicit grouping and every resolved output. A record
  ID used as an implicit output name must be one filename component; path-like
  IDs, IDs containing ASCII control characters, and Windows-reserved device,
  stream, or wildcard names require an explicit output prefix.
- Diagram, multi-format, batch, and sidecar targets are preflighted as complete
  sets before rendering. Existing regular files require current overwrite
  permission; directories, special files, dangling symlinks, invalid parents,
  and targets that appear after preflight fail closed. Package-root
  `Diagram.save()` and strict typed exports create final files exclusively, and
  session sidecars use a same-directory staged commit. Multi-format generation
  remains sequential rather than transactional, so completed formats survive a
  later conversion failure.
- Canonical `renderRequest` schemas 5 and 6 persist explicit `single`, `grid`,
  or `batch` grouping. Both schemas store one output
  object for a single diagram or grid and an output array for Circular batch.
  Record loading is mode-neutral; planners own topology warnings and mode,
  comparison, and cardinality policy.

Active and public runtime collinearity configuration uses
`LosslessCollinearityParameters`; canonical request schemas 1 and 2 privately
migrate legacy `standard` parameter payloads while preserving their effective
fields. Current schema 6 accepts only the lossless form.

Phase 2 completes the internal state/planner consolidation:

- Frozen `CircularRenderProfile` and `LinearRenderProfile` values resolve the
  active typed configuration once before canvas and drawing-configurator
  construction.
- `CircularRecordRenderContext` carries resolved radial layout state.
  `LinearRecordRenderContext` carries the Linear profile, resolved track layout,
  and Depth availability; its active feature-track layout, strandedness, and
  axis-ruler state are derived read-only values.
- Record placement, annotation slot planning, and common track-slot parsing now
  have one shared owner. Reusable Circular radial contracts live in
  `gbdraw.layout`, and feature render groups require prepared
  `FeatureBuildResult` layers.
- The temporary shared `DiagramOptions`, `TrackOptions`, and `OutputOptions`
  implementation bridge is removed. Builders accept only the mode-specific
  typed option contracts.

## Python/Web session version 40

- gbdraw 0.14.0b0 writes session version 40 and canonical `renderRequest`
  schema 6. Readers accept session versions 27–33 and 39–40; the public typed
  bridge accepts versions 31–33 and 39–40.
- Version 40 stores file bytes once under `resources`, with `webFiles` binding
  them to active and inactive inputs. Version 39 sessions with legacy embedded
  `files` remain readable.
- Version 40 stores one base SVG Result per logical diagram and a schema-3
  feature catalog under `editorState.featureCatalog`. Readers collapse paired
  base and interactive Results from supported older sessions.
- Linear comparison draft intent lives under `config.linearComparisonPlan`
  with mode `none`, `adjacent`, or `selected`; placement alone remains under
  `config.linearRecordLayout`. Version 40 writers omit global `blastSource`,
  nested layout comparisons, and per-record BLAST and LOSAT filename fields.
  Supported pre-40 Web sessions migrate directly to the plan.

## Compact LOSATP runtime handles

- Session version 39 introduced compact runtime handles with canonical
  `renderRequest` schema 5. Session version 40 and request schema 6 retain them.
  Public typed
  conversion is available for canonical session versions 31–33 and 39–40;
  versions 27–30 remain CLI replay inputs.
- Generated protein FASTA, raw LOSATP QUERY/SUBJECT fields, protein maps, and derived comparison references now use deterministic session-global handles with the form `h_[a-z2-7]{26}`. The handles bind a record instance to the complete CDS feature identity without repeating long feature hashes or readable aliases in every hit row. Upload filenames, modification times, session resource names, display aliases, and reopen time do not determine the handle or raw-cache identity.
- Current protein raw cache entries use schema 4, derived comparison payloads use schema 3, and the protein identity manifest uses schema 2. The manifest remains the authority for complete feature identity and display metadata and separates runtime binding from display binding. Nucleotide raw cache entries intentionally remain schema 2, and mixed protein/nucleotide caches are validated by entry type.
- **Save Raw LOSAT TSV** hydrates generated protein results at download time: it resolves every internal handle through the manifest and replaces only QUERY and SUBJECT with readable, percent-encoded protein or feature aliases. Duplicate aliases receive deterministic short ordinals. Row order, columns 3–12, numeric text, comments, and line endings are preserved; an unresolved or wrong-binding handle aborts the download instead of exposing an internal ID. User-uploaded comparison TSV is never rewritten.
- Versions 27–33 retain their existing schema-2 protein candidate path and derived schema-1 evidence. Save-before-Generate preserves pending candidates, verified results are copied to raw schema 4, and an unverifiable candidate becomes only a pair-local miss.

## Unified feature arrow geometry and underlay rendering

- Feature rendering accepts `arrow`, `rectangle`, or `underlay` through `--feature_shape`, Python `feature_shapes`, and the Web feature editor. Existing directional feature defaults remain `arrow`.
- `--arrow_head_length_ratio` configures arrow head length as `auto` or a positive finite ratio of full feature thickness. Auto keeps the mode-specific length at full shaft width and lengthens the head by the thickness removed as the shaft narrows. Numeric head ratios remain independent of shaft width. `--arrow_shaft_width_ratio` configures the centered shaft in `(0, 1]` and defaults to `1.0`. Both settings apply globally to every feature type rendered as `arrow`, are available to Python and Web requests, and survive current session round trips.
- Numeric head lengths resolve in display space in both diagram modes. A shaft ratio of `1.0` preserves the legacy full-width outline; smaller values add visible head shoulders and narrow the shaft. Short arrows fall back to triangles, multipart bodies use the configured shaft width, and undefined-strand features remain rectangular.
- The Web rendering-options capability advances to schema 3 and rejects runtimes that do not advertise the complete `arrow`, `rectangle`, and `underlay` list and unified geometry contract. The render protocol and persisted request/session schemas are unchanged.
- New configurations render `repeat_region` as an underlay: the interval covers the full feature band behind foreground glyphs and is excluded from overlap lanes and feature labels. Use `repeat_region=rectangle` to restore the previous appearance.
- Underlays are generic to any feature type and retain resolved colors, feature legends, interactive metadata, search/edit behavior, and protein-comparison eligibility. Rendering assignments do not change feature visibility.
- Automatic feature underlays are private render-time highlights, not saved region annotations. Custom track stacks require exactly one enabled feature slot when a visible underlay exists.
- Session version 40 and canonical request schema 6 record the new default.
  Supported older sessions and schema 1/2 requests with no repeat assignment
  migrate to `repeat_region=rectangle` so visual replay remains stable.

## Python/Web session version 33

- Session version 33 updates Linear custom track geometry to schema v2. Feature
  slot `height` now sets a minimum reservation and `spacing` sets clearance to
  the adjacent outer track. New Python and Web sessions preserve these values in
  canonical Linear slot lists and restore the slot order and axis position.
- Version 32 sessions remain readable. Because feature-slot `height` and
  `spacing` were no-ops in schema v1, migration removes only those two feature
  values to preserve the old effective rendering. Numeric-track geometry,
  renderer parameters, enabled state, side, order, and axis position are kept.
  The Python typed session bridge applies the same compatibility rule before
  decoding canonical requests.
- The Web app now validates and migrates the complete session before replacing
  live state. Invalid canonical data, slot geometry, or references reject the
  import without clearing the current work.

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

### Outer diagram composition now follows visible content

Circular and Linear figures now place docked legends and plot titles against
the measured visible bounds of the completed primary diagram. The shared outer
composer uses 24 px between primary content and a docked legend, 20 px between
primary content and a docked title, 20 px between a legend and title on the
same side, and 16 px around the final figure. These are internal defaults; this
release adds no layout setting. Default canvases contain less unused space, and
titles center on primary content rather than on the width after legend docking.

Changing a docked legend or title position now changes only the outer group
transforms and final viewBox. Circular radii and track bands, Linear record and
track geometry, feature labels, and comparisons retain their record-local
geometry. Corner legend positions remain overlays and keep collision-aware
placement with 8 px of obstacle clearance.

Circular multi-record layout now packs records from authoritative bounds
reported by record assembly instead of serializing and reparsing intermediate
SVG. Generated SVGs carry internal composition metadata schema 1. The Web
editor uses that metadata for legend and title reflow, then applies any saved
user drag offset. Supported older SVG results without the metadata pass through
one explicit legacy adapter when loaded. That layout change retained request
schema 5 and session schema 40.

### Web session restore now follows canonical authority

For session versions 31–33, the Web app now reconstructs render settings and
semantic tables from `renderRequest` and its canonical resources. Stored GUI
config, UI layout fields, feature edits, and prior SVG results no longer
silently overwrite canonical render semantics; allowlisted navigation metadata
and preview artifacts remain restorable. Session imports validate before commit,
suppress upload-only watcher side effects during hydration, and restore the
pre-import state if commit fails.

Linear Pairwise Match Height now accepts only Auto or a positive finite number
at Web save/run boundaries, the CLI boundary, and the Python canvas config
boundary. Historical version 31–33 Web sessions containing an invalid value are
repaired to Auto without a warning.

Specific-color TSV import is now strict and transactional. Exact duplicate rules
are deduplicated, caption/color conflicts reject the upload without partial
changes, and file-owned SVG legend entries are reconciled by desired state.
Repeated uploads no longer create numbered duplicate captions, while recolor,
rename, and removal update only file-owned legend entries.

### Sparse depth tracks keep their logical series

Repeated `--depth_track` groups and the corresponding Python depth-track matrices can now omit a file for an individual record without changing the series identity. Use `''`, `-`, `none`, or `null` as a CLI placeholder when supplying one value per record. A group that is empty for every record remains invalid.

For a missing record/series cell, gbdraw draws no depth area, quantitative axis, or ticks. It does not substitute another record's file or treat the missing cell as zero coverage. Labels, colors, heights, shared-axis ranges, legends, and custom-slot `track_index` values remain attached to the original logical series. Linear default and custom slots and Circular multi-record diagrams use the same sparse binding behavior. In Linear mode, a missing cell reserves no vertical geometry, so that record's stack compacts without changing the logical series identity.

Web canonical requests and sessions keep Depth files as a record-major matrix in
both Linear and Circular mode: one row per displayed record and one column per
logical series. `null` cells and column order survive save/load. In a Linear
custom stack, enabled Depth slots are authoritative for the legend; entries
follow slot order, honor `legend_label`, and omit unselected logical series.

The Depth CLI, Python API matrix, settings, and canonical resource shapes are unchanged. Dense depth inputs retain their logical binding and data semantics.

### Linear tracks use record-specific measured occupancy

Linear feature lanes, labels, axes, numeric tracks, Depth tracks, annotations,
and spacers now contribute measured paint and reserve bands to one per-record
vertical plan. `--resolve_overlaps` can add feature lanes without colliding with
GC, skew, or Depth tracks; only the records that need more room are expanded.
Default tracks and custom slots use the same packer.

Feature-slot `height` is now a minimum reservation: the effective band is
`max(measured, configured)`. It does not change glyph thickness, which remains
controlled by `--feature_height`. Feature-slot `spacing` is clearance to the
adjacent outer track. Missing Depth data produces no paint, reserve, axis, or
ticks; logical identity remains independent of record-local geometry.

Comparison links now begin and end at the outer edges of the records' painted
exclusion bands without an additional endpoint pad. `--comparison_height`
remains a minimum corridor; record spacing expands around taller occupancy while
keeping links outside the track stacks. Geometry metadata schema v2 reports
record-specific paint, reserve, comparison-exclusion, and canvas bands.

Row spacing now composes body, comparison, and definition constraints with an
X-aware `CollisionBand` solver. It takes the maximum eligible clearance instead
of adding independent reservations. Comparison clearance applies only at row
boundaries crossed by a comparison, and non-overlapping left-column definitions
do not enlarge the plot corridor. Single- and multi-record rows use the same
constraint policy.

In the web app, the Custom Track Slots caret controls only whether the editor is
open. **Use custom stack** controls rendering, disabling and re-enabling it
preserves the stack, and **Reset** is the only action that rebuilds the stack
from the simple controls.

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

Direct imports of `gbdraw.render.export.save_figure` remain callable for
compatibility, but the function now emits `DeprecationWarning`. It is scheduled
for removal in gbdraw 0.16. Use `save_figure_to` for file output or
`render_to_bytes` for in-memory output.

### Typed execution rejects options for the wrong mode

Previously, a non-default option for the other drawing mode in the shared
`DiagramOptions` bundle could be silently ignored. Typed requests now use
`CircularDiagramOptions` or `LinearDiagramOptions`, so the incompatible fields
are absent from the request contract. Root `CircularOptions` and `LinearOptions`
also avoid those fields at construction.

## Added integration capabilities

The following remain available from `gbdraw.api` for CLI, web, session, and custom
integration work:

- Circular multi-record typed layout: `CircularMultiRecordOptions`.
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
- `CircularDiagramOptions` and `LinearDiagramOptions` fields for the DataFrame
  and file forms of label whitelist, qualifier-priority, and label-override
  inputs.
- Region annotation models, TSV loading, coordinate/feature resolution, Circular and Linear annotation track rendering, and interactive SVG annotation metadata.

New drawing code should prefer the top-level interface described above.

## Lower-level contract

- Low-level assemblers and builders remain implementation modules, but are no
  longer re-exported from `gbdraw.api`. New integrations should use the root
  drawing facade or typed request/render contracts.
- The existing `collinearity_anchor_mode="rbh"` default remains unchanged.
- For each label-table input, pass either a DataFrame or a file path. Passing both
  forms raises `ValidationError` instead of choosing one silently.
- Catch `GbdrawError` for all expected gbdraw library failures, or catch
  `ValidationError` and `ExportError` separately when recovery differs.

## Session API boundary

The gbdraw 0.14.0b0 public session bridge accepts canonical documents from
versions 31–33 and 39–40.
`load_session_document`, `build_session_document`, `materialize_session`,
`session_to_request`, and `render_session` are exported from `gbdraw.api` and use
the typed `renderRequest` payload rather than CLI argument names or positions.

gbdraw 0.14.0b0 writes session version 40 and canonical `renderRequest` schema
5. Version 39 introduced the canonical `ui.layoutPreferences` tree for
Circular-single, Circular-multi, and Linear legend/title preferences; version
40 retains it, and supported older fields migrate on load. Schema 5 persists
explicit Circular `single`, `grid`, or `batch` grouping; batch output is an
array with one resolved entry per record, while other requests use one output
object. `renderRequest.output.prefix` is the sole output-prefix owner; schema 1
and 2 readers ignore the legacy nested prefix. Versions 27 through 30 can be
regenerated with `gbdraw circular --session` or
`gbdraw linear --session`. Public typed conversion rejects them with
`SessionVersionError` instead of reconstructing a request from `cliInvocation` or
GUI state. Paths decoded by `materialize_session` remain valid only inside the
active materialization context.

[Documentation home](./DOCS.md) | [About](./ABOUT.md) | [Reference](./REFERENCE/README.md) | [Session and request compatibility](./REFERENCE/session-and-request-compatibility.md)
