[Documentation home](./DOCS.md) | [Installation](./INSTALL.md) | [Changelog](../CHANGELOG.md) | [Beta history](./RELEASE_NOTES_0.14.0b0.md)

# gbdraw 0.14.0 release notes

**Status: unreleased.** These notes describe the implemented changes planned for
the final 0.14.0 release. The source version is `0.14.0` (internal final candidate);
0.14.0 has not been published to PyPI. Installation availability is documented in
[Installation](./INSTALL.md).

## Highlights

- Rotate the display start of a complete circular record in Circular or Linear
  diagrams without changing its sequence or source coordinates.
- Place a whole feature on Main or an available secondary lane, with its labels
  and comparison endpoints following the placement.
- Resume work with preserved session resources, and use **Run Info** to obtain
  a Source recipe or an Exact replay of the successful Result.
- Pan large previews reliably and see actual processing stages during Generate.
- Use the package-root Python drawing API introduced during the beta, with
  `Diagram` results and mode-specific options.
- Install self-contained wheel or sdist packages with the local GUI's palette
  data and browser runtime assets included.

## Web app improvements

**Generate Diagram** reports runtime and input preparation, comparison work,
rendering, and finalization as those stages occur. Comparison status reflects
the work being performed, including eligible cache reuse. It is stage feedback,
not a predicted percentage or completion time. A cancelled or failed run leaves
the previous successful Result available.

Preview dragging now follows pointer movement over both blank space and
comparison ribbons, including large diagrams. Pan, zoom, Fit, and Reset change
the view without regenerating the biological diagram.

**Region Annotations → Download TSV** saves the current editor draft as
`annotations.tsv`, including changes made since the last Generate. The file
works with Web **Import TSV** and the existing Python/CLI annotation-table
readers. Download works offline and preserves effective row targets and styles,
including explicit no-fill. Empty drafts cannot be downloaded. See the
[annotation-table reference](./REFERENCE/input-formats-and-tsv-schemas.md#annotation-table-fields)
for the TSV round-trip limits.

The hosted app at [gbdraw.app](https://gbdraw.app/) and local `gbdraw gui` share
the browser interface. Only hosted gbdraw.app uses Google Analytics 4 for
aggregate page-usage metrics; gbdraw does not send uploaded genome files or
generated diagrams to Google Analytics. Local Web assets and the browser wheel
have no hosted analytics injection. See [Installation](./INSTALL.md#1-hosted-web-app)
for the local/offline and hosted Gallery boundaries.

## Circular / Linear layout and editing

**Display start** accepts a 1-based source coordinate on a complete circular
record. That base appears at 12 o'clock in Circular mode or the left edge of
the wrapped record in Linear mode. An unset start differs from explicit 1 after
reverse complementation. Rotation can split a displayed feature or match into
multiple fragments; it does not create additional biological features or alter
source coordinates. Cropping and an explicit display start cannot be combined.

**Feature placement** applies to an entire feature, including multipart
features. Auto removes the manual override. Main fixes it on its nominal lane;
directional lane 1 is available only in compatible layouts. Fixed placements
work with the overlap resolver on or off. The base-pair overlap tolerance
defaults to 0, and conflicting fixed placements fail with an explanation.
See the [placement reference](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#manual-feature-placement)
for lane availability, conflict rules, and dormant placements.

Rotation and placement changes remain editable drafts until Generate succeeds.
Undo/Redo and Save/Load preserve those drafts separately from the last successful
Result. Feature labels, leaders, and feature-associated comparison ribbons
follow the final placement.

Linear **Arrange in rows** applies a source card's row to all records selected
from that source, preserving card and source order. Measured feature, label,
and track occupancy determines row spacing. Sparse Depth inputs keep their
logical series: a missing cell draws no coverage and is not treated as zero.
Circular grids, docked legends, and titles use visible diagram bounds to reduce
unused canvas space.

Arrow head length and shaft width can be controlled independently. New
configurations draw `repeat_region` as an underlay behind foreground features;
use `repeat_region=rectangle` to request the earlier appearance. Supported older
sessions retain their previous effective repeat shape.

## Session / replay / save compatibility

Current writers emit session version 42 and canonical `renderRequest` schema 7.
Save Session also preserves settings before the first source is loaded. Supported
older Sessions and legacy settings JSON remain readable; settings-only Sessions
need a biological source before rendering.
These persisted-format numbers are separate from the package version. The
[session and request compatibility reference](./REFERENCE/session-and-request-compatibility.md)
owns the accepted-reader table and migration details.

**Save Session** preserves the resources needed by the committed Result along
with editable state, including saved comparison data. Saving after edits but
before Generate retains both the newer draft and earlier Result; loading does
not silently treat that draft as an already generated result.

**Run Info** separates a Source recipe using original inputs and public CLI
settings from an Exact replay using the committed canonical session and analysis
artifacts. Its downloadable helpers remain tied to the successful generation
through supported history operations. Keep the original files for a Source
recipe; **Download reproducibility files** supplies the listed generated helpers
and Exact replay session. Neither command represents ungenerated control edits.
See [Replay boundaries](./REFERENCE/session-and-request-compatibility.md#replay-boundaries)
for the target Result and saved-editor limits.

## Python API / CLI changes

The package-root interface, already available in the beta, exports
`read_genbank()`, `read_gff()`, `draw_circular()`, and `draw_linear()`.
`CircularOptions` and `LinearOptions` keep mode-specific settings separate.
One Circular function accepts a single record or a collection; `CircularLayout`
selects a grid. A `Diagram` provides `to_svg()`, `to_bytes()`, and `save(path)`.
Saving a non-SVG format writes the requested file without an extra base SVG.

Integrations that need typed requests, tables, resource materialization, or
session replay continue to use `gbdraw.api`. See the [Python API](./REFERENCE/python-api.md)
and [typed request reference](./REFERENCE/typed-requests.md) for the supported
entry points.

The CLI adds record display and placement controls through
`--record_topology`, `--display_start_coordinate`, `--feature_placement_table`,
and `--feature_overlap_tolerance_bp`. Use a records table for multiple display
targets. The [command-line reference](./REFERENCE/command-line.md) includes a
reproducible rotation/placement example and links to the generated option inventory.

Explicit output prefixes retain dots; Circular batch output uses numbered
suffixes for multiple records. Output targets are checked before rendering,
and existing files require explicit overwrite permission. Python export failures
raise `ValidationError` or `ExportError` (both `GbdrawError` subclasses);
`save_figure_to()` returns only files actually written. Multi-format exports
remain sequential: earlier completed files survive a later conversion failure.

## Packaging / installation

Bioconda remains the recommended local installation route. The PyPI Trusted
Publishing workflow is prepared, but publication has not occurred. The future
PyPI command and current source-install route are distinguished in
[Installation](./INSTALL.md); no final package availability is implied here.

Wheel and sdist installation has been verified in isolated Linux environments
on Python 3.10, 3.11, and 3.12, including CLI, Python API, session replay, and
non-SVG exports. Package contents exclude development tests, private artifacts,
and hosted Gallery examples while retaining the GUI palette data and required
local browser assets. The local GUI was also verified from an installed package.
SVG needs only the base package; other formats need the `export` extra and
platform-appropriate Cairo libraries.

## Breaking or renamed interfaces

Fresh CLI commands and Python/configuration inputs reject retired spellings.
Supported saved documents use the dedicated compatibility readers instead.

| Earlier interface | Current replacement |
| --- | --- |
| `--show_gc`, `--suppress_gc` | `--gc`, `--no-gc` |
| `--show_skew`, `--suppress_skew` | `--skew`, `--no-skew` |
| `--depth`, `--show_depth` | Repeat `--depth_track` for each series |
| `--depth_tick_interval` | `--depth_large_tick_interval` |
| `--feature_table` / Python `feature_table` | `--feature_visibility_table` / `feature_visibility_table` |
| `--collinear_max_gene_gap` | `--collinear_max_unit_gap` |
| Circular `--multi_record_size_mode sqrt` | `--multi_record_size_mode auto` |
| Linear `--label_placement on_feature` | `--label_placement above_feature` |
| Linear `--track_layout spreadout` / `tuckin` | `--track_layout above` / `below` |
| Flat `show_labels` / `allow_inner_labels` configuration | Mode-specific `labels.circular.*` / `labels.linear.*` leaves |
| Circular slot `spacing`, `strict`, `compress`, `reserve` | Explicit `inner_gap_px` / `outer_gap_px`; geometry determines reservation and compression |
| `gbdraw.api` shared `DiagramOptions`, `TrackOptions`, `OutputOptions` | Root mode-specific options or typed mode-specific request options |
| Low-level canvas/configurator/assembler re-exports and `plot_*_diagram` save wrappers | Root `draw_circular()` / `draw_linear()`, or typed request/render helpers |
| `OutputOptions.output_prefix` | `RenderOutputRequest.output_prefix` in typed integrations |

The thin `gbdraw.api.canvas`, `gbdraw.api.configurators`, and
`gbdraw.circular_diagram_components` modules are removed. Undocumented SVG ID
spellings are not an integration contract; use the documented
[semantic hooks](./REFERENCE/interactive-svg-and-semantic-hooks.md).
The internal `gbdraw.render.export.save_figure` compatibility function emits
`DeprecationWarning`; use `save_figure_to()` or `render_to_bytes()`.

## Migration notes

1. For 0.13 scripts, update retired names using the table above and check the
   installed CLI help or Python reference. `ComparisonRingOptions` and
   `comparison_rings` are the preferred Circular names; the older
   `Conservation*` names and `conservation` option remain compatibility aliases.
2. Review rendering defaults when comparing a new figure: Circular shows GC
   content/skew by default, Linear hides them, and new repeat regions use
   underlays. Explicitly select settings when a prior appearance matters.
3. Keep the original session before opening and saving it in the new version.
   Let the reader migrate supported formats; never edit version numbers or
   resource identities by hand. Legacy factor-based Circular slot spacing can
   replay but needs explicit pixel gaps before lossless saving to the current format.
4. Use Exact replay for a successful generation with its saved analysis artifacts,
   and Save Session to resume editable work. Keep the same gbdraw version when
   comparing output; SVG bytes and text metrics can differ across versions.

Two layout/identity corrections can affect older results: overlapping
undefined- and negative-strand Auto features share the negative pool when
separate strands and overlap resolution are enabled, and colliding GFF feature
IDs are disambiguated using the complete original source order. Also, an
explicit non-default `collinearity_anchor_mode` is now honored; set `rbh` if
you relied on the previously forced default.

## Compatibility

The supported Python versions for this release are 3.10, 3.11, and 3.12. Linux
isolated-install evidence does not establish Windows/macOS or later-Python
validation. Hosted Web and the packaged local GUI are supported entry points;
saved interactive SVG files remain a separate, self-contained output.

See [Session and request compatibility](./REFERENCE/session-and-request-compatibility.md)
for supported legacy readers, save/load behavior, and replay limits, and
[Output formats and export](./REFERENCE/output-formats-and-export.md) for export
requirements.

## Known limitations / deferred work

- Manual feature placement is limited to Main and supported lane 1 directions;
  higher lanes, arbitrary pixel dragging, and per-exon placement are not offered.
- Rotation requires a complete circular record with known length. Gapped
  comparison fragments use endpoint interpolation, not reconstructed alignments.
- Source recipe is unavailable when the committed semantics cannot be expressed
  losslessly as a CLI recipe; Run Info explains the reason.
- Generate remains explicit for rotation and placement drafts. A new automatic
  redraw scheduler and new zoom-to-selection navigation are not release features.
- The hosted Gallery is not bundled with local installs. Windows/macOS installed
  package verification remains outstanding.

[Documentation home](./DOCS.md) | [Beta history](./RELEASE_NOTES_0.14.0b0.md)
