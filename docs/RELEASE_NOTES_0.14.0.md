[Documentation home](./DOCS.md) | [Installation](./INSTALL.md) | [Changelog](../CHANGELOG.md) | [Beta history](./RELEASE_NOTES_0.14.0b0.md)

# gbdraw 0.14.0 release notes

These notes describe version `0.14.0` and migration from 0.13. Check
[GitHub Releases](https://github.com/satoshikawato/gbdraw/releases) for publication
dates and [PyPI](https://pypi.org/project/gbdraw/) for available packages.
See [Installation](./INSTALL.md) for distribution and source-install routes.

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

## Linear Similarity Group alignment

The Web feature popup and Similarity Groups drawer each provide one **Align…**
action for an exact reference feature. Resolved plans apply automatically with
record directions preserved; ambiguity or **Review alignment options…** opens
the review palette. Python preselects usable anchors with visible
reasons; ambiguous recommendations use a unique representative or stable
candidate 1 as a convenience heuristic. Each target can replace its anchor
or choose **Skip**. Exclusive **Keep current directions**, **All selected arrows
right →**, **All selected arrows left ←**, and **Custom** choices include the
exact reference and selected known-direction anchors. Whole records reverse;
source +/− strands remain unchanged. Unknown, skipped, missing and unusable
anchors keep their directions with reasons. Local choices start no Worker job;
Apply validates once per attempt. Changed final facts refresh the preview and
require another Apply. Failures retain editable choices and the previous artifact.

A successful Apply stores a fully resolved schema-2 plan containing anchors and
Skip decisions. Record Reverse settings own orientation. Ordinary Generate,
stable reorder, and manual Reverse preserve the plan. **Reset alignment…**
restores immediate pre-Align positions and optionally the directions actually
changed by the latest Align. Its preview identifies later manual edits that
combined Reset replaces. Both scopes consume restoration evidence; Undo is
needed before trying the other scope. Save/fresh Load preserves valid evidence;
missing old evidence and an empty current change list have distinct reasons.
Undo/Redo restores the complete artifact. The typed Python request accepts a
resolved plan, while the CLI accepts an exact reference only when target choices are
unambiguous. See [Web alignment](./REFERENCE/web-app.md#similarity-group-alignment-in-linear-view),
[CLI behavior](./REFERENCE/command-line.md#strict-similarity-group-alignment),
and [typed Python usage](./REFERENCE/python-api.md#typed-linear-similarity-group-alignment).
Collinear alignment controls, anchor TSV, scored inference, support-count
ranking, and multi-hop automatic selection are unsupported.

## Session / replay / save compatibility

Current writers emit session version 44 and canonical `renderRequest` schema 9.
Save Session also preserves settings before the first source is loaded. Supported
older Sessions and legacy settings JSON remain readable; settings-only Sessions
need a biological source before rendering.
Schema 8 stores Linear Similarity alignment as an exact nested schema-2 plan
and finite per-record X/Y base translations. The Session version and request
schema did not change for this nested plan update. The old protein-setting
string remains a reader-only compatibility input and is not written by current
Sessions.
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
and `--feature_overlap_tolerance_bp`. `--feature_override_table` sets the
Feature visibility, Label visibility, and label text of individual features. Use a records table for multiple display
targets. The [command-line reference](./REFERENCE/command-line.md) includes a
reproducible rotation/placement example and links to the generated option inventory.

The CLI and Python API run LOSAT directly. `gbdraw linear --losat
losatn|tlosatx|losatp` and `gbdraw circular --losat losatn|tlosatx` search
nucleotide, translated, and protein comparisons without BLAST+ or the Web app.
Python uses `LosatSearchOptions` (`LinearComparisonOptions(losat=...)`,
`LinearDiagramOptions.losat_search`) for Linear and `ComparisonRingOptions(losat=...)`
with `CircularDiagramOptions.losat_search` for Circular rings. New controls are
`--losatn_task`, `--losat_gencode` (and the records-table `losat_gencode`
column), the comparisons-table `source` column that mixes searched and uploaded
edges, and `--losat_output_dir` for raw TSVs and a reusable `comparisons.tsv`.
Circular comparison genomes may be FASTA, GenBank, or DDBJ
(`--conservation_sequence`). CLI, Python, and Web searches share raw cache keys,
so a saved Session replays without LOSAT. LOSATP now uses the Web source-file
database scope: one file is one genome, and a record never searches itself
unless requested; E-values change only for records from a multi-record file.
CLI Sessions record the native runtime (kind, version, source, path, program,
CLI dialect), and Run Info lists the search runtime of each displayed result.
See the [command-line reference](./REFERENCE/command-line.md) for the options.

Explicit output prefixes retain dots; Circular batch output uses numbered
suffixes for multiple records. Output targets are checked before rendering,
and existing files require explicit overwrite permission. Python export failures
raise `ValidationError` or `ExportError` (both `GbdrawError` subclasses);
`save_figure_to()` returns only files actually written. Multi-format exports
remain sequential: earlier completed files survive a later conversion failure.

## Packaging / installation

Bioconda remains the recommended local installation route. Each package index
provides the versions published there; the source version alone does not establish
package availability. [Installation](./INSTALL.md) explains how to check the
distribution and install a checkout when the desired version is unavailable.

Wheel and sdist installation has been verified in isolated Linux environments
on Python 3.10, 3.11, and 3.12, including CLI, Python API, session replay, and
non-SVG exports. Installed wheel contents exclude development tests, private artifacts,
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
| `--protein_blastp_mode pairwise` / `orthogroup` / `collinear` | `--losat losatp --losatp_mode pairwise` / `similarity_groups` / `collinear`; `none` is omitted |
| `--losatp_bin`, `--ncbi_blastp_bin`, `--losatp_threads` | `--losat_bin`, `--ncbi_blast_bin`, `--losat_threads` |
| `--protein_blastp_max_hits`, `--protein_blastp_candidate_limit` | `--losatp_max_hits`, `--losatp_max_target_seqs` |
| `--align_orthogroup_feature` | `--similarity_alignment_feature` |
| `--protein_blastp_output FILE` | `--losat_output_dir DIR` (writes `DIR/losatp.raw.tsv`) |
| `LinearComparisonOptions(protein_mode=..., blastp_executable=..., candidate_limit=..., orthogroup_member_max_hits=...)` | `losat=` with `losatp_mode=`, `ncbi_blast_executable=`, `max_target_seqs=`, `member_max_hits=` |
| `LinearDiagramOptions` LOSATP fields (`protein_blastp_mode`, `protein_comparison_pairs`, `losatp_bin`, ...) | `losat_search=LosatSearchOptions(...)` with `LosatRuntimeOptions` |
| Circular `--conservation_fasta` | `--conservation_sequence` (FASTA, GenBank, or DDBJ) |
| `--conservation_table` column `comparison_fasta` | `comparison_sequence` |
| `CircularDiagramOptions(conservation_fasta_files=...)` | `conservation_sequence_files` |

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
5. For LOSAT, rename the flags and fields in the table above; a retired CLI flag
   exits with status 2 and names its replacement, and a retired Python field
   raises `TypeError` (no alias). Recorded 0.12/0.13 Session arguments are
   rewritten on replay; the complete list is under
   [Retired inputs](./SESSION_COMPATIBILITY.md#retired-inputs).

| Earlier input | Replacement |
| --- | --- |
| `--protein_blastp_mode orthogroup` | `--losat losatp --losatp_mode similarity_groups` |
| `--protein_blastp_mode collinear` / `pairwise` | `--losat losatp --losatp_mode collinear` / `pairwise` |
| `--losatp_bin X` | `--losat_bin X` |
| `--ncbi_blastp_bin X` | `--ncbi_blast_bin X` |
| `--losatp_threads N` | `--losat_threads N` |
| `--protein_blastp_max_hits N` | `--losatp_max_hits N` |
| `--protein_blastp_candidate_limit N` | `--losatp_max_target_seqs N` |
| `--align_orthogroup_feature ID` | `--similarity_alignment_feature ID` |
| `--protein_blastp_output FILE` | `--losat_output_dir DIR` |
| `--conservation_fasta` | `--conservation_sequence` |
| `comparison_fasta` column | `comparison_sequence` |
| `conservation_fasta_files` | `conservation_sequence_files` |

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
- Very dense Circular external-label layouts can require a long computation
  until completion or cancellation. There is no completion-time guarantee.
  Existing label selections and controls remain available; cancellation keeps
  the previous successful Result.
- The hosted Gallery is not bundled with local installs. Windows/macOS installed
  package verification remains outstanding.

[Documentation home](./DOCS.md) | [Beta history](./RELEASE_NOTES_0.14.0b0.md)
