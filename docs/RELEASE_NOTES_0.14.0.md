[Documentation home](./DOCS.md) | [Installation](./INSTALL.md) | [Changelog](../CHANGELOG.md) | [Beta history](./RELEASE_NOTES_0.14.0b0.md)

# gbdraw 0.14.0 release notes

gbdraw `0.14.0` lets you rotate a circular record to start at any base, place
individual features on a chosen lane, line up Similarity Groups across Linear
records, and run LOSAT comparisons from the command line and Python. If you
are upgrading from 0.13, read [Renamed and removed options](#renamed-and-removed-options)
and [Upgrading from 0.13](#upgrading-from-013) first: several options were
renamed. [GitHub Releases](https://github.com/satoshikawato/gbdraw/releases)
gives the publication date, [PyPI](https://pypi.org/project/gbdraw/) lists the
published packages, and [Installation](./INSTALL.md) explains each install route.

## Highlights

- Start a complete circular record at any base, in Circular or Linear
  diagrams. The sequence and its coordinates do not change.
- Place a whole feature on the main lane or an available secondary lane. Its
  labels and comparison links follow it.
- Save a Session and pick up where you stopped, with its comparison results.
  **Run Info** gives you a command that rebuilds the figure from your original
  files, or the files to replay it exactly.
- Watch each stage of **Generate Diagram** as it runs, and drag a large
  preview over blank space and comparison ribbons alike.
- Draw from Python with `draw_circular()` and `draw_linear()`, which return a
  `Diagram` you can save in any format.
- Install from a wheel or source package that includes the local web app's
  palettes and browser files.

## In the web app

In Linear mode, **Generate Diagram** with a pair set to **Upload BLAST TSV**
but no file opens a **BLAST TSV missing** dialog. It offers **Choose BLAST TSV
for #i → #j…**, **Set to No comparison and Generate**, or **Cancel**.

**Generate Diagram** shows each stage as it runs: preparing the runtime and
inputs, comparing, rendering, and finishing. When a comparison reuses a cached
search, the status says so. The status shows stages, not a percentage or a
time estimate. If you cancel, or the run fails, the previous Result stays.

Circular and Linear mode now keep their own settings, edits, and Result.
**Generate Diagram** replaces only the Result of the mode you are in, and
switching modes changes no setting; a mode without a Result says **No Circular
Result yet** or **No Linear Result yet**. Input files, LOSAT run settings, and
the rich feature popup are shared by both modes. **Reset Settings** resets both
modes, and one **Undo** restores both.

**Load Session** asks **Replace the current work?** before it replaces unsaved
changes. On gbdraw.app, **Load an example** opens any Gallery example, and the **Load
example** button in the header does the same after a Result is shown; `gbdraw
gui` has no Gallery and no such buttons. The Gallery examples were rebuilt:
the Vibrio example starts each chromosome at its replication initiator, and the
BGC and majanivirus_orthogroup examples show alignments.

Each change of a checkbox, radio button, or button is one **Undo** step. Undo
and Redo are unavailable while Generate runs, and the header says why. Error
messages name the field, track row, Depth series, table line, or setting
involved. A value gbdraw cannot use, such as a GC window of 0, is reported
instead of being replaced by a default without notice. Every visible control
has an accessible name, every **?** tip is a **Help** button that opens on
hover, focus, click, or tap, and **Escape** closes only the top layer.

Large Sessions (thousands of features) open popups and apply color choices in
seconds instead of freezing the page.

Dragging the preview now follows the pointer over blank space and over
comparison ribbons, also in large diagrams. Pan, zoom, **Fit**, and **Reset**
change only the view; they do not redraw the diagram.

**Region Annotations → Download TSV** saves the annotations you are editing as
`annotations.tsv`, including changes made since the last Generate. Load the
file again with **Import TSV** in the web app, or pass it to the command line
or Python. Download works offline and keeps each row's target and style,
including an explicit no-fill. An empty editor has nothing to download. The
[annotation-table reference](./REFERENCE/input-formats-and-tsv-schemas.md#annotation-table-fields)
lists what a round trip keeps.

The hosted app at [gbdraw.app](https://gbdraw.app/) and `gbdraw gui` share one
browser interface. Only hosted gbdraw.app uses Google Analytics 4, for
aggregate page-use counts; gbdraw never sends your genome files or diagrams to
it. The local web app and the browser wheel contain no analytics. See
[Installation](./INSTALL.md#1-hosted-web-app) for what works offline and what
only the hosted Gallery provides.

## Layout and editing

**Display start** sets which base of a complete circular record comes first:
it appears at 12 o'clock in Circular mode and at the left edge of the record in
Linear mode. Enter a 1-based coordinate of the source record. Leaving it unset
is not the same as entering 1 when the record is reverse-complemented. A
feature or match that crosses the new start is drawn in pieces; gbdraw adds no
features and changes no source coordinates. You cannot crop a record and set
its display start at the same time.

**Feature placement** moves a whole feature, including every part of a
multipart feature. **Auto** removes your placement. **Main** keeps the feature
on its usual lane. A directional lane 1 is offered only in layouts that have
one. Placements work whether overlap resolution is on or off. **Feature overlap
tolerance (bp)** defaults to 0. When two fixed placements conflict, Generate
fails and says which features conflict. The
[placement reference](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#manual-feature-placement)
lists which lanes are available and how conflicts are resolved.

Rotation and placement changes take effect when you Generate. Until then they
are pending edits: **Undo**, **Redo**, **Save Session**, and **Load Session**
keep them apart from the last Result. Feature labels, leader lines, and the
comparison ribbons of a feature follow its final placement.

In the feature popup, **Edit** groups its controls by when they apply:
**Appearance** changes update the current Result, and **Layout** changes
(feature placement and record rotation) apply on Generate. To start a record at
a feature, open **Rotate record using this feature** and choose **Put this
feature at**: **Start of the record** (the default), **End of the record**, or
**Custom position**. **Apply on Generate** collects the rotations of several
records for one **Generate Diagram**; **Apply and regenerate** redraws that
record at once.

The Features list and **Search features** keep hidden features, each with a
**Visibility** checkbox, so you can find a hidden feature and show it again.
Hiding a feature also hides its label. When you set **Label visibility** to
**On** for a feature whose label cannot be drawn, for example because the
feature is hidden, a dialog says why and lets you choose. The label text
settings, **QUALIFIER PRIORITY**, and the Circular **LABEL GEOMETRY** stay
visible while any feature's **Label visibility** is **On**, even when **Show
Labels** or **Label Mode** is **None**. Per-feature visibility and label edits
follow their feature through cropping, reverse complement, record reordering,
and copies. **Export Feature Edits TSV** and **Load Feature Edits TSV**
exchange them with the command line's `--feature_override_table`.

In **Search features**, **All** no longer searches nucleotide or amino-acid
sequences or `/translation` values. Locations are shown and searched as
1-based INSDC locations, with every part of a split location. In the Legend
editor, removed rows are listed under **Deleted items** with **Restore** and
**Restore all**, and renamed and reordered rows keep their place through
Generate and **Save Session**.

Every editor edit (stroke, fill, palette, color rule, feature visibility,
Legend rename, recolor, delete, Restore, and sort) shows on the displayed Result
as the next **Generate Diagram** draws it, also after Undo, Redo, a Result
switch, **Reset Settings**, and **Load Session**. **Reset Settings** returns
Legend renames to the default captions. Renaming a Legend row onto a deleted
row's caption opens the Legend Name Conflict dialog. A Default colors value that
differs from the selected palette counts as your color, as `-d` does over `-p`
on the command line: switching the palette or pressing Default colors **Reset**
asks before it discards your colors.

In Linear mode, **Arrange in rows** puts every record selected from a source
file on that file's row, in file and record order. Row spacing follows the
features, labels, and tracks drawn. In a Depth track, a record with no
value draws no coverage; it is not drawn as zero. Each Linear record gets its own
definition from its `/organism` and `/strain`. Circular grids, docked legends,
and titles are sized to the visible diagram, so less canvas is left empty.

Arrow head length and shaft width can now be set separately. A
`repeat_region` is drawn behind the other features in new figures; set
`repeat_region=rectangle` to get the earlier look. An older Session keeps the
repeat shape it was saved with.

## Line up Similarity Groups in Linear view

To line up records on one reference feature, choose **Align…** in the feature
popup or in the Similarity Groups drawer. When each record has one clear
anchor, gbdraw applies the alignment at once and keeps each record's
direction. When the choice is ambiguous, or when you choose **Review alignment
options…**, a review palette opens:

- Each target record shows its suggested anchor and the reason. Replace it, or
  choose **Skip**.
- Choose the directions: **Keep current directions**, **All selected arrows
  right →**, **All selected arrows left ←**, or **Custom**.

When several anchors fit equally, gbdraw suggests a unique representative or
the first stable candidate. This is a convenience, not a biological judgement.
Reversing a record reverses the whole record; the +/− strands in your source
file do not change. A record whose anchor is unknown, skipped, missing, or
unusable keeps its direction, and the palette says why.

Choices in the palette run no search. **Apply** checks the plan once. If the
final directions or reference position differ from the preview, the preview
refreshes and you apply again. If Apply fails, your choices and the previous
Result stay.

The alignment is saved with every anchor and **Skip** choice. A record's
direction is its **Reverse complement** setting. Generate, reordering
records, and a manual **Reverse complement** keep the alignment. **Reset alignment…** moves records back
to where they were just before the latest Align, and can also undo the
direction changes that Align made; its preview lists later manual edits that
this reset replaces. Each Align can be reset once: to try the other kind of
reset, use **Undo** first. Saving and loading a Session keeps what Reset needs.
When an older Session lacks it, or when Align changed no direction, the dialog
says which. **Undo** and **Redo** restore the whole Result.

In Python, pass the resolved plan in a typed request. On the command line,
`--similarity_alignment_feature` takes an exact reference feature only when
each target record has one unambiguous anchor. See
[Web alignment](./REFERENCE/web-app.md#similarity-group-alignment-in-linear-view),
[command-line alignment](./REFERENCE/command-line.md#strict-similarity-group-alignment),
and [typed Python alignment](./REFERENCE/python-api.md#typed-linear-similarity-group-alignment).

## Sessions, Run Info, and reproducing a figure

gbdraw 0.14.0 saves Sessions as session version 46, with
canonical `renderRequest` schema 9. These format numbers are separate from
the gbdraw version. A Session now keeps
the settings and edits of Circular and Linear mode separately, and **Save
Session** works before you load any source file. Sessions saved by earlier
supported versions, and settings JSON files, still open; a Session that holds
only settings needs a source file before you Generate. Similarity alignments
are saved as a resolved plan with each record's X/Y offset in bases. The
[Session and request compatibility reference](./REFERENCE/session-and-request-compatibility.md)
lists every version that opens and what happens to older files.

A Session that the command line or the Python API saves opens in the web app
with the tables it records (`-t`, `-d`, `--feature_visibility_table`,
`--label_table`, `--label_whitelist`, `--qualifier_priority`, and
`--feature_override_table`), and the next **Generate Diagram** keeps them.

**Save Session** stores what the last Result needs, including its comparison
results, together with your pending edits. If you save after editing but
before Generate, the Session keeps both the newer edits and the earlier
Result; loading it does not present those edits as already generated.

**Run Info** describes the last successful Result and offers two ways to
reproduce it:

- **Source recipe**: a command that uses your original input files and the
  public CLI options. Keep your original files for it.
- **Exact replay**: the saved Session and comparison results. **Download
  reproducibility files** gives you the Session and the generated helper files
  it lists.

Neither includes control edits you have not generated. See
[Replay boundaries](./REFERENCE/session-and-request-compatibility.md#replay-boundaries)
for which Result each one reproduces.

## Command line and Python

From Python, `read_genbank()`, `read_gff()`, `draw_circular()`, and
`draw_linear()` are available from the `gbdraw` package itself (they first
appeared in the beta). `CircularOptions` and `LinearOptions` hold the settings
of each mode. `draw_circular()` takes one record or several; `CircularLayout`
arranges several in a grid. The returned `Diagram` has `to_svg()`,
`to_bytes()`, and `save(path)`. Saving to PNG, PDF, or another non-SVG format
writes only that file, with no extra SVG.

`CircularDiagramOptions` and `LinearDiagramOptions` take
`feature_override_table`, the Python form of `--feature_override_table`.
`read_gff()` without `features` now returns CDS and other features linked by
`Parent` at the record level. For typed requests, tables, saved resources, or
Session replay, use `gbdraw.api`. See the [Python API](./REFERENCE/python-api.md) and the
[typed request reference](./REFERENCE/typed-requests.md).

New command-line options set record display and feature placement:
`--record_topology`, `--display_start_coordinate`,
`--feature_placement_table`, and `--feature_overlap_tolerance_bp`.
`--feature_override_table` sets Feature visibility, Label visibility, and
label text for individual features. To set the display of several records, use
a records table. The [command-line reference](./REFERENCE/command-line.md)
has a rotation and placement example you can run, and links to the full
option list.

The command line and Python now run LOSAT themselves, without BLAST+ or the
web app:

- `gbdraw linear --losat losatn|tlosatx|losatp` runs nucleotide, translated, or
  protein comparisons between Linear records.
- `gbdraw circular --losat losatn|tlosatx` runs the comparison rings. A
  comparison genome can be FASTA, GenBank, or DDBJ (`--conservation_sequence`).
- In Python, use `LosatSearchOptions`: `LinearComparisonOptions(losat=...)` or
  `LinearDiagramOptions.losat_search` for Linear, and
  `ComparisonRingOptions(losat=...)` with `CircularDiagramOptions.losat_search`
  for Circular rings.
- New options: `--losatn_task`, `--losat_gencode` (also a records-table column,
  `losat_gencode`), a `source` column in the comparisons table that mixes
  searched and uploaded comparisons, and `--losat_output_dir`, which writes the
  raw TSVs and a reusable `comparisons.tsv`.

On the command line and in Python, Linear comparison tables (`-b`,
`--comparisons_table`, and `linear_comparisons`) are now read in the frame of
the search: the selected, cropped record on its source
strand. A forward BLAST table therefore draws the right region on a
reverse-complemented record. If you wrote a table for the reversed display,
convert that record's coordinates with `L + 1 - x`. A row outside its cropped
record stops the run. gbdraw reads the first 12 columns by position, so extra
columns (`-outfmt "6 std qlen slen"`) work, and a `#` inside an ID no longer
cuts the row short. A malformed table is reported with its line number. In the
web app, Circular comparison rings also accept GenBank and DDBJ files.

The command line, Python, and the web app compute the same raw cache keys, so
a saved Session replays its saved searches without running LOSAT again. LOSATP now treats each input file
as one genome, as the web app does: a record is not searched against itself
unless you ask for it. E-values change only for records from a file that holds
several records. A Session saved from the command line records which LOSAT it
used: its kind, version, source, program, and command-line dialect, and for a
LOSAT outside the package only the executable's name, not its full path.
**Run Info** lists the search runtime of each displayed result.

An explicit output prefix keeps its dots. A Circular batch numbers the output
of each record. gbdraw checks every output file before it draws, and an
existing file is replaced only when you allow it (`--overwrite`). In Python,
an export failure raises `ValidationError` or `ExportError` (both subclasses of
`GbdrawError`), and `save_figure_to()` returns only the files it wrote. Formats
are written one after another: when a later conversion fails, the files
already written stay.

## Installation

Bioconda remains the recommended way to install gbdraw locally. Each package
index lists the versions published there; a version in the source code does
not mean a package exists yet. [Installation](./INSTALL.md) shows how to check
and how to install from a checkout when your version is not published.

Wheel and source packages were tested in clean Linux environments on Python
3.10, 3.11, and 3.12: the CLI, the Python API, Session replay, and non-SVG
export. The local web app was also tested from an installed package. An
installed wheel includes the web app's palettes and browser files, but not the
development tests or the hosted Gallery examples. SVG output needs only the
base package; other formats need the `export` extra and the Cairo libraries for
your platform.

## Renamed and removed options

New commands, Python code, and configuration files no longer accept the
earlier names. A Session
saved by an earlier release still opens: gbdraw rewrites the old names when it
reads the file.

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
| Web app **Add legend item** (reachable only from the browser console) | A Specific color rule with a Legend caption |
| `--align_orthogroup_feature` | `--similarity_alignment_feature` |
| `--protein_blastp_output FILE` | `--losat_output_dir DIR` (writes `DIR/losatp.raw.tsv`) |
| `LinearComparisonOptions(protein_mode=..., blastp_executable=..., candidate_limit=..., orthogroup_member_max_hits=...)` | `losat=` with `losatp_mode=`, `ncbi_blast_executable=`, `max_target_seqs=`, `member_max_hits=` |
| `LinearDiagramOptions` LOSATP fields (`protein_blastp_mode`, `protein_comparison_pairs`, `losatp_bin`, ...) | `losat_search=LosatSearchOptions(...)` with `LosatRuntimeOptions` |
| Circular `--conservation_fasta` | `--conservation_sequence` (FASTA, GenBank, or DDBJ) |
| `--conservation_table` column `comparison_fasta` | `comparison_sequence` |
| `CircularDiagramOptions(conservation_fasta_files=...)` | `conservation_sequence_files` |

The modules `gbdraw.api.canvas`, `gbdraw.api.configurators`, and
`gbdraw.circular_diagram_components` are removed. SVG `id` spellings can
change between releases; to select elements, use the documented
[semantic hooks](./REFERENCE/interactive-svg-and-semantic-hooks.md). The
internal function `gbdraw.render.export.save_figure` now warns with
`DeprecationWarning`; use `save_figure_to()` or `render_to_bytes()`.

## Upgrading from 0.13

1. Update your 0.13 scripts with the table above, and check the installed CLI
   help or the [Python reference](./REFERENCE/python-api.md). For Circular
   comparison rings, `ComparisonRingOptions` and `comparison_rings` are the
   preferred names; the older `Conservation*` names and the `conservation`
   option still work.
2. A retired CLI flag exits with status 2 and names its replacement. A retired
   Python field raises `TypeError`; there is no alias. Command lines recorded
   in a 0.12 or 0.13 Session are rewritten when the Session is replayed.
   [Retired inputs](./REFERENCE/session-and-request-compatibility.md#retired-inputs)
   lists every rewritten name.
3. Check the defaults when you compare a new figure with an old one: Circular
   shows GC content and GC skew by default, Linear hides them, and new repeat
   regions are drawn behind other features. Set these explicitly when the
   earlier look matters.
4. Keep a copy of an old Session before you open it and save it again in
   0.14.0. Let gbdraw update the file; do not edit version numbers or resource
   IDs by hand. A Circular slot whose spacing an earlier release stored as a
   factor still draws, but set explicit pixel gaps before you save it, or the
   spacing cannot be kept exactly.
5. To reproduce a finished figure, use **Exact replay** from **Run Info**; to
   keep editing, use **Save Session**. Compare output only between runs of the
   same gbdraw version: SVG bytes and text measurements can differ between
   versions.

Some command-line input that 0.13 accepted now stops the run with a message:

- `--window`, `--step`, `--depth_window`, and `--depth_step` must be positive,
  and `--scale_interval` must be greater than 0.
- `currentColor`, `inherit`, `url(#id)` and other paint references, `icc-color()`,
  and an empty color in a configuration file are no longer accepted as colors
  in the CLI and Python; Default colors are checked against one documented
  domain, also in the web app. A Session Load drops a stored Legend entry color
  outside that domain and names the entry in the Load notice.
- `-n`/`--nt` takes two letters from A, C, G, T, and U (U counts as T).
- A missing, unreadable, or malformed `-b` file, or more `-b` files than
  adjacent record pairs, is an error.
- The styling tables (Default colors, Specific colors, Qualifier priority,
  Label whitelist or blacklist, Label overrides, and Feature visibility) and
  the annotation table are no longer read with CSV quoting: a `"` is part of
  the cell value, so remove quotes your files relied on. In these tables a line
  that starts with `#` is a comment, and a `#` after other text is part of the
  value. In the styling tables, a row with more columns than the table has is
  an error that names the file and line. The feature override and feature
  placement tables are unchanged.

Some older figures can change:

- Every Circular figure drawn with the web app's default **Multi-Record
  Canvas** changes: it no longer reserves an empty depth slot, and its legend
  now lists custom slot labels, skew slots, annotations, and Depth labels.
- Label, Legend, and title widths use the font's kerning, so text measures as
  browsers draw it. PNG, PDF, EPS, and PS export place circular tick labels and
  mixed italic and roman captions where browsers do. Web PDF pages are 75% of
  their former size.
- A manual tick interval that is not a whole unit prints the decimals it needs:
  a 500 bp interval no longer reads "1 kbp, 1 kbp".
- Depth TSV positions follow crops and reverse complement. Before,
  `--region chr:601-800` drew the rows of positions 1-200, and a
  reverse-complemented record drew its coverage mirrored.
- The GFF3 CDS phase sets the reading frame. A non-pseudo CDS without
  `/translation` now starts with `M` when its 5' end is complete, its reading
  frame starts at the first base, and its first codon is a start codon of its
  `transl_table`, as in INSDC `/translation`.
  This changes **Copy aa FASTA**, interactive SVG metadata, and LOSATP input.
- On the command line and in Python, the color name `seashell` is `#FFF5EE`; it
  was drawn as seagreen.
- With separate strands and overlap resolution on, an Auto feature with no
  strand that overlaps a negative-strand feature now shares the negative-strand
  lanes.
- Colliding GFF feature IDs are now told apart using the full original order of
  the source file.
- An explicit, non-default `collinearity_anchor_mode` is now used; before, the
  default was always forced. Set `rbh` if you relied on it.

## Supported platforms

gbdraw 0.14.0 supports Python 3.10, 3.11, and 3.12. The clean-install tests
ran on Linux; Windows, macOS, and later Python versions were not tested for
this release. The hosted web app and the web app in the installed package are
both supported. An interactive SVG you saved is a separate, self-contained
file.

See [Session and request compatibility](./REFERENCE/session-and-request-compatibility.md)
for which older files open and how replay works, and
[Output formats and export](./REFERENCE/output-formats-and-export.md) for what
each export format needs.

## Known limitations

- Feature placement offers only the main lane and lane 1 in each available
  direction. Higher lanes, dragging to an arbitrary position, and placing
  single exons are not available.
- Rotation needs a complete circular record of known length. Where a
  comparison fragment has gaps, its ends are interpolated; the alignment is not
  rebuilt.
- **Source recipe** is unavailable when the command line cannot express the
  Result exactly; **Run Info** gives the reason.
- Rotation and placement edits still need **Generate**: there is no automatic
  redraw and no zoom to the selection.
- A very dense Circular figure with external labels can take a long time to
  lay out, with no upper bound. You can keep using the label controls while it
  runs; **Cancel** keeps the previous Result.
- Similarity Group alignment has no Collinear-mode controls, no anchor TSV
  input, no scored inference or ranking by support count, and no automatic
  selection across several steps.
- The hosted Gallery is not included in local installs.

## Known issues

- Editing the color of a group of features in the feature popup writes one rule
  per feature (each keyed by the feature's hash) instead of one rule for the group.
- The command line and the browser can differ in the last digits of collinear
  support scores. The figure is unchanged.
- The web app does not record a LOSAT runtime version in a Session, and a FASTA
  comparison ring gets a different default label in the web app than on the
  command line.
- The web app and the command line judge color names with non-ASCII whitespace or
  digits differently.
- A Session saved from Python with listed LOSATP `pairs` opens read-only in the
  web app.
- A Session 31 to 39 whose Label whitelist has a row with a blank keyword loads
  in the web app, but replaying it on the command line fails.
<!-- WS-D findings pending: OV-330, OV-338 -->

[Documentation home](./DOCS.md) | [Beta history](./RELEASE_NOTES_0.14.0b0.md)
