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
- Choose which records of a multi-record file are drawn (**Choose records…**, **Draw this record**).
- Draw diagram text in any font: **Font** in the web app, `--font_family` and
  `--font_file` on the command line. SVG exports embed the characters they use;
  see [Fonts](./REFERENCE/output-formats-and-export.md#fonts).
- Install from a wheel or source package that includes the local web app's
  palettes and browser files.
- Draw a line under feature labels with a leader and a white box or halo
  behind their letters (`--label_line`, `--label_protection`; **Label Line**
  and **Label Protection** in the web app).
- Switch the web app between English, Bahasa Indonesia, Español, Français,
  Português, Tiếng Việt, ไทย, 日本語, 简体中文, 繁體中文, and 한국어 with the
  header **Language** menu.

## In the web app

**Generate Diagram** shows each stage as it runs: preparing the runtime and
inputs, comparing, rendering, and finishing. When a comparison reuses a cached
search, the status says so. The status shows stages, not a percentage or a
time estimate. If you cancel, or the run fails, the previous Result stays.

In the Result preview, **×** at the end of the **Search features** row hides the
search bar, and **Search** in the bottom toolbar shows it again with the same
query and matches. A click on a Legend row or a record definition in the
preview opens the Legend item popup or the record popup while **Layout edit**
is off.

A loaded Circular Session shows its saved Result without reading the embedded
source again: **Source records** says **Records not inspected** until you select
**Inspect source records** or **Generate Diagram**.

Circular and Linear mode now keep their own settings, edits, and Result.
**Generate Diagram** replaces only the Result of the mode you are in, and
switching modes changes no setting; a mode without a Result says **No Circular
Result yet** or **No Linear Result yet**. The Circular inputs and the Linear File
cards stay separate inputs that a mode switch keeps; LOSAT run settings and the
rich feature popup are shared by both modes. **Reset Settings** resets both
modes, and one **Undo** restores both, except the rich feature popup.

**Load Session** asks **Replace the current work?** before it replaces unsaved
changes. On gbdraw.app, **Load an example** opens any Gallery example, and the **Load
example** button in the header does the same after a Result is shown; `gbdraw
gui` has no Gallery and no such buttons. The Gallery examples were rebuilt: the
Vibrio chromosomes start at their replication initiator, and the aminoglycoside
BGC example draws gene arrows with Shaft Width Ratio 0.6 and aligns the clusters
on the neoU similarity group. The majanivirus example keeps its rows centered
and reverse-complements record 6 (LvMJNV) with that record's **Reverse complement**
setting, so its links to the neighboring genomes no longer cross.

In Linear mode, **Generate Diagram** with a pair set to **Upload BLAST TSV**
but no file opens a **BLAST TSV missing** dialog. It offers **Choose file…** for each
missing pair, **Skip these comparisons and Generate**, or **Cancel**; **Generate**
appears once every file is chosen.

Choose which records of a file are drawn. Each Linear record card has **Draw
this record**, and a file with several records has **Choose records…**, a list
that searches and sorts by file order, record ID, or length and turns records
on or off in bulk. **Choose records…** also sits beside **Record** in Circular,
where Multi-Record Canvas and separate diagrams draw the records that are on. A
record that is off is left out of Generate, comparisons, row layout, and the
Legend, keeps its settings (also in a saved Session), and returns in its place
when turned on; **Delete settings** removes its card settings, feature edits
and Feature placements, annotations that target it, and comparison pairs that
name it. LOSAT keeps the whole subject file as its E-value database. After an
upload, a file with more than 20 records opens the list, and turning off the
last drawn record of a file asks whether to remove the file. In 0.13.0, a
Linear File card drew every record or the one record named in Record ID, and
Multi-Record Canvas drew every record.

A row of icon buttons above the settings lists the sections of the current
mode. A button opens its section and brings it to the top, so a section far
down the panel (**Labels** in Linear mode with several files) is one click
away. The section title shows on hover and on keyboard focus.

**Collapse settings to icons** turns the settings panel into a narrow column of
icons with one section beside it, and **Close section** leaves only the column,
which gives the diagram about 265 px more width in a 1600 px window. **Expand
settings** brings the full panel back. The browser remembers this choice per site
address (with `gbdraw gui`, until the local server stops).

The Linear file list is shorter after you load a Session: each File card folds to
one line with its file name. Click the line, or press Enter or Space on it, to
open the card; **Show in settings** opens the card it points into. A folded
line still shows a file that could not be read, a missing file, and Gencode
required. The fold is not saved in a Session, and Undo does not change it. A
long file name in a file input now ends with an ellipsis instead of being cut
at both ends.

Click a drawn element to reach its settings. With Layout edit off, the Circular
scale ticks and the Linear scale show **Axis & Scale**, a GC content or GC skew
track shows **Dinucleotide content/skew**, and the plot title shows **Titles &
Record Labels**, each in a popup beside the figure. The hint **Click to open
its settings** appears when the pointer is on one of these elements. The popup
holds the section itself: an edit there is the same setting and the same
**Undo** step as in the settings panel, and it applies at the same time.
Opening or closing the popup changes nothing in the figure. Popups of the preview now
move back inside the window when you resize it.

Each change of a setting's checkbox, radio button, or button is one **Undo**
step. **Rich Feature Popup** and the search's **Regular expression** are not
Undo steps, and Undo and Redo do not change them. Undo and Redo are unavailable while Generate runs, and
the header says why. They are also disabled while a mode-switch step waits
for a label rerender, **Generate Diagram**, or a Session save or load,
instead of a click that does nothing. The notice about discarded Undo history appears only when the memory
limit drops steps, not for the routine 30-step window. Error
messages name the field, track row, Depth series, table line, or setting
involved. A value gbdraw cannot use, such as a GC window of 0, is reported
instead of being replaced by a default without notice. Every visible control
has an accessible name, every **?** tip is a **Help** button that opens on
hover, focus, click, or tap, and **Escape** closes only the top layer.

A Legend row of color None shows an empty struck-through swatch and **None**
in the Legend tab and in the Legend item popup, as the Palette pane does,
and later live color edits of such a row, also after a **Generate Diagram**
that drew it, no longer fail.

The header **Language** menu switches the interface text at once, without a
reload, between English, Bahasa Indonesia, Español, Français, Português,
Tiếng Việt, ไทย, 日本語, 简体中文, 繁體中文, and 한국어. English is shown on a
first visit whatever the browser's language, and the choice is remembered
per browser. The figure, every export, the Session, History, and the command
line text of **Run Info** are the same in every language; numbers keep the
English format (1,000.5), and warnings from Python, the Gallery tutorials, and
the documentation stay English. Translated wording may contain errors until it is
corrected; the Bahasa Indonesia, Español, Français, Português, Tiếng Việt,
ไทย, and 한국어 wording was drafted and checked without a native-speaker
reviewer.

Count sentences in the web app (for example "1 source record inspected") use
the singular and plural forms of each language, and several English texts now
use the control names the UI shows.

Documentation: the GUI images and the Gallery tutorial media show the header
**Language** menu and the final label look.

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
**Custom position**. **Apply when I Generate** collects the rotations of several
records for one **Generate Diagram**; **Apply now** redraws that
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

Editor edits (stroke, fill, palette, color rule, feature visibility,
Legend rename, recolor, delete, Restore, and sort) show on the displayed Result
as the next **Generate Diagram** draws them, also after Undo, Redo, a Result
switch, **Reset Settings**, and **Load Session**. **Reset Settings** returns
Legend renames to the default captions. Renaming a Legend row onto a deleted
row's caption opens the **Legend Name Conflict** dialog. A Default colors value that
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
options…**, the **Review alignment** panel opens beside the figure (below it in
a narrow window) with one card per record. Each card draws its record as it
will lie after **Apply**, and rings on the figure mark the candidates; a click
on a ring chooses that gene:

- A record with several candidates lists them with how each links to the
  reference, marks one **Recommended**, and offers **Skip**. A record with one
  candidate shows a compact card with a **Skip** checkbox.
- **Arrow direction** offers **Keep as is**, **All right →**, **All left ←**,
  or **Custom**. Each card shows its arrow before and after Align.

When several anchors fit equally, gbdraw recommends a unique representative or
the first stable candidate. This is a convenience, not a biological judgement.
Reversing a record reverses the whole record; the +/− strands in your source
file do not change. A record whose anchor is unknown, skipped, missing, or
unusable keeps its direction, and its card says why.

Choices in the dialog run no search. **Apply** checks the plan once. If the
final directions or reference position differ from the preview, the preview
refreshes and you apply again. If Apply fails, your choices and the previous
Result stay.

While the panel is open, the figure is fitted again whenever its room changes:
when you resize the window, drag the settings panel wider or narrower, switch
it to icons or back, open or close a section beside the icons, the review
panel moves between the side of the figure and below it, or a notice above the
Result appears or is dismissed. A zoom or position you set by hand during the review
stays until the next such change, and the view you had before the review
returns when the panel closes. Dragging the settings panel's width now follows
the pointer from the first move.

The alignment is saved with every anchor and **Skip** choice. A record's
direction is its **Reverse complement** setting. Generate, reordering
records, and a manual **Reverse complement** keep the alignment. Changing a
pair's source or file, the comparison program, or the LOSATP mode clears it,
except that choosing **Upload BLAST TSV** for a pair keeps it until a file is
chosen. Result filters such as **Minimum identity** keep it, even when they
change the Similarity Groups.
**Reset alignment…** moves records back
to where they were just before the latest Align. Choose **Positions and
directions** to also undo the direction changes that Align made; it lists the
records it restores and marks later manual direction edits that this reset
replaces. Each Align can be reset once: to try the other kind of
reset, use **Undo** first. Saving and loading a Session keeps what Reset needs.
When an older Session lacks it, or when Align changed no direction, the dialog
says which. **Undo** and **Redo** restore the whole Result.

The feature popup's **Similarity group** section and the similarity-group and
collinear-block ribbon popups now lead with one summary line, and so does the
Editor's **Similarity groups** tab while an alignment is active. Each member
has a **FASTA** menu, the group's FASTA copy and download are in one **Copy**
menu (**More** in the feature popup), and **Reset Rename** is now **Reset
name**.

In Python, pass the resolved plan in a typed request. On the command line,
`--similarity_alignment_feature` takes an exact reference feature only when
each target record has one unambiguous anchor. See
[Web alignment](./REFERENCE/web-app.md#similarity-group-alignment-in-linear-view),
[command-line alignment](./REFERENCE/command-line.md#strict-similarity-group-alignment),
and [typed Python alignment](./REFERENCE/python-api.md#typed-linear-similarity-group-alignment).

In Linear mode, **Paralog links per group** also appears in LOSATP **Similarity
groups** mode, and a Similarity groups Generate draws with the value it shows.
Before, every such Generate drew with 2.

## Sessions, Run Info, and reproducing a figure

gbdraw 0.14.0 saves Sessions as session version 48, with
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

A region annotation drawn as an overlay can span several tracks: add
`to_slot=<slot id>` next to `anchor_slot` on the annotation slot, and the
overlay runs from one slot's band to the other's, for example one tint per
region from the GC track to the feature track. It works in Circular and Linear,
from the command line, Python, a Session, and the web app, where the
annotation slot's **Extend to** list sets it. See
[Region annotations](./CLI_Reference.md#region-annotations).

A region label without its own `label_font_size` is drawn at the feature label
font size of its diagram. In Circular, labels with `label_orientation` `auto`
(the default) or `arc` are drawn along their circle; a label longer than its
circle is cut with an ellipsis and a warning. In Linear, a region label that
meets feature labels or leader lines is moved along its region, moves those
feature labels to another row, or is drawn in front of the leader lines; when
none of these clears it, it warns.

Small label defects of horizontal Circular labels, up to 100 per record, are
now repaired by moving outer and inner labels slightly: overlapping labels,
crossing leader lines, and inner labels over the definition text in the centre.
A record is changed only when the repair clears all of them. After a repair,
the neighbouring labels make room around the repaired spot, so the leader lines
there become shorter and straighter. A repaired map takes about 1 to 8 s
longer to draw, and its canvas can grow by up to about 24 px; maps that fit
without a repair are drawn as before.
Labels that would still overlap are no longer drawn: Circular outer, inner, and
embedded labels in either placement, and Linear labels. The command line exits
with an error and Python raises `LabelsDoNotFitError`; the message lists
settings that fit, such as a larger outer-label offset, a smaller font, all
labels outside, radial labels, or `--resolve_overlaps` for labels on features
drawn on top of each other. Radial labels that need more room than the
circular tracks leave stop the same way instead of being drawn with crossing
leaders or failing with a track error. Leader lines that cross, and labels over
the definition text in the centre, are drawn as before with the same settings as
a warning. With 50 or more outer labels, `--circular_label_spacing` now also
sets the gap between the packed rows when those rows fit. See [Labels that do
not fit](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#labels-that-do-not-fit).

Circular labels keep the features' angular order across the origin (the seam
of the circle), for outer and inner labels, so fewer leaders cross. Some
figures that stopped with "labels do not fit" now draw. When set track radii
(`r=`) leave region labels on other tracks, the label warning counts them
(`N region labels on other tracks`) and the figure is unchanged.

A feature label with a leader line has a line and a protection setting per
diagram mode: `--label_line {none,under}` and
`--label_protection {none,halo,box}` (configuration `labels.<mode>.line` and
`.protection`; web app **Label Line** and **Label Protection** under
**Labels**). Linear defaults to `under` and `box`: an underline joined to the
leader, under horizontal labels only (a turned label keeps its vertical
leader and gets none), and a white box at 80 % opacity behind the letters.
Circular defaults to `none` and `none`; `--label_protection halo` draws a
white outline around the letters. Leaders and underlines are 0.07 em wide, at
least 1 px. In Linear the underline takes no room in its label row, labels
keep their rows, and region labels keep off the underline and box as off the
letters. Neighbouring Linear feature labels in one row keep at least 0.25 em
apart. The web app's PDF export draws labels without their halos and
draws centred and end-aligned text where the SVG draws it.

In Linear `--label_placement above_feature`, a horizontal label now starts one
`auto` row off its own feature, with a leader from the feature, and gets the
underline and the white box as in `auto`; a label that must move off another
takes the next `auto` row. A turned label (`--label_rotation`) keeps its place
on its feature and, with `--label_protection box`, gets the white box. Figures
with horizontal `above_feature` labels change height: most get taller, and a
figure whose labels were stuck between rows can get lower; `auto` figures
change only where a label lay over a record's definition text, which feature
labels now keep off. `--label_rotation` of a whole turn (360, -360, 720) places labels as
horizontal.

Figure size changes: a Circular annotation overlay that spans split feature
lanes, or spans nothing, now reserves room for its labels, so the canvas can
grow when it lies at the outer edge (OV-1359, OV-1360), and a `to_slot` span
from split lanes to an inside slot no longer reserves room past that slot.
Linear definition text keeps its `definition_gap` from the tick labels of a
depth or GC% axis drawn at a record's left (OV-1399), and the figure's left
edge covers those labels (OV-1403): such figures are about 20 to 35 px wider,
and a GC% track beside a depth axis can change height by about a pixel.

A CLI Session records `--label_blacklist` as keywords, not as a Blacklist
File; an older CLI Session keeps its binding.

`gbdraw render --session FILE` draws the diagrams saved in a Session;
`--list_drawings` lists them and `--drawing ID` picks one. In Python,
`SessionDocument.drawings` lists them and `render_session_drawings()` draws
several at once.

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
`GbdrawError`), and `save_figure_to()` returns only the files it wrote. Each
diagram writes its formats all or nothing: when a conversion fails, no existing
file changes and no new file is left, with or without overwrite. A Circular
batch writes record by record, so the records before a failing one stay
written.

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
| Configuration `labels.stroke_width` (leader line width) | Removed: leaders are 0.07 em wide; an override fails with the unknown-path error |
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
   A figure whose labels overlapped in 0.13 now has those labels moved
   slightly or stops with an error; rerun it with one of the settings the
   message lists.
4. Keep a copy of an old Session before you open it and save it again in
   0.14.0. Let gbdraw update the file; do not edit version numbers or resource
   IDs by hand. A Circular slot whose spacing an earlier release stored as a
   factor still draws, but set explicit pixel gaps before you save it, or the
   spacing cannot be kept exactly.
5. To reproduce a finished figure, use **Exact replay** from **Run Info**; to
   keep editing, use **Save Session**. Compare output only between runs of the
   same gbdraw version: SVG bytes and text measurements can differ between
   versions.
6. For cropped or reverse-complemented records, SVG feature IDs
   (`data-gbdraw-feature-id`, element `id`, `data-gbdraw-stable-feature-id`) and
   `hash=` selectors now use the feature's hash in its source record, as for
   untransformed records. In 0.13.0 they used the hash after the crop or reverse
   complement. Scripts that match those IDs, and TSV tables or CLI options with
   `hash=` rows written for such records, need the new values. Loading or
   upgrading a Session saved before 0.14.0 renames its **This feature only**
   color rules on cropped or reverse-complemented records to the source-record
   hash. `--session` replay of a Session 31–44 renames the same rules, and the
   label text edits in its label override table. A Session 27–30 replayed with
   `--session` keeps these rules as saved and warns. Other saved `hash=` values
   load as saved, and every value that cannot be renamed is reported.
7. SVG feature IDs: of two features with the same type and location, the one
   left drawn when the other is hidden keeps its `__instance_<source index>`
   rendered ID (element `id`, `data-gbdraw-rendered-feature-id`,
   `data-label-feature-id`), so the live Web Result and Generate agree. When
   the first feature of a record is one of them, its instance is `0` instead of
   `item`. Scripts that match these IDs need the new values (OV-1485).

Some command-line input that 0.13 accepted now stops the run with a message:

- `--window`, `--step`, `--depth_window`, and `--depth_step` must be positive,
  and `--scale_interval` must be greater than 0.
- `currentColor`, `inherit`, `url(#id)` and other paint references, `icc-color()`,
  and an empty color (in a configuration override or a color option such as
  `--depth_color ""`) are no longer accepted as colors in the CLI and Python.
- Font sizes must be greater than 0 and stroke widths 0 or greater.
- In the Specific colors table (`-t`), `None`, `NA`, and `null` are values, not
  blanks, and hex colors with alpha are rejected.
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

In the web app, **Default colors file** import checks Default colors the same
way, and **Load Session** drops a stored Legend entry color outside that set and
names the entry in the Load notice.

Some older figures can change:

- Every Circular figure drawn with the web app's default **Multi-Record
  Canvas** changes: it no longer reserves an empty depth slot, and its legend
  now lists custom slot labels, skew slots, annotations, and Depth labels.
- Legend rows that share a name are numbered by one rule. Two depth tracks
  labelled `Coverage`, which 0.13.0 drew as `Coverage` and `Coverage (2)`, now
  draw `Coverage (1)` and `Coverage (2)` when their colors differ and one
  `Coverage` row when the colors are equal. Two GC content tracks always share
  the palette's GC content color, so they now draw one `GC content` row
  instead of `GC content` and `GC content (2)`. Legend row edits of a 0.13.0
  Session stay on the rows of the same tracks; edits saved for a second row
  that now shares the first row are not used.
  In the web app, **Keep both** in the **Legend Name Conflict** dialog keeps
  the name you typed and numbers every row of that name, for example
  `CDS (1)` and `CDS (2)`, instead of giving only the renamed row a number.
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
- A Circular comparison ring without `--conservation_labels` is labelled with its
  file name without the last extension (`NC_002333.2.fna` draws `NC_002333.2`),
  as in the web app. A saved Session keeps its labels.
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

- Coloring several selected features at once in the feature popup, or renaming a
  palette row in the Legend editor, adds one Specific color rule per feature
  instead of one rule for the group, so the rules table and the Session can grow
  large.
- After a palette switch, or **Apply to all** on a palette Legend row, a color
  rule in the row's earlier color can show as one Legend row until the next
  **Generate Diagram**, which draws two.
- The command line and the browser can differ in the last digits of collinear
  support scores. The figure is unchanged.
- **Run Info** and the Session record no LOSAT version for comparisons the web
  app runs; the browser build of LOSAT carries no version number.
- In the web app, the LOSATP comparison of a Session saved from Python with
  listed `pairs` cannot be edited; **Inherit saved comparison** reuses it as saved.
- A Session 31 to 39 whose Label whitelist has a row with a blank keyword loads
  in the web app, but replaying it on the command line fails.
- Similarity Groups assignment reasons and path IDs can differ between the web app
  and the command line for a reverse-complemented record, because inference depends
  on the order of LOSAT hits. The drawn diagram is the same.
- A Session that `gbdraw circular` saves for a batch or grid with a cropped or
  reverse-complemented record stores that record as drawn, so replaying it gives its
  features different IDs than the original render. The drawing is the same.
- After you load a Session saved by 0.13.0, its Result keeps the 0.13.0 feature IDs
  on cropped or reverse-complemented records until you click **Generate Diagram**. Colors
  are unaffected. An SVG downloaded before **Generate Diagram** carries the old IDs.
- A label whitelist row that the former **Whitelist only** action saved before 0.14.0,
  for a feature without a unique qualifier, names the feature by its old drawn ID.
  On a cropped or reverse-complemented record it no longer matches that feature,
  and its label is no longer whitelisted. To show the label, set **Label visibility**
  to **On** in the feature popup, or add a matching row under **Labels** →
  **LABEL FILTERING** → **Whitelist**.

[Documentation home](./DOCS.md) | [Beta history](./RELEASE_NOTES_0.14.0b0.md)
