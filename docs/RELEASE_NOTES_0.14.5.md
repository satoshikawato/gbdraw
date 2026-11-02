[Documentation home](./DOCS.md) | [Installation](./INSTALL.md) | [Changelog](../CHANGELOG.md) | [0.14.4 release notes](./RELEASE_NOTES_0.14.4.md) | [0.14.3 release notes](./RELEASE_NOTES_0.14.3.md)

# gbdraw 0.14.5 release notes

gbdraw `0.14.5` is a patch release of [0.14.4](./RELEASE_NOTES_0.14.4.md). The
Web app now draws a LOSATP link or Collinear block that uses a CDS across a
circular record's origin over the CDS's own arc, as the command line and the
Python API have done since 0.14.4. The release also adds `--wrap_length`, which
draws one long Linear record over several rows (command line and Python), four
built-in color tables, a fixed identity color scale, italic gene names, a transparent
Web PNG by default, and two Linear alignment options (`--auto_orient` and `--anchor_feature`, with
Web buttons),
and reads more real-world GenBank, EMBL, and GFF3 files. The Web app's LOSAT
engine is now a build of a LOSAT development commit, and threaded LOSAT searches
no longer stop on a memory fault.
The Web UI is available in 17 languages (English and 16 translations).
Sessions are now saved as session version 50 (0.14.4 writes 49; the canonical
`renderRequest` schema stays 9): **a Session saved by 0.14.5 opens only in
0.14.5 or later.** Sessions saved by 0.14.4 and earlier open as before.
[GitHub Releases](https://github.com/satoshikawato/gbdraw/releases) gives the
publication date, [PyPI](https://pypi.org/project/gbdraw/) lists the published
packages, and [Installation](./INSTALL.md) explains each install route.

## Fixed

- **Web app: LOSATP links and Collinear blocks across the origin.** In 0.14.4
  the Web app still read a CDS whose location crosses a circular record's origin
  as an interval over the whole record, and drew every LOSATP Pairwise link,
  Similarity group link, and Collinear block that used it over the whole record.
  It now draws them over the CDS's own arc, as one piece on each side of the
  origin, the same rows as the command line and the Python API. **Generate
  Diagram** again on a Session saved by 0.14.4 to redraw it: the derived LOSATP
  entries are rebuilt once from the raw searches. Records without such a CDS
  draw as before.
- **Upgrading a Session with a retired Collinear anchor mode.** Upgrading an
  older Session in Python, or re-saving it with `--session_output`, writes
  `rbh` for a retired mode such as `one_to_one`, in its requests and drafts. A
  later Web **Load Session** of such a file shows no replacement row, and a
  stored Collinear result computed with `one_to_one` is labelled `rbh` without a
  notice. A hand-edited spelling of a retired mode, such as `ONE-TO-ONE`, now
  reads and upgrades instead of being refused.
- **Real-world GenBank, EMBL, and GFF3 files.** These inputs crashed, were
  skipped, or were misread; they now read or stop with a named error:
  - A LOCUS line with spaces in the name (`LOCUS       Contig 1 edited
    by me  3000 bp`) reads; each run of spaces in the name becomes `_`.
  - An EMBL flat file reads as EMBL. The Web file pickers accept `.embl` and
    `.emb`.
  - A qualifier line indented with tabs reads as the standard columns; before,
    it made a bogus feature and the CDS lost its qualifiers.
  - A `CONTIG` record without sequence stops in both modes with an error that
    says how to download the sequence; before, Circular printed a traceback and
    Linear wrote an empty SVG.
  - A GFF3 feature type such as `cds` is drawn as `CDS`; before, the type filter
    dropped it and nothing was drawn.
  - GFF3 attributes written as `ID=g1; Name=x ; sequence=y` read; a remaining
    parser assertion is a named error that says which column to change.
  - A `--gbk` file without a record next to a good file is no longer skipped
    silently: the error names the file and says whether it is empty, FASTA, or
    neither flat-file format.
- **Web app: threaded LOSAT searches stop on "memory access out of bounds"
  (OV-2369).** A threaded LOSAT search could stop with that message and run the
  whole pair again (the console said "running the pair again"), or fail after
  the third stop. In Chromium a copy into memory that another LOSAT thread had
  just grown could trap. The Web app now checks each such copy and first brings
  the thread's view of the memory up to date. Results are unchanged.
- **LOSATP links among tied hits are the same on the command line, in Python,
  and in the Web app (OV-2370, OV-2380).** LOSATP Pairwise and Similarity groups
  broke ties between equal-score protein hits on an internal handle that differs
  between the surfaces and between Web pages, so the same records and hits could
  draw a few other links among identical paralogs. Every such tie now goes to
  the first CDS by record and location. Inputs without ties draw as before.
  **Generate Diagram** again on a saved Session: the derived LOSATP entries are
  rebuilt once from the raw searches; no search runs again.
- **`--losatp_mode pairwise` on a Linear row with two or more records
  (OV-2371).** Command line and Python only. It exited with "Ordered
  BLAST/protein comparison inputs are ambiguous". It now searches every record
  pair between adjacent rows, as Similarity groups and Collinear blocks do.
  `-b/--blast` on such a layout is still refused; the message now says to list
  the tables in `--comparisons_table`.
- **Web app: a hover popup after a pan.** After a pan or wheel zoom of the Result
  preview, the hover was replayed at the point where the pointer was released,
  and a popup could open and stay for a pointer that had already left the
  drawing (OV-2281). A pointer that stays over a match or feature still shows
  its hover.
- **UTF-8 byte order mark.** A GenBank, EMBL, GFF3, or FASTA file that starts
  with a byte order mark no longer loses its first record, or stops a one-record
  file with an error. A Session saved by 0.14.4 or earlier with such a file of
  several records counted its records without the first one, so each saved
  record index now points one record earlier; see
  [Session and request compatibility](./REFERENCE/session-and-request-compatibility.md).
- **Web app: recoloring comparison matches live.** A palette change on a Linear
  comparison, or on Collinear blocks colored by identity, now gives the same
  match colors as **Generate Diagram**, except where rounding differs by one in
  a channel value.
- **Web app: a Depth row and its last file.** A Linear Depth row now leaves with
  its last Depth file after **Generate Diagram** too, and in a Session saved
  earlier, so the next Generate draws instead of refusing the stack. A moved row
  leaves too, and a file added again puts it at the default place; a height typed
  on the row stays on its series. A Depth row you turned off now stays in the
  stack, off, when its last file goes and when a file is added again, in Circular
  and Linear; before, it left and came back on.
- **Web app: text left in English or hard to read in other languages.** The Undo
  and Redo tooltips of preserved settings, imported comparisons, and the
  **Labels do not fit** dialog, Auto values such as `above` and `ignored`, and
  that dialog's row titles and setting names stayed in English in every language;
  they now follow the language setting. Tooltips and dialogs that insert an action,
  a control, or a track row now read as sentences in more languages (Japanese
  「ファイルを削除」を元に戻す). In English, a tooltip no longer names the path or
  the row (`Undo Reset preserved setting`), the lane warning quotes the control,
  several rows read `Track rows 3 and 4`, and `Auto (subject (auto))` reads
  `Auto (subject)`. The sentences in the **Labels do not fit** dialog's notes and
  reasons, and its line with each record's label count and font size, are still
  in English.
- **Web app: corrected translations in the ten earlier languages.** Japanese,
  Simplified and Traditional Chinese, Korean, Spanish, French, Portuguese, Thai,
  Indonesian, and Vietnamese had wrong meanings, broken grammar, and inconsistent
  terms corrected throughout the UI, the same process the six new languages went
  through. The Align button of the similarity group drawer, English in every
  language until now, is translated too. All translations are still LLM drafts.
- **Web app: Reset Settings and a Result not shown.** After **Reset Settings**, the
  other mode's Result and the other Results of a batch now show the reset feature
  and label edits when you display them: features and labels hidden when they were
  generated return, an edited label text returns to its own text, and after a
  Session load the Legend rows of removed color rules leave. Before, these stayed
  until the next Generate. This redraws that Result once, with **Auto Reflow** on
  or off, when it was drawn with such an edit, even where the redraw changes
  nothing, as for a label text on a feature whose label it does not draw. A failed
  redraw reports its failure and does not run again at the next display. These stay
  until the next redraw of that Result or **Generate Diagram**:
  - a Result whose redraw failed (a later color edit on it can also correct it);
  - the label text a Result was drawn with, when **Reset Settings**, Save, and Load
    pass before it is shown again (it still draws the features it lacks), and its
    **Label visibility** too when the Reset was made in the other drawing (in the
    reset drawing, the redraw applies it);
  - a label text on a Result not shown when the edit left, returned by an Undo of a
    later **Generate Diagram**.

  Settings marked **Applies on Generate** still wait for **Generate Diagram**.
- **Web app: text typed in the Features list.** Text typed in a label input of the
  Features list and not yet committed stays when a redraw finishes.

## Changed behavior

- **Built-in color tables.** `-t/--table` and `FeatureOptions.color_table`
  accept `builtin:cog`, `builtin:phrog`, `builtin:plastid`, and
  `builtin:mitochondrion`; an unknown name fails with the list of names, and a
  file path works as before. In the Web app, **Colors > Preset scheme** lists
  "COG categories (Bakta)", "PHROG categories (Pharokka)", "Plastid gene
  classes", and "Mitochondrial gene classes" instead of "Pharokka
  (function/phrog/vfdb)" and "Bakta (COG via note)"; **Apply** replaces the
  specific rules with the table's rows. The COG preset keeps the old Bakta
  preset's CDS default `#cccccc`, Legend box, and font size 12. The old Bakta
  preset matched `/note` only and colored nothing on Bakta's default GenBank
  output, which writes `/db_xref="COG:L"`; the `cog` table reads `db_xref`, GFF3
  `Dbxref`, and `note`. Nothing is colored by a built-in table unless you name
  it, and a Session stores the rows that a preset or `builtin:` table added, not
  the name.
- **Web app LOSAT engine.** The engine is the LOSAT development commit
  `dde737f18` (version 0.2.0 before its release; it was a build of LOSAT 0.1.0),
  not a LOSAT release. Its searches give the rows of the `losat` command of that
  commit. LOSATP gives the same rows as before; on the Vibrio check its serial
  searches took 26-38 % of the old engine's time.
  LOSATN with task `blastn` gives different E-values and hits than before (on the
  lambda x DE3 check, 21 rows instead of 19); task `megablast` and TLOSATX give
  the same rows.
  - Task `dc-megablast` now runs dc-megablast (the old engine ran `blastn` for
    it) on threaded execution: **Threaded**, and **Auto** on a cross-origin
    isolated page. On serial execution (**Serial**, or **Auto** on a page
    without isolation) it stops before the search with a message that
    `dc-megablast` runs only with threaded LOSAT.
  - **Auto** now runs LOSATN and TLOSATX threaded too, and LOSATP of any size,
    wherever the browser supports threads; the label reads `Auto (threaded)` or
    `Auto (serial)`. Before, **Auto** threaded only LOSATP above 500,000 FASTA
    characters.
  - Run Info names the engine as `LOSAT 0.2.0+dde737f188dc (wasm)`.
  - Saved searches: a search saved by the command line or by this build is
    reused. A search saved by an earlier Web app is reused for LOSATP and LOSATN
    task `megablast`, which both engines give alike, and runs again once at the
    next **Generate Diagram** for LOSATN `blastn`, `dc-megablast`, and TLOSATX.
    No Gallery example searches again.
  - `gbdraw/web/wasm/losat/losat-build.json` records the build.
- **Web app PNG background.** A PNG downloaded from the Web app now keeps the
  drawing's transparent background, as SVG and PNG from the command line and
  Python already do. Check **White background**, next to the PNG DPI choice, to
  paint it white as earlier Web PNGs were. The checkbox is a page setting: a
  Session does not store it.
- **Web app: removals that other track rows depend on.** Removing a Depth series,
  the last Depth file of a series, or a custom stack row that another row starts
  at or extends to (Circular and Linear), removing or clearing a Linear File that
  holds such a last file, removing a Circular comparison source whose ring row
  leaves with it, and **Hide GC Content** or **Hide GC Skew** with **Use custom
  stack** off when it turns off such a row now ask first, with the dialog a row's
  switch uses. It offers to make the change and turn those rows off as one Undo
  step, or Cancel, which changes nothing. With no row reaching what changes,
  nothing asks and the change applies at once, as before.
- **Tobacco chloroplast tutorial and Gallery entry** use `-t builtin:plastid`
  (Web: **Plastid gene classes** in **Preset scheme**) instead of a downloaded
  color table. The figure keeps its layout and labels.

## Added

- **Web UI: German, Italian, Polish, Russian, Turkish and Ukrainian.** The
  Language menu offers these six languages, each with a complete translation of
  the UI. Polish and Turkish letters use the bundled Inter Latin Extended
  subset, and Russian and Ukrainian letters the bundled Inter Cyrillic subset.
  The figure, the exports, Session files, and the command line stay English.

- **`--wrap_length` for Linear diagrams.** `--wrap_length AUTO|LENGTH` (for
  example `50000`, `50kb`, `0.5Mb`, or `auto`) and the config key
  `canvas.linear.wrap_length` draw one long record over several rows at one
  scale, on the command line and in Python. The record stays one record in the
  SVG, the interactive catalog, and the Session. Each row has its own height;
  the ruler sits on each row's axis with absolute coordinates; a feature cut by
  a row has an arrowhead only at its true end and one label. `auto` picks the
  smallest 1-2-5 x 10^n length that keeps the longest record within 5 rows (at
  least 10 kb, at most 50 kb). With no record longer than the row, the figure is
  what the same command draws without the option. Not available while a record
  is cut: `--track_layout middle`, a `side=overlay` feature slot, a
  `--scale_style` other than `ruler`, `--normalize_length`, comparisons, and
  several records in one grid row. See
  [Wrap one record over several rows](./REFERENCE/command-line.md#wrap-one-record-over-several-rows).
  The Web app has no control for it yet: it draws and keeps a length held by a
  loaded Session, and lists it under the settings that are not editable, with
  **Reset**.
- `--scale_interval` takes a unit, for example `5kb`.
- A PNG over 32,767 px per side fails with an error that names the limit. A
  wrapped PNG over that limit, or a wrapped PDF over 14,400 pt, stops before
  writing, and the error names a longer row length that should fit. SVG is never
  refused. See
  [limits](./REFERENCE/output-formats-and-export.md#png-and-pdf-size-limits).
- **Fixed identity color scale.** `--identity_scale MIN MAX` (Linear and
  Circular), `identity_scale` on `LinearComparisonOptions` and
  `ComparisonRingOptions`, and **Identity scale min (%)** and **Identity scale
  max (%)** in the Web comparison settings fix the identity range that colors
  matches or similarity-ring hits, so figures drawn with different thresholds
  share colors and Legend labels. A match below MIN or above MAX takes the end
  color, and the Legend's two labels show MIN and MAX. `0 <= MIN < MAX <= 100`; a
  blank Web maximum is 100%. Without it, matches are colored from the identity
  threshold to 100% as before. The config keys are
  `objects.blast_match.identity_scale` (Linear) and
  `objects.conservation.identity_scale` (Circular).
- **Italic gene names.** `--italic_gene_labels`
  (`labels.<mode>.italic_gene_labels` in a config file) draws in italics every
  feature label whose text comes from the `gene` qualifier; the Web app has an
  **Italic gene names** checkbox in the labels card. Off by default.
- **`--auto_orient` and the Web button "Orient records from comparisons" for
  Linear diagrams.** `--auto_orient` reverse-complements a record when its
  reverse-strand nucleotide hits against the row above weigh more than half of
  its hits. Also `auto_orient` on `LinearComparisonOptions` and
  `LinearDiagramOptions`. It reads `-b/--blast` files, `--comparisons_table`
  rows, and `--losat losatn` / `tlosatx` results; protein comparisons are not
  read. The first row keeps its orientation, a record that `--reverse_complement
  1` or a region's `:rc` reverses stays reversed, and the log names each record
  it reverses. In the Web app, **Orient records from comparisons** (Linear,
  Comparisons, Record orientation) applies the same rule to the comparisons of
  the last Generate: it sets Reverse complement of every record below the top
  row (checked when the hits lie mostly on the reverse strand, unchecked
  otherwise, also where you set it yourself) as one History step, and the next
  **Generate Diagram** draws it. It is disabled until a Generate has drawn a
  nucleotide comparison. A Session saved with `--auto_orient` stores the
  resulting Reverse complement flags, not the option, and replays the same
  figure.
- **`--anchor_feature KEY=VALUE` and the Web action "Align at feature" for Linear
  diagrams.** `--anchor_feature gene=COX1` moves every record that has the named
  feature so that the feature's first drawn base sits under the first such
  record's. Also `LinearOptions.anchor_feature` and
  `LinearDiagramOptions.anchor_feature`. The selector is an exact `key=value`
  match. A record without the feature keeps its place with a `WARNING:`; a
  selector that matches nothing, or matches features at different starts in one
  record, stops with a message. It is refused with records that share a row, a
  `--wrap_length` that cuts a record, and `--similarity_alignment_feature`. In
  the Web app, **Align at feature** at the end of **Record Layout** takes a
  `KEY=VALUE` selector, lists which records keep their place and which move by
  how many px, and **Apply** writes the record translations as one History step;
  later changes do not align the records again. A Session saved with
  `--anchor_feature` stores the resulting record translations, not the option,
  and replays the same figure.
- **Input summary.** `gbdraw circular` and `gbdraw linear` log one `INFO:` line
  per drawn record and per Linear comparison table. A `WARNING:` line follows when
  a record has no feature to draw, or when the thresholds leave out every row of
  a comparison table; the figure is still written and the exit status stays 0.
  The Web app shows the same warnings under the Generate status line (**Part of
  the input was not drawn.**) and the same lines in an **Input summary** section
  of **Run Info**.
- **`--origin_feature KEY=VALUE` for Circular diagrams (command line and Python
  only).** Rotates each record that has the named feature, such as
  `--origin_feature gene=psbA`, so that the feature's start as drawn is the
  display start. Also `origin_feature` in the Python Circular options. The
  selector is an exact `key=value` match. A minus-strand gene starts at its low
  coordinate and reads counterclockwise; reverse-complement the record to read
  it clockwise from its other end. A join starts on the part its path begins
  with, and a feature across the origin starts on its piece before the origin.
  A record without the feature keeps its start with a `WARNING:`; a selector
  that matches nothing, or matches features at different starts in one record,
  stops with a message. It is refused for a record that already has a display
  start (`--display_start_coordinate`, a records-table `display_start`), a cropped
  record, a record drawn linear, and a feature with no first base in drawn order
  (unstranded, mixed strands, a fuzzy end, `order()`, or overlapping parts). A
  Session saved from such a run holds the resulting `startCoordinate`, not the
  option, and replays the same drawing.
- **Web UI strings.** The comparison notice for a split figure points to the
  `canvas.linear.wrap_length` row under **Preserved session settings**, and the
  new strings are translated in the ten non-English languages.

The drawn figures are the same as before, except where listed above.

## Installation and compatibility

Installation from PyPI and Bioconda is unchanged; see
[Installation](./INSTALL.md). gbdraw 0.14.5 supports the same Python versions
as 0.14.4. 0.14.5 saves Sessions as version 50 and opens Sessions saved by
0.14.4 (version 49) and earlier; gbdraw 0.14.4 and earlier refuse a Session 50
file. See
[Session and request compatibility](./REFERENCE/session-and-request-compatibility.md).

## Known issues

- A Similarity group member popup and the Similarity alignment panel still show
  the coordinates of a CDS across the origin as `1..` the record length. The
  drawing is correct; only this text is wrong.
- With **Auto Reflow** off, on Circular and Linear diagrams, a label text edit leaves
  the label where it is drawn: for a text that **Generate Diagram** places away from
  its feature with a leader line, the label stays without the leader, and a label
  drawn with a leader keeps it for a text that fits, until the next placement of the
  labels or **Generate Diagram**. An Undo of a Generate returns such a label too, also
  when **Auto Reflow** is on by then. 0.14.4 does the same.
- An Undo or a Redo back to **Label visibility** On draws no label on a figure that
  lacks it, until the next redraw or **Generate Diagram**.
- With **Auto Reflow** on, an Undo of a feature hide or of **Label visibility** Off
  changes the figure only when the redraw it starts finishes, not at once. A Redo
  applies at once. 0.14.4 does the same.
- With **Auto Reflow** on, on Circular and Linear diagrams, after you drag the Legend
  and then edit a label text, the Undo that should return the Legend does nothing and
  uses up its step, so the Legend stays where you dragged it. 0.14.4 does the same.
- After you set a stroke on a Legend row, delete the row, press **Generate Diagram**,
  then remove or recaption its color rule or give its feature another Legend name
  (for a row with no color rule, only when its whole row was renamed first and this rename applies to the whole row), and press **Restore**, the next **Generate
  Diagram** fails with "The generated result could not be accepted", on Circular and
  Linear diagrams. A Legend row without a stroke does not fail; set no stroke on a
  row you delete. 0.14.4 does the same.
- Edits made on a Linear File card that draws every record of its file, in a
  Session saved by gbdraw 0.13.0, do not follow its records after the first
  **Generate Diagram** splits the card into one card per record.
- After you choose to turn off a track that other rows reach in the dialog that
  asks first, an Overlay row you turn back on by hand, while its track is back,
  returns in the Result at once, while the note beside the track stack still says
  **Applies on Generate**. The Result then equals the next Generate.
- When two Legend rows share a name, a color rule edit or removal can number them
  in the opposite order until the next **Generate Diagram**. 0.14.3 does the same.
- After you delete a Legend row, giving another color rule that row's name, or a
  new color to a rule of that name, can place the returning row last in the Legend
  until the next **Generate Diagram**. 0.14.3 does the same.

Carried from 0.14.1:

- A Session saved by 0.14.0 with a crop and ring LOSAT still holds unconverted
  rows. Its ring draws as in 0.14.0 on a forward crop, and at another
  unconverted place on a reverse-complemented crop, until the next **Generate
  Diagram**, which reuses the raw LOSAT cache and converts the rows.

Carried from 0.14.0:

- After a palette switch, or the **Legend item** choice on a palette Legend
  row, a color rule in the row's earlier color can show as one Legend row until
  the next **Generate Diagram**, which draws two.

See the [0.14.0 known limitations](./RELEASE_NOTES_0.14.0.md#known-limitations),
which still apply.

[Documentation home](./DOCS.md) | [0.14.4 release notes](./RELEASE_NOTES_0.14.4.md)
