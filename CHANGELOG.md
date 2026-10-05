# Changelog

All notable changes to gbdraw are documented here, in the style of
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/). This project uses
[Semantic Versioning](https://semver.org/); pre-1.0 minor versions
(`0.MINOR.0`) may include breaking changes.

Detailed, per-release notes (migration steps, session/schema compatibility,
and full feature descriptions) live under `docs/RELEASE_NOTES_*.md`. This
file is the short, chronological index; follow the links below for the full
write-up of a release.

## [Unreleased]

- Label whitelist, Label overrides, and Feature visibility tables: a `#` after
  other text in a line is now part of the cell value instead of the start of a
  comment. A Feature visibility row whose value contains `#` (for example
  `foo#bar`) failed with "Missing values", and a Label override text or
  whitelist keyword such as `Gene #1` was silently cut to `Gene `. Only a line
  whose first non-blank character is `#` is a comment, as before.
- Annotation table (`--annotation_table`, `read_annotation_table()`): a `"` is
  now part of the cell value, as in the styling tables; the file is no longer
  read as CSV with quoting. A label that starts with `"` and is never closed
  (such as `"lead`, which the web app downloads as typed) failed with
  "unexpected end of data", and a quoted cell lost its quotes. **Behavior
  change for CLI files:** an annotation table that wrapped a cell in CSV quotes
  (for example to hold a tab) now keeps the quotes in the value; remove them,
  and replace a tab inside a cell with a space, as the web app does.
- Label whitelist and Qualifier priority (web app): a tab or line break typed
  or pasted into a rule cell, such as a whitelist keyword, no longer adds a
  column or row to the table that **Generate** writes. Each run of such
  characters becomes one space and the cell is trimmed, as Feature visibility
  and Label override cells already were. A tab used to shift the cells, so the
  rule was read with the wrong feature type and qualifier, and a line break
  split the rule into two rows and failed with "Missing values".
- Default colors, Specific colors, Qualifier priority, Label whitelist or
  blacklist, Label overrides, and Feature visibility tables: a `"` is now part
  of the cell value; these tables are no longer read as CSV with quoting. A
  cell that starts with `"` used to lose its quotes (`"quoted"` was read as
  `quoted`), and a value that starts with a `"` that is never closed (such as
  `"lead`, which the web app writes as typed) failed with "unexpected end of
  data". **Behavior change for CLI files:** a user table that wrapped a field
  in CSV quotes now keeps those quotes in the value, so remove them from such
  files. The Feature override and Feature placement tables are unchanged.
  files. The Feature override, Feature placement, and Annotation tables are
  unchanged.
- Default colors, Specific colors, Qualifier priority, Label whitelist or
  blacklist, Label overrides, and Feature visibility tables: a row with more
  columns than the table has is now an error that names the file and line
  ("Malformed line ... expected N columns"). Before, a Label whitelist,
  Qualifier priority, or Default colors file whose first row had an extra
  column was read shifted: `CDS`, `product`, `two`, `words` became feature type
  `product`, qualifier `two`, keyword `words`, and no error was raised. Remove
  the extra cells, or the tab inside a value, from such files.
  files. The Feature override and Feature placement tables are unchanged.
- Default colors, Specific colors, Qualifier priority, and Annotation tables:
  a line whose first non-blank character is `#` is now a comment and is
  skipped, as in the Label whitelist, Label overrides, and Feature visibility
  tables. These tables used to read such a line as a row and fail (for example
  "Missing values" or a wrong column count), while the web app import already
  skipped it. A `#` after other text in a line is still part of the cell value.
- Labels (web app): applying **Label visibility** **On** now asks **Feature Is
  Hidden** also for a feature that a **Feature Visibility** rule hides, for
  example a rule with a record ID, another qualifier, or a `hash`, `location`,
  or `record_location` value, as a loaded visibility TSV can give. The dialog,
  the popup note, and the preview after Undo, Redo, or a Result switch decide
  whether a feature is drawn as Generate does, with Python's regular
  expressions; the preview used to apply only the popup's exact product and
  protein ID rules. Such an **On** used to be saved without a dialog and was
  not drawn.
- Linear File order (web app): the File up and down buttons now work when
  each File uses its own consecutive rows, including a CLI Session that draws
  each record of a multi-record file on its own row. A move exchanges the
  File's whole block of rows with the adjacent File's block. Such Sessions no
  longer show "File order is unavailable because Record Layout is custom";
  the notice remains for Files that share a row or whose rows are separated by
  another File's row.
- Feature popup (web app): record rotation adds **Apply on Generate** next to
  **Apply and regenerate**. It stores the previewed start and orientation for
  the feature's record without redrawing, so rotations of several records are
  drawn together by one **Generate Diagram**; each is one Undo step, and the
  rotation section shows **Pending for Generate:** with the staged value until
  then. **Apply and regenerate** still redraws only its own record and leaves
  values staged for other records pending. Closing or switching the popup
  while **Apply and regenerate** runs no longer writes its status into the
  next popup, and a target that becomes stale before Apply states its reason
  once.
- Feature popup (web app): when a feature's input file differs from the one
  the current Result was drawn from, **Rotate record using this feature** now
  says so and asks for **Generate Diagram**. It used to say "The popup feature
  source changed after the popup opened.", which was wrong when the input had
  changed before the popup opened.
- Feature popup (web app): **Edit** now groups its controls by when they
  apply. **Appearance · updates the current Result** holds Fill Color, Stroke,
  Label text and visibility, Feature visibility, and Legend name; **Layout ·
  applies on Generate** holds Feature placement and the closed **Rotate record
  using this feature** section; the Similarity group section follows. The
  header no longer repeats the fill color input or the similarity group, and
  shows `<record ID>: <location>`. Record rotation asks **Put this feature
  at**: **Start of the record** (new default, places the feature first for
  either strand), **End of the record**, or **Custom position** (5′ end,
  midpoint, 3′ end, or just after the feature, shifted by a signed offset).
  **Place this feature at the end** and the Anchor select are replaced; saved
  Sessions are unchanged.
- Feature popup (web app): after **Load Session**, expanding the record
  rotation section reads the feature's source records and shows the preview of
  the new display start, so a record rotates without **Generate Diagram**
  first. It no longer says "The popup feature target is stale or ambiguous."
  Load itself still reads no record bytes. The section shows **Reading
  records…** during the read, a failed read shows its own reason, and a reason
  that leaves no placement available (for example, a non-circular record)
  appears once in place of the controls instead of under each control.
- Gallery: the *Vibrio nigripulchritudo* TUMSAT-TG-2018 Session is now built
  from its declared command and stores its GBFF file as one GenBank resource
  whose six records are selected by record ID, instead of one resource per
  replicon. Load and Save Session keep it as one File; the figure is unchanged.
- Sessions (CLI `--session_output`, `gbdraw.api.save_session_document`): when
  every record of an input file is drawn unchanged, the request reads that
  file's own bytes, the same resource the web app's File binding names, instead
  of a rewritten copy. After **Load Session**, the feature popup's record
  rotation shows its preview without **Generate Diagram** and no longer says
  "The popup feature source changed after the popup opened." Each such file is
  stored once. The *Vibrio harveyi* group and *V. nigripulchritudo*
  TUMSAT-TG-2018 Gallery Sessions now store their GBFF files byte for byte;
  their figures are unchanged.
- Sessions (CLI `--session_output`, `gbdraw.api.save_session_document`): a
  Linear record drawn with `--region`, `--reverse_complement`, or the records
  table `region` or `reverse_complement` column now also reads its input file,
  and the request stores the crop (in source coordinates) and the orientation,
  as Web Save does. After **Load Session** the Linear rows keep the crop and
  orientation, so **Generate Diagram** draws the record the CLI drew instead of
  the full forward record, with the same source coordinates. A reversed record
  rotates without **Generate Diagram**; a cropped record shows "Record rotation
  is unavailable for a cropped record." A `-b` table that touches a reversed
  record is stored unchanged, in its search frame. Circular batch and grid
  records, for which the web app has no per-record crop, and Sessions written
  before this change keep their drawn copy.
- Sessions (CLI `--session_output`, `gbdraw.api.save_session_document`): a
  Linear Session that draws only some records of a multi-record file (for
  example `--record_id`, or a records table that names some of its records)
  now reads that file and selects each drawn record by record ID (by `#n` when
  another record of the file has the same ID), instead of a copy of the drawn
  records. After **Load Session** the file is one File that draws only those
  records, so **Generate Diagram** no longer adds the records the CLI did not
  draw, and record rotation works without **Generate Diagram**. Linear rows of
  a CLI Session now take the record ID from the request instead of `#n`, so
  they no longer show "Selected record was not found in the current file."
  CLI Sessions written before this change load as before.

Fixes from the 2026-09-30 Web GUI audit of `dev`. The plan and the approved
decisions are in
[`docs/internal/web-gui-audit-20260930/`](./docs/internal/web-gui-audit-20260930/03_IMPLEMENTATION_REFERENCE.md).

- Labels: a per-feature override (a label-override `hash` row) now shows or
  hides its label whatever the label display scope selects, in Linear
  (`none`, `first`, `orthogroup_top`) and Circular (`none`) diagrams. In the web
  app, **Label visibility** **On** in the feature popup shows a label on a later
  record under **First Record Only**, or alone under **None**; Generate no longer
  fails with an unclassified render error after such an edit. Multi-record
  Linear label overrides now match their features in the renderer, so **On**
  takes effect there. **Enable Labels** is removed: **On** no longer changes
  **Show Labels** or the label filter. A label text edit on a feature without a
  label asks whether to show the label (**Label Not Shown**).
- Labels (web app): applying **Label visibility** **On** now asks when the
  diagram cannot draw the label: for a hidden feature (**Show feature and
  label** or **Keep feature hidden**), for a feature drawn as **Underlay**
  (**Keep without label**), and, after the label rerender, for a label that does
  not fit with **Label Rendering** set to **Embedded Only** (**Keep without
  label** or **Cancel**). A kept **On** applies when the label can be drawn, and
  Generate no longer fails on it. The feature popup says why a feature has no
  label instead of always suggesting **On**.
- GFF3: a color-table or feature-visibility-table row whose `feature_type` is
  `*` no longer removes the CDS and other features linked to a gene by
  `Parent`. Loading every feature type now flattens those features as a
  type-filtered load does, so the diagram draws the features the feature popup
  lists. In the web app, a manual **Feature visibility** rule with the default
  **Feature Type** `*` removed every such CDS. `read_gff()` without `features`
  now returns these features at the record level too, so `draw_circular()` and
  `draw_linear()` draw them. For a reverse-complemented GFF3 record, the
  feature popup lists features in the drawn start order (OV-15).
- GFF3: a per-feature override that shows a feature of a type the type filter
  dropped no longer parses the GFF3 file a second time; one parse serves every
  type filter, so the case takes about 30% less time. In the web app, changing
  the selected feature types of a GFF3 record no longer parses the file again
  either. The diagram is unchanged.
- Feature placement: a lane placement made in one mode no longer breaks
  Generate in the other mode. Each mode keeps its own placements through mode
  switches, Undo/Redo, and Save/Load Session, and they apply again after you
  switch back. When the current feature slot has no such lane (for example,
  Outward lane 1 after **Track Preset** Tuckin or Spreadout, or Above lane 1
  with **Separate Strands**), Generate names the feature and says to set its
  placement to Auto or Main or to use a slot with that lane, instead of an
  unclassified render error. The CLI and Python API message names the record
  key and biological feature ID (OV-08, OV-09).
- Web Feature placement: a change to **Track Preset**, **Track Layout**,
  **Separate Strands**, or a custom features row's lane or placement that would
  leave lane placements of the current mode without their lane now asks first.
  **Reset N placements to Auto** applies the change and resets exactly those
  placements in one undoable step; **Cancel change** keeps the setting and the
  placements (OV-09).
- Typed requests: canonical `renderRequest` schema 9 adds
  `diagramOptions.featureOverrides`, per-feature Feature visibility, Label
  visibility and label text addressed by record key and biological feature ID,
  and the `featureIdentity` annotation target (`FeatureIdentitySpan`). These
  keep naming the same feature after crop, reverse complement, reordering and
  record duplication, and decide before the visibility and label tables; text
  alone never shows a label. Edits whose feature is not drawn, including exact
  Feature placements whose feature the source does not have, no longer fail the
  render: they are reported as `feature_identity_notices` (Python API, Web
  metadata, one CLI log line each). Schema-8 requests remain readable.
- CLI and Python API: `gbdraw circular` and `gbdraw linear` accept
  `--feature_override_table`, and `CircularDiagramOptions` and
  `LinearDiagramOptions` accept `feature_override_table` (DataFrame) or
  `feature_override_table_file`. Each row names one original-source feature
  with an exact `feature_selector` and sets its Feature visibility, Label
  visibility, or label text, as a `featureOverrides` row does. Run Info's
  **Source recipe** writes a request's `featureOverrides` as this table instead
  of being unavailable; it is unavailable when an edit or placement names a
  feature that the source does not have, which no table row can name.
  Resolving `hash=` rows of feature placement and feature override tables no
  longer reads every feature of the record once per row (4,318 rows of
  `MG1655.gbk`: 145 s to 2.0 s).
- Web per-feature edits: Feature visibility, Label visibility and label text
  edits are kept by the feature's source identity and sent as
  `diagramOptions.featureOverrides` rows, so an edit stays on its feature after
  crop, reverse complement, record reordering and record copies, in the live
  preview, after Generate, and through Save and Load. Session version 45 stores
  them as `features.featureOverrides`; Session 44 and older Web Sessions move
  their rendered-ID edits onto the feature they name, and Load reports how many
  edits it could not match. A Generate that replaces a source removes the edits
  of features the new source does not have; other unmatched edits stay until
  **Remove N unmatched feature edits** (OV-01, OV-02, OV-04, OV-11, OV-12).
- Region Annotations (web app): **Selected features** annotations name each
  feature by its source identity (`featureIdentity` target), so Generate draws
  them on the selected feature after crop, reverse complement, and record
  copies. They used to name the source hash, which Generate compared with drawn
  hashes: on a cropped record the annotation landed on another feature or was
  skipped, and a feature of one copy of a duplicated record failed Generate. The
  editor shows **Selected feature:** and the feature's caption. **Download TSV**
  writes such an annotation as the current Result draws it (`record=#<n>`,
  `feature_selector=hash=<drawn hash>`) and says so (OV-03).
- Web per-feature edits: **Export Feature Edits TSV** and **Load Feature Edits
  TSV** in the Features list write and read these edits as a
  `--feature_override_table`, so a table from the web app runs with the CLI and
  the Python API, and the reverse. Load replaces the edits of the current
  diagram's records as one Undo step, counts rows whose record or feature the
  diagram does not have without applying them, and rejects a malformed table
  with the row it names.
- Comparison tables: every BLAST outfmt 6/7 reader (CLI `-b` and
  `--comparisons_table`, Web uploads, Circular similarity rings, and the LOSATP
  parser) now reads the first 12 columns by position and validates their types.
  Tables with extra columns, such as `-outfmt "6 std qlen slen"`, are no longer
  misread. Only lines that start with `#` are comments, so a `#` inside an ID
  no longer truncates the row. The web app reports a malformed table as a
  comparison-input error with its line number instead of an unclassified
  validation error (CO-05).
- **CLI behavior change:** a missing, unreadable, or malformed `-b` file now
  stops the run with a non-zero exit status instead of being skipped, which
  shifted later tables onto the wrong record pair (N-03).
- Linear comparisons now reject a table whose query or subject IDs name the other
  endpoint or another displayed record: the CLI names the conflicting ID, and
  the web app reports a comparison-endpoint error and keeps the previous
  Result. Unknown IDs keep the positional pair (the CLI logs a warning), and
  SVG record-ID metadata always names the endpoint records (CO-06).
- Python API: `SimilarityAlignmentReference(feature_id="CAG38695.1")` in
  `LinearDiagramRequest.similarity_alignment` or
  `LinearComparisonOptions.similarity_alignment` aligns the Linear records on
  that exact protein ID or feature SVG ID after the requested orthogroup
  analysis, as `--similarity_alignment_feature` does. Python and the CLI share one
  resolver, `gbdraw.api.resolve_similarity_alignment_plan()`, with the same
  errors for a Similarity Group ID, an unmatched ID, or an ambiguous record, and
  Sessions store only the resolved plan. The Python LOSATP tutorial, which still
  passed the removed `align_orthogroup_feature` option, uses it and draws the
  same SVG as the CLI tutorial. `gbdraw linear --similarity_alignment_feature` now
  checks its output paths before the search (OD-2).

- Circular Multi-Record Canvas output (the Web default) no longer reserves an
  empty depth slot when there is no depth input; a one-record canvas now has the
  same slot geometry, legend and center definition as a single-record diagram,
  so every Web-default Circular figure changes (PV-08, N-01).
- A Multi-Record Canvas legend now uses the same builder as a single-record
  diagram: custom slot labels and colors, added skew slots, region annotation
  and depth slot `legend_label` values now appear (TR-01, N-04).
- When the center definition blocks the inside tracks, the species line is
  wrapped at word boundaries and the tracks are placed again. An explicit
  `center_reserved_radius` keeps one line. `definition_font_size` 18 is the
  default, so 18 counts as not explicit and may wrap; any other value is
  explicit and never wraps. The remaining failure names the slot and the reserved definition
  radius, and the Web shows it as `TRACK_LAYOUT` / `DEFINITION_RESERVED`; a
  failure caused by an explicit `center_reserved_radius` names that radius
  instead (`CENTER_RESERVED`) (PV-08, PD-OI-078).
- A specific-color rule whose caption equals a generated legend row, such as a
  rule captioned `rRNA` on tRNA features, keeps its own color row with a hex
  suffix (`rRNA [#ff0000]`) instead of being dropped (N-06).

- Non-pseudo CDS translated without `/translation` now start with `M` when the
  5' end is complete, the reading frame starts at the first base, and the first
  codon is a start codon of `transl_table` (for example `GTG` and `TTG` in
  table 11), as in INSDC `/translation`. This changes **Copy aa FASTA**,
  Interactive SVG feature metadata, and LOSATP and similarity-group protein
  inputs for GFF3 input and for GenBank CDS without `/translation` (FE-07).
- GFF3 CDS phase now sets the reading frame. A CDS with phase 1 or 2 was
  previously translated out of frame, or skipped by the feature popup when its
  length was not a multiple of 3 (N-05).
- LOSATP now skips a CDS whose `codon_start` is invalid instead of translating
  it from the first base. Saved LOSATP rows for a record whose proteins changed
  are not reused; the next search recomputes them.

- The CLI, the Python API, and typed Web or Session requests reject the same
  invalid values and name the option or setting. `--window`, `--step`,
  `--depth_window`, and `--depth_step` take positive integers; before, `0` drew
  empty GC content and skew tracks (X-02).
- A dinucleotide (`-n/--nt`, a track slot's `nt`) is two letters from `A`, `C`,
  `G`, `T`, and `U` in any case. `U` is counted as `T`, so `AU` matches `AT`.
  `-n XY` no longer draws flat tracks and `-n G` no longer raises `IndexError`
  (D-26, N-13).
- Font sizes must be greater than zero and stroke widths zero or greater;
  `--block_stroke_width -1` no longer ends in a traceback. Offsets, spacing,
  `track_axis_gap`, and label rotation keep their current ranges (D-27).
- A Circular definition font set without an interval in a Web or Session
  request, or in Python API `config_overrides`, uses the font size plus 2 as the
  line interval, as the CLI and 0.13.0 do (GE-03).
- Qualifier Priority and whitelist edits apply to Sessions written by the CLI
  or saved from `main`. Label maps compiled by an earlier render are no longer
  kept as preserved settings or reused when a table is attached (SE-08).
- Python validation errors reach the Web with a code, field, reason, and
  Track row, Depth series, table line, or setting path instead of an
  unclassified error (X-01).

- Web validation failures report a recognized code with the field and the
  Sequence, Line, Track row, Depth series, setting, or available band that
  applies, instead of an unknown error. Session import keeps the Worker's
  diagnosis: a file that is not JSON or not a Session names the format problem
  (X-01, SE-09).
- Generate sends numeric settings as typed and no longer replaces a rejected
  value with Auto or a default: a GC window of 0, a decimal step, a non-numeric
  e-value, or a comparison filter outside its range is rejected with the field
  named, and Generate no longer rewrites the settings it reads (X-02: TR-04,
  GE-08, GE-04, CO-09).
- Save before the first Generate explains the missing input with the same check
  Generate uses (IN-08). After loading a Session older than Session 40, Save asks
  for one Generate and offers it; Generate then Save writes the current format
  (D-25, SE-05). The error panel offers **Generate** when that is the correction
  and no **Save Session** after a failed Save (N-12).

- Each Linear record gets a **Definition** inferred from its own `/organism`
  and `/strain`; the first record's value is no longer the file default for all
  records. A definition is the record's own value, then the file default you
  typed, then the inferred value, and **Using file default** marks only a value
  you typed. Choosing another record in a row infers that record's definition.
  Sessions save it as `inferred_definition`; older and CLI Sessions are not
  re-inferred when loaded (IN-06, D-12).
- **Reset Settings** also clears each Linear record's **Definition**,
  **Subtitle**, crop, and reverse complement, and the comparison alignment plan;
  Files, selected records, file defaults, and Depth assignments stay, and
  **Undo** restores all of it, alignment plan included (SE-10, D-15).
- A GenBank record whose ACCESSION or VERSION line is empty (Prokka format) is
  named by its LOCUS instead of `KEYWORDS` in the record choices, the LOSAT FASTA
  extraction, and the match sequences, and an empty ORGANISM line no longer gives
  `Unclassified.` as the organism. The browser reads the header in one place and
  hands any file it cannot read exactly to the Python reader; a file that starts
  with a UTF-8 byte-order mark is rejected as the command line rejects it instead
  of being counted as one record (IN-02, N-02, D-33).
- GFF3 + FASTA record discovery lists only the sequences that have GFF3 rows, in
  FASTA order, like the command line; a FASTA-only sequence is no longer offered
  as a record that Generate cannot draw (IN-03, D-13).
- Replacing a cropped Linear File with a multi-record File expands it into one
  row per record carrying only File-level values, instead of keeping one row
  with the old crop and definition (IN-04).
- A Circular source turned inactive keeps its Multi-Record Canvas record order
  across a switch to Linear and back (IN-05).
- Circular **Region Annotations** offer the records the next Generate draws and
  require a target record whenever it draws several, whether on one canvas or
  as separate diagrams; before, separate diagrams accepted a target the Python
  request could not bind (FE-05).
- **This feature only** always writes the stable feature hash, also when two
  records share a record ID; before, it wrote a rendered ID Python could not
  match. For identical duplicate features the color applies to both (FE-09, D-14).
- On a cropped or reverse-complemented Linear record, a **This feature only**
  color and a legend rename now show in the current Result as Generate draws
  them. Before, the live preview matched the rule against the source feature
  instead of the drawn one, so it kept the old color or dropped the renamed
  legend row until Generate (D-14).
- Hiding a feature also hides its label in the current Result, as Generate
  does; with **Auto Reflow** on, the remaining labels are placed again.
- A legend row that a specific-color rule draws with a hex suffix, such as
  `rRNA [#ff0000]`, now edits that rule in the Legend editor, and a live rule
  edit keeps that row instead of recoloring or failing on the generated `rRNA`
  row (N-06).

- History records one Undo step for each checkbox, radio button, or button
  change, also when it is made with its label text or the keyboard, or while a
  text field has focus. Such changes were previously dropped or merged into
  another step, so edits can now produce more Undo steps (SE-02, SE-03, N-18).
- Ctrl+Z, Ctrl+Shift+Z, and Ctrl+Y (Cmd on macOS) now work while a select has
  focus; text fields keep the browser's text undo (SE-04).
- **Undo** and **Redo** are unavailable while a diagram is generating, and the
  header names the reason. They previously reverted the draft behind a running
  Generate. Settings edits remain available (GE-06).
- Undo or Redo of a step such as a feature color that adds a legend entry, or
  **Reset Settings**, no longer empties the feature list after a mode switch;
  these History steps no longer copy the feature metadata (SE-01, N-19, N-20).
- Cancel during Generate preparation no longer stops the loaded diagram
  engine, so the next Generate reuses it (GE-09).
- Undo, Redo, and the rollback of a failed Session Load restore settings
  exactly. Before, a round trip set unset track slot sides, the Features lane
  direction, and both track axis indexes, pinned an unset Circular
  multi-record legend position, and rewrote the Result's named stroke color as
  hex (F-1).

- A CLI Linear BLAST Session keeps its comparison read-only, and **Inherit
  saved comparison** then Generate succeeds: each Linear file takes the record
  identity of the committed request instead of the CLI binding uid, and the
  saved comparison of an older CLI Session (version 42) is promoted to the
  current request schema before reuse (SE-06, N-17, D-36).
- A CLI Session keeps its `--legend` position through loading and the first
  Generate. The legend position was written into the wrong layout slot and
  replaced by the Web default (SE-07).
- A CLI Session stores the records of one source file as one GenBank resource
  named after that file and selects them by record ID (by index when an ID
  repeats), as the Web does for one multi-record File. The *Vibrio
  parahaemolyticus* and *V. alginolyticus* Gallery Session is now built from its
  declared command and holds its two GBFF files as two multi-record Files instead
  of four single-chromosome Files, so loading it no longer reports a custom
  Record Layout and File order moves are available.
- **Source recipe** is unavailable, with a reason, when a track slot legend
  label contains `,` or ` #` (the CLI would cut or reject it), and for a Linear
  scale font without a ruler-label font while ruler labels are drawn (the CLI
  ruler labels would follow the scale font) (TR-07, GE-03).

- Circular **Species**, **Strain**, plot title text, position, and font size,
  **Keep Full Definition with Plot Title**, and **Default font size** now apply
  on Generate, as in Linear. Editing them after Generate no longer rebuilds the
  record definitions from the whole input file, which replaced the region
  length, GC%, and record label in the Result, Save Session, and exports
  (IN-01).
- The global block, line, axis, and scale stroke colors and widths now apply on
  Generate. A cleared or invalid width no longer reaches the current Result,
  and the global setting no longer overwrites a feature's own stroke before
  Generate. Stroke edits on selected features and legend entries stay live
  (GE-02).

- Label text and label visibility edits are no longer cleared when another
  Result, record, or hidden feature changes the displayed diagram. After
  Generate replaces a source file, only the edits of features that no longer
  exist are removed (FE-01).
- With one Result per record, the Features list and its **Edit** actions follow
  the displayed Result, and its color, visibility, legend, and label edits reach
  the other Results when they are displayed, including their export and saved
  Session (FE-02, FE-03, PV-09).
- Redo, or an unrelated Undo, keeps a feature hidden by **Exact product** or
  **Exact protein ID**; the editor matches these rules like Generate, ignoring
  case (FE-04).
- **Reset fill color** uses the default color of the feature it resets, also
  after the Reset dialog of another feature was canceled (FE-10).
- A label rerender no longer applies settings that wait for Generate, such as
  **Species** or a block stroke width (N-16).

- A legend order made with **Sort** or the move buttons, and a renamed legend
  entry without features such as **GC content**, now survive Generate and Save
  Session; entries that appear later follow the ordered ones (PV-02, PV-03,
  D-08).
- Canvas padding now survives Generate and reaches every Result of a batch,
  without being applied twice (PV-07, D-09).
- Renaming a legend entry with features to the caption of another entry of a
  different color offers **Merge**, **Suffix**, or **Cancel** instead of
  failing with an unrecognized error (PV-04, D-06).
- A Linear **Legend position** change now applies on Generate, like Circular;
  the in-place move did not match the generated layout. Changing it after a
  Generate without a legend no longer raises an error (PV-10, D-30, GE-07).
- Closing the **Editor** with Close or Escape returns keyboard focus to the
  Editor toggle (PV-11), and the **Legend** tab help names the edits it offers
  (PV-12).

- Web LOSAT keeps one job per source file pair, but a record no longer
  searches a database that contains itself unless that search was requested.
  Two records in one file now give the same links as two separate files with
  **Max target seqs** `1`, where self hits previously removed every link.
  Comparisons within one file need more jobs. E-values for multi-record files
  still differ from the CLI, which searches each record pair (CO-04, D-40).
- The Settings LOSAT job count now comes from the plan that Generate runs
  (N-10). Run Info states the E-value database.
- A LOSATP Generate after a Feature visibility change searches again instead of
  reusing comparisons that still included hidden proteins (CO-02).
- A saved Similarity groups or Collinear conversion is reused only for the same
  displayed record pairs in the same direction. A plan change that only flips
  which direction of a record pair is displayed, or a Collinear run with search
  scope **All** and a CLI grid row layout, now converts again instead of
  reusing the previous result. Conversions saved before this change are not
  reused (F-4).
- Threaded LOSAT without cross-origin isolation fails with a LOSAT diagnostic
  instead of UNKNOWN, and the option reads **Threaded (unavailable here)**
  (CO-01).
- Rotating a Linear record with **Orient feature forward** now sets the File
  card's reverse-complement checkbox, so LOSATN ribbons stay on the homologous
  region and the checkbox controls the orientation again (CO-03, N-09).
- A renamed Similarity group keeps its name only on the same members after
  regrouping; other names are saved until those members return (CO-08).
- The web app now shows **Comparison table rows were placed by position.**
  when an uploaded comparison table has sequence IDs that match no displayed
  record. The CLI already logged this warning; the browser discarded it
  (CO-06 carry-over, PD-OI-074). Sessions store the notice in
  `runMetadata.comparisonWarnings`.

- **Changed (command line and Python API):** comparison tables (`-b/--blast`,
  `--comparisons_table`, `linear_comparisons`) are now read in the search
  frame, the selected and cropped record on its source strand. Before, a table
  was read after `--reverse_complement` or a reversed region was applied, so a
  forward BLAST table drew the wrong region on a reversed record. A table
  written for the reversed display needs its coordinates of that record
  converted (`L + 1 - x`). Output changes only for comparisons that touch a
  reverse-complemented record (N-08, CO-07, PD-OI-073).
- A comparison row outside its cropped record stops the run with a
  comparison-input error instead of being drawn beyond the record (N-07).
- **Save Raw LOSAT TSV** of a reverse-complemented record uploads again to the
  same ribbons; the web app no longer converts LOSAT rows to the displayed
  orientation (CO-07).
- Linear match popups and match FASTA headers report input-file coordinates,
  as feature popups do, with the table interval added for a cropped record. The
  SVG record group of a cropped or reverse-complemented Linear record carries
  `data-gbdraw-record-source-start`, `-end`, and `-step` (CO-10, PD-OI-076).
- Sessions saved from `main` with a reverse-complemented Linear record load
  with the same ribbons: the CLI and the web app convert their stored
  comparison rows to the search frame once at Load, and the web app rewrites
  the stored table bytes. A CLI sidecar of `-b` with `--reverse_complement`
  replays the ribbons of the original run. An empty Similarity alignment
  target (`""`) in a `main` Web Session is read as no target. A comparison row
  outside its record reports the new `SEARCH_FRAME` reason of
  `COMPARISON_INPUT` (CO-07, PD-OI-073).

- Web Custom Track Slots: a Depth row is added only when a logical Depth series
  gets its first file and removed when the series loses its last file, in both
  Circular and Linear. Unrelated toggles and **Add Depth TSV series** no longer
  re-create, re-enable, or move Depth rows, and removing a Depth file with the
  uploader **Remove** leaves a stack that generates (TR-02, TR-03).
- Linear: an enabled Depth row whose series has no file shows the same row issue
  as Circular and stops Generate on that row (TR-03).
- Turning off **Hide GC Content** or **Hide GC Skew** while the custom stack is
  off restores the rows that Hide disabled; rows you disabled stay disabled
  (TR-06).
- A disabled track row no longer shows another row's resolved **(auto)**
  geometry (TR-08).
- **Reset to Tuckin**, **Reset to Middle**, and **Reset to Spreadout** follow
  **Show Coordinate Scale**, like **Reset** (TR-09).

- Feature Search and Interactive SVG search: **All** no longer matches
  nucleotide or amino-acid sequences or `/translation` values; use the
  **Nucleotide** and **Amino acid** fields for sequence search. The raw 0-based
  **Start**/**End** search values are removed; **Location** matches the 1-based
  INSDC location (FE-06, FE-11).
- The feature list, feature popup, hover summary, match popup feature
  sections, and Interactive SVG popups show every part of split and
  origin-spanning locations, 1-based, with the summed length (FE-11, N-11).
- A one-feature label, color, or visibility rule no longer spreads to a
  feature whose qualifier value differs only in case (FE-08).
- Specific-color tables accept `none` in any case in the Web app and the CLI.
  Both reject hex colors with alpha and unknown names; the CLI reports the line
  and no longer reads `None`, `NA`, or `null` cells as blank (FE-12).
- Web PDF pages convert CSS px to pt, so they are 75% of their former size and
  match the CLI PDF. Curved and tick labels keep their spaces in the PDF text
  layer (PV-05, PV-06).

- Every visible Web form control has an accessible name. Controls with a
  visible label use that label, including **Window**, **Step**,
  **Dinucleotide**, **GC Content Mode**, and the depth and color-rule fields;
  Custom Track Slots row controls are named by slot id in both modes.
  Placeholders and state-dependent titles are no longer names (TR-10).
- Every **?** help tip is a **Help** button that hover, keyboard focus, a click,
  or a tap opens and **Escape** closes. Tips sit outside labels, and the
  control a tip explains references its text as the accessible description
  (TR-10, PD-OI-057).
- Custom Track Slots rows put the move, duplicate, and remove buttons on their
  own line in both modes, so the slot id and renderer stay readable (TR-11).
- The Vibrio harveyi group Gallery tutorial links to the current Linear record
  layout reference, and a packaging test checks every tutorial link (TR-12).
- Web Run Info and Session request GenBank record counts now use the shared
  reader's record-start rule (a line that starts with `LOCUS` and seven spaces),
  so a malformed `LOCUS` line no longer adds a record the Python loader does
  not read.
- The feature popup **Location** line reads only the location the popup click
  builds with the shared 1-based INSDC formatter; the unreachable fallback that
  rebuilt `start+1..end` from the feature envelope is removed (B7).
- A blank Linear record **Definition** field now shows the text Generate draws
  for it as its placeholder: the File default, else the definition inferred
  from the record, else the example text (B5).
- Web record discovery reuses the settled result for the same source files.
  Generate reports an error the reader already returned for those files without
  a new Worker read, so one upload and one Generate of a rejected GFF3 + FASTA
  pair read it once instead of three times (Linear) or twice (Circular). A new
  upload reads again, and so does a Worker start-up or transfer failure.
- The Web Legend editor's **Add entry** measures the new caption at the DPI the
  displayed Result was rendered with (the committed request's configuration, 96
  by default) instead of 72, so a new entry's width and wrapping match what
  **Generate** draws.
- A Generate or Similarity alignment Apply that fails inside the diagram engine
  without a classified cause is reported as a render failure (`RENDER_FAILED`)
  naming the Python exception class, and offers **Save Session** instead of
  **Retry Generate**, because the same inputs fail the same way. It no longer
  says "Input validation failed".
- PNG, PDF, EPS, and PS exports from the CLI and Python API now place text where
  browsers draw the SVG. CairoSVG ignored `dominant-baseline` on circular tick
  labels, so lower-half tick labels sat about one font height closer to the
  center than in the SVG, and upper-half labels one font descent closer.
  `hanging` and `middle` text (Linear scale and legend labels, record
  definitions, feature labels) also sat slightly off. SVG output is unchanged.
  Gallery thumbnails and the published export examples are regenerated (B8).
- **Load Label TSV** in a Circular batch now applies its rows to the labels of
  every Result, live when each Result is displayed and after Generate, and
  **Undo**/**Redo** treat the import as one step for every Result (B6).
- **Inherit saved comparison** then **Generate Diagram** now works for a CLI
  Linear Session made from one multi-record GenBank file. Each record of that
  file loads as its own Linear row (keyed by its committed record key, selected
  by `#n`), so the saved BLAST ribbons are drawn between the same records again;
  single-record CLI Sessions are unchanged (B3).
- A CLI Linear Session written without `-b` or `--protein_blastp_mode` now loads
  with **No comparison** and without selecting LOSATP, so **Generate Diagram**
  redraws the CLI figure without ribbons instead of starting a LOSATP run (which
  failed with `COMPARISON_INPUT` on records without CDS). CLI Sessions with `-b`
  or a protein mode and Web Sessions are unchanged (B13).
- A display start beyond the record length and a GenBank or GFF3/FASTA file that
  cannot be parsed at render time now show an input error with the actions that
  fix it, instead of the render-failure panel that offers only Save Session
  (B10).
- CLI and Python PNG, PDF, EPS, and PS exports and Gallery thumbnails draw the
  italic and roman parts of a centered or right-aligned caption (for example the
  plot title "*Vibrio nigripulchritudo* TUMSAT-TG-2018, complete genome") side
  by side where browsers draw them; CairoSVG had aligned each part on its own
  width, so the parts overlapped. SVG output is unchanged (B12).
- CLI: `gbdraw linear` with more `-b` files than adjacent pairs of loaded
  records now stops before drawing with `Too many -b/--blast files (expected at
  most N)`. Before, the extra file was silently ignored in the figure and
  written into the Session as a comparison with a record that was not loaded, so
  the Session could not be replayed. With `-b` and two or more input files each
  file contributes one record, so `--gbk multi.gb other.gb -b a.tsv b.tsv` now
  reports the extra table (B14).
- The **Live edit failed** note above the Result names the cause and suggests
  retrying the live edit only when the same request can succeed, such as after a
  Worker failure. A failure that repeats for the same edit, such as
  `RENDER_FAILED` or `INPUT_INVALID`, asks to change the edit, or change the
  settings and use Generate.
- A Linear Session without a stored comparison plan loads with **No
  comparison**. A CLI Session written with `-b` no longer loads with a **Run
  LOSAT for all adjacent pairs** replacement draft, so **Replace with current
  controls** waits until a comparison is set up; **Inherit saved comparison**
  still reuses the CLI ribbons. A CLI Session written with
  `--protein_blastp_mode`, including a 0.12 or 0.13 sidecar, keeps the adjacent
  LOSATP comparison it drew.
- In a Circular batch, **Undo** and **Redo** of a **Layout edit** drag or a
  position reset restore the Result it was made on, also while another Result is
  displayed, and keep the displayed Result. Before, the Undo wrote the dragged
  Result's positions into the displayed Result. Displaying another Result no
  longer records an **Undo** step (B17, B18).
- In a Circular batch, a legend entry that only one Result draws (for example a
  **This feature only** color entry) keeps its place after a legend **Sort** or
  **Move** when another Result is displayed and that Result is shown again
  (B18).
- In a Circular batch, **Undo** and **Redo** of a legend step made on another
  Result no longer copy a legend entry that only that Result draws into the
  displayed Result, and the restored legend order reaches each Result when it is
  displayed. **Sort by default** now also reaches a Result that was shown in an
  earlier sorted order (B19, B20).
- Adding a Circular **Pairwise Comparisons** ring file with **Add Seq**, or
  choosing an uploaded BLAST row's **Comparison sequence (optional)**, is now
  one **Undo** step. Before, these hidden file inputs recorded no step, so
  **Undo** could not remove the added file (B21).
- In a Session whose Circular rings replay saved LOSAT rows, **Undo** and
  **Redo** of **Remove series** or **Add Seq** now restore those rows with the
  ring rows. Before, with **Use custom stack**, the next Generate silently drew
  one ring fewer or failed with a track-settings error (B23).

LOSAT CLI/API: the CLI and Python API run LOSATN, TLOSATX, and LOSATP directly
and share the Web app's raw search keys. The design and the approved decisions
are in
[`docs/internal/LOSAT_CLI_API_DESIGN_PROPOSAL_2026-10-03.md`](./docs/internal/LOSAT_CLI_API_DESIGN_PROPOSAL_2026-10-03.md).
Retired names and their replacements are listed under
[Retired inputs](./docs/SESSION_COMPATIBILITY.md#retired-inputs).

- **Breaking:** LOSATP CLI flags and Python fields were renamed. Use
  `--losat losatp` with `--losatp_mode similarity_groups|collinear|pairwise`
  (replaces `--protein_blastp_mode`), `--losat_bin`, `--ncbi_blast_bin`,
  `--losat_threads`, `--losatp_max_hits`, `--losatp_max_target_seqs`,
  `--similarity_alignment_feature` (replaces `--align_orthogroup_feature`), and
  `--losat_output_dir` (replaces `--protein_blastp_output`); `--losatp_member_max_hits`
  and `--collinear_infer_orthogroups` are new. Retired flags exit with status 2
  and name the replacement, 0.13.0 Session argv is rewritten on replay, and the
  typed LOSATP settings move to `LosatSearchOptions` and `LosatRuntimeOptions`
  (`protein_mode`, `blastp_executable`, `candidate_limit`, and
  `orthogroup_member_max_hits` are replaced as listed under Retired inputs).
  Persisted Session names are unchanged. Docs and Gallery commands use the new
  names (#730).
- Added: `gbdraw linear --losat losatn|tlosatx` and the Python
  `LosatSearchOptions(program="losatn"|"tlosatx")` /
  `LinearComparisonOptions(losat=...)` run LOSATN and TLOSATX directly, with
  `--losatn_task`, `--losat_gencode`, the records-table `losat_gencode` column,
  the comparisons-table `source` column (mix searched and uploaded edges), and
  `--losat_output_dir` raw TSVs plus a reusable `comparisons.tsv`. Results match
  the web app's search and raw cache keys, and saved Sessions replay without
  LOSAT (#731).
- Added: Circular similarity rings can run LOSATN or TLOSATX from the CLI
  (`gbdraw circular --losat losatn|tlosatx --conservation_sequence ...`) and the
  Python API (`ComparisonRingOptions(losat=...)`,
  `CircularDiagramOptions.losat_search`). Comparison genomes may be FASTA,
  GenBank, or DDBJ (#732).
- **Breaking:** `--conservation_sequence` replaces `--conservation_fasta`, the
  `--conservation_table` column `comparison_sequence` replaces `comparison_fasta`,
  and `CircularDiagramOptions(conservation_sequence_files=...)` replaces
  `conservation_fasta_files`. The retired flag and column are rejected with the
  replacement named, and the Python field has no alias (#732).
- Added: native LOSAT runtime handling for LOSATN, TLOSATX, and LOSATP lives in
  one owner (`gbdraw.comparisons.losat_runtime`). CLI Sessions record the runtime
  (kind, version, source, path, program, CLI dialect) in each new `losatCache`
  entry (#726).
- Changed: CLI and Python LOSATP searches use the Web source-file database scope
  (one file is one genome; a record never searches itself unless requested), so
  raw cache keys equal the Web keys. LOSATP E-values change only for records that
  come from a multi-record file (#733).
- Fixed: CLI and Python Sessions store the records of one multi-record source
  file in one GenBank resource with record-index selectors (the Web layout), so
  replay and the Web app reuse their LOSATP raw searches. Rendered output is
  unchanged (#733).
- Fixed: `gbdraw linear --session` replays the 0.13.0 Gallery Session for
  BGC0000708-BGC0000713 again; similarity alignment members saved for
  reverse-complemented records now bind to their source features (#734).
- Added: Web Circular similarity rings accept GenBank and DDBJ comparison files.
  The diagram worker reads them with the Python reader of
  `--conservation_sequence`, so a Web ring has the CLI raw cache key for every
  format and FASTA layout. The ring row's "Subject gencode" is now "Comparison
  gencode" (#736).
- Added: Web LOSAT raw cache entries record the search runtime (`wasm`), and Run
  Info lists the search runtime of each displayed result, including the path that
  a CLI Session recorded (#736).
- Fixed: default Web Linear raw LOSAT TSV names replaced the letter "s" instead
  of whitespace; they now follow the CLI naming rule (#736).
- Fixed: `gbdraw ... --session ... --save_session` keeps nucleotide BLAST and
  LOSATN/TLOSATX comparisons as `nucleotideBlast` items with their table bytes
  instead of re-saving them as protein comparison tables. It also no longer fails
  when two input files share a name: ring files with one stem (`X.fna` and
  `X.gb`) or files with one basename in different directories
  (`--conservation_sequence a/X.fna b/X.fna`, `--conservation_blast`, `-b`). The
  later file is saved as `X.2.fna` or `X.circular_conservation.<program>.2.tsv`,
  the numbering `--losat_output_dir` already uses (#738).
- Fixed: Web Circular similarity ring rows added from a GenBank or DDBJ file
  without a typed label are labelled with the first record's DEFINITION, or the
  organism when the DEFINITION is empty, like `gbdraw circular
  --conservation_sequence` and the Python API. FASTA rows keep the file name
  without its extension (#740).
- Fixed: a Web Session with two Circular LOSAT rings of one sequence now loads.
  The Session stores one LOSAT cache entry per raw key, as the CLI does; a Session
  that repeats a key reports a Session diagnostic instead of an unknown error
  (#743).
- Fixed: Redo of a Web Circular similarity ring row added with **Add Seq** from
  a GenBank or DDBJ file restores its DEFINITION or organism label instead of
  the file name. The label read is part of the add's one History step; Undo,
  Redo, and Session save and load are unavailable until it answers (#747).
- Fixed: in a standalone interactive SVG, the popup of a Linear collinear block
  with several anchors now copies and downloads its Query span, Subject span,
  and both spans. It showed "Match feature endpoint identity is invalid." under
  each span and offered no sequence actions; a one-anchor block was not
  affected. Gallery SVGs carry the fix after their next refresh.
- Fixed: in a standalone interactive SVG, the popup of a Linear collinear block
  with several anchors lists, for each Similarity group, only the anchors that
  belong to that group under Query member and Subject member, and lists one row
  per anchor under Query and Subject. It listed every anchor of the block for
  every group and joined the anchors into one row; the Web popup was not
  affected. Gallery SVGs carry the fix after their next refresh (OV-18).

## [0.14.0](./docs/RELEASE_NOTES_0.14.0.md)

- Added circular-record display-start rotation and manual feature lane placement.
- Improved Run Info / Exact replay, Save Session resource preservation, preview
  navigation, and Generate processing-stage feedback.
- Includes the beta's package-root Python API and current session/request
  compatibility, plus isolated wheel/sdist installation and local GUI packaging fixes.
- See the [full release notes](./docs/RELEASE_NOTES_0.14.0.md) for migration,
  compatibility, and installation availability. Publication dates are recorded
  in [GitHub Releases](https://github.com/satoshikawato/gbdraw/releases).

## [0.14.0b0](./docs/RELEASE_NOTES_0.14.0b0.md) — unreleased (beta)

- Added a small top-level Python interface (`read_genbank`, `read_gff`,
  `draw_circular`, `draw_linear`, mode-specific `CircularOptions` /
  `LinearOptions`, and a first-party `Diagram` result) alongside the
  existing typed `gbdraw.api` request/session/table contracts.
- One Circular function now handles both single- and multi-record input.
- Removed obsolete low-level convenience re-exports from `gbdraw.api`.
- See [the full release notes](./docs/RELEASE_NOTES_0.14.0b0.md) for the
  complete list of changes, including architecture/API and session-format
  updates.

## Earlier releases

Releases before 0.14.0b0 predate this changelog and were not recorded with
per-version release notes. Their tags and dates are listed below for
reference; see `git log <tag>` or the
[GitHub tag list](https://github.com/satoshikawato/gbdraw/tags) for the
commits each one contains.

| Version | Date |
| --- | --- |
| 0.13.0 | 2026-07-05 |
| 0.12.1 | 2026-06-27 |
| 0.12.0 | 2026-06-26 |
| 0.11.0 | 2026-05-07 |
| 0.10.0 | 2026-04-29 |
| 0.9.2 | 2026-04-08 |
| 0.9.1 | 2026-04-06 |
| 0.9.0 | 2026-04-06 |
| 0.8.0 | 2025-12-18 |
| 0.7.0 | 2025-10-26 |
| 0.6.0 | 2025-10-09 |
| 0.5.3 | 2025-09-29 |
| 0.5.2 | 2025-09-09 |
| 0.5.1 | 2025-09-08 |
| 0.5.0 | 2025-09-01 |
| 0.4.0 | 2025-08-07 |
| 0.3.0 | 2025-07-24 |
| 0.2.0 | 2025-05-25 |
| 0.1.1 | 2025-05-18 |
| 0.1.0 | 2025-05-14 |
