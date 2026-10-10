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

- Sessions (CLI and Python): a Session saved from a Linear diagram with a
  reverse-complemented record now stores that record's feature-bound comparison
  rows (saved or fresh LOSATP rows) in the search frame, as the web app does,
  also with `--similarity_alignment_feature`. Before,
  `gbdraw linear --session ... --session_output` and Python Session saves wrote
  them after the reverse complement with extra `*_view_feature_svg_id` columns;
  the drawing was the same (OV-399).
- Gallery publication: a Session saved in the web app after Load, without
  **Generate Diagram**, is accepted when its committed request states a default
  that the web draft leaves implicit: a definition-line `font_weight` of
  `normal` or a feature rendering at its default (OV-398).
- Comparison (web app): **Generate Diagram** in Linear mode with pairs set to
  **Upload BLAST TSV** but without a file opens the **BLAST TSV missing** dialog
  instead of failing. It lists the pairs and offers **Choose BLAST TSV for
  #i → #j…**, which opens the first pair's file chooser, **Set to No comparison
  and Generate**, one undoable step, and **Cancel**. Before, Generate failed with
  "Choose a BLAST TSV for this pair, or set the pair to No comparison or Run
  LOSAT." Another problem in the same pair plan still reports its own error.
- Sessions (web app): **Generate Diagram** on a loaded Linear Session whose saved
  comparison waits for **Inherit saved comparison**, **Replace with current
  controls**, or **Clear comparison** now says "Choose how to handle the saved
  comparison in the Comparison panel, then Generate again." Before, it showed the
  generic "The operation failed without recognized diagnostic information" (OV-303).
  A saved comparison whose resource is missing still reports the comparison input
  error.
  After **Replace with current controls**, an incomplete pair in the comparison
  controls reports its own comparison input error, and controls that set no
  comparison say "Set up a comparison in the Comparison panel, or choose Clear
  comparison, then Generate again." Both showed the generic error before (OV-301).
- Sessions and Legend (web app): a stored Legend entry color and a feature color
  override are checked with the same rule as the **Override File (-d)** import. A
  Legend entry whose stored color is outside that rule, such as `ButtonFace`, is
  dropped at Session Load and named in the Load notice. Before, Load kept it, and
  when the Result did not draw that entry, every **Generate Diagram** failed with a
  generic error (OV-300). A stored color
  name is now saved as its table hex, such as `#FF0000` for `red`.
- Default colors (CLI, Python API, and web app): the web app's **Override File (-d)**
  import accepts the documented Default colors forms, including `rgb()`, `rgba()`,
  `hsl()`, and `hsla()` rows that it dropped before (OV-302), and rejects any other
  value with its line number. The web app no longer asks the browser to resolve a
  color name outside the color-name table that `gbdraw` uses, such as the system
  color `ButtonFace`: where it needs a hex color it shows its default swatch, and
  **Generate Diagram** reports the invalid color (OV-272). Compatibility: the CLI
  and the Python API no longer accept `currentColor` and `inherit` as a user color
  (`-d` rows, configuration override colors, depth track colors, and the CLI stroke
  and label color options). Through 0.13.0 they accepted both without documenting
  them. They also no longer accept the other values of svgwrite's paint type as a
  user color: a paint reference such as `url(#id)`, `icc-color()`, and an empty
  configuration color; the web app never accepted them. A blank `-d` cell and an
  unset depth track color still keep their defaults.
- Sessions (web app): a Linear Session that the Python API saved from a CLI LOSATP
  comparison loads with the comparison plan the CLI drew, so the next **Generate
  Diagram** draws the same comparison. Before, the web app took the file bindings
  and the comparison plan of a Session without a Web draft from the request only
  when the CLI had written it, so the Python API's copy generated without its
  comparison (OV-296).
- Sessions (web app): after **Generate Diagram** on a Linear Session whose record rows
  come from one input file (a GenBank file with several records), every row still
  reads its file. Before, every such row but the last read as empty once Generate had
  handed the shared file to the renderer, so later steps that read those rows' file
  saw no content (OV-304).
- Gallery tooling: `tools/publish_gallery_session.mjs prepare` accepts the Session that
  `gbdraw circular` writes, a `gbdraw linear -b` Session (its read-only comparison is
  inherited, as **Generate Diagram** does after **Inherit**), and CLI Sessions with `-d`
  or the visibility, label, whitelist, and qualifier priority tables. The published
  Session keeps those tables when the web app loads it (OV-266, OV-267, OV-268, OV-299).
- Sessions (CLI, Python API, and web app): a Session 40 or 41 file that an unreleased
  `main` build of the CLI wrote opens from `--session`, the Python API, and the web app
  with its tables, and the next **Generate Diagram** draws the CLI figure. Before, the
  web app rejected these files or dropped their tables, and the earliest Session 40 files
  also failed in the CLI and the Python API (OV-269, OV-273). Session 40 files that the
  web app of those builds saved also open; they did not record moved legend, title, or
  scale positions, so the next **Generate Diagram** lays those out automatically.
- Sessions (web app): a Session that the CLI or the Python API writes keeps its
  color, feature visibility, label, whitelist, qualifier priority, and feature
  override tables when the web app loads it, so the next **Generate Diagram**
  draws the figure the Session was saved with. Before, Load dropped those
  tables and Generate drew a different figure (OV-221).
- Labels (web app): the Labels panel shows the label text settings (**Label Font
  Size**, **Label Rendering**, placement, rotation, and spacing), **QUALIFIER
  PRIORITY**, and the Circular **LABEL GEOMETRY** while any feature of the shown
  drawing has **Label visibility** **On**, also when **Show Labels** (Linear) or
  **Label Mode** (Circular) is **None**. They were hidden there although the On
  labels use them. Label filtering keeps its condition: an On label ignores the
  whitelist and blacklist (OV-222).
- Specific colors (web app): a Specific-colors table imported in the web app accepts the
  color names that `gbdraw -t` accepts, in any letter case, and rejects the others with
  the line number. Before, the browser also read `currentColor` and system colors such
  as `ButtonFace` and imported them as a browser-chosen hex color (OV-271). In the CLI,
  a mixed-case name such as `DarkGrey` in `-t` or `-d` no longer stops drawing with an
  svgwrite `TypeError` (OV-270).
- Custom Track Slots (web app): stepping a Depth row's track index past the loaded
  Depth series (for example ArrowUp then ArrowDown) no longer adds a series that made
  Generate and **Save Session** fail with "Depth series 3"; a label edit on a row whose
  index names no series changes nothing (TK-03).
- Color Change Scope (web app): **Cancel** closes the dialog at once and records no
  Undo step, also right after a Session load, when it used to wait seconds (OV-161).
- Output Prefix (web app): the typed **Output Prefix** reaches Python as typed, so
  `../../x` fails Generate with an error that names **Output Prefix** instead of
  drawing a Result named `.._.._x.svg` that the browser saved as `_.._x.svg`. A
  Result named from a record ID is named as the browser saves it: a record such
  as `gi|1|ref|X` draws `gi_1_ref_X.svg` instead of failing Generate (FL-10).
- Layout (CLI and Python API): Legend, multi-record and definition coordinates no
  longer depend on the Python version. The side-by-side rows of a Circular Legend, the
  multi-record grid, the rows of a multi-record Linear diagram, and the definition text
  block summed floats with the built-in `sum()`, which Python 3.12 made compensated, so
  Python 3.10 and 3.11 placed them up to about 1e-13 px differently from 3.12 and later,
  and from the web app.
- Legend editor (web app): every Legend edit (adding, removing, renaming, or moving a
  row, Sort, a side change, Undo and Redo) lays the Legend out as Python does, with the
  same text measurement and the font, size, and wrap width the diagram was drawn with,
  so the Legend and the canvas you see equal those of the next Generate, in Linear and
  Circular. Before, the web app used its own layout: rows could overlap after an edit
  (for example an added row and a renamed row in Circular), Generate could move rows
  and change the canvas size, and in Linear a renamed row that was not the last one
  moved at Generate.
- Legend editor (web app): Reset Stroke, Reset all strokes, and the **Auto** stroke
  color give each feature the stroke Generate draws when the diagram has an automatic
  underlay (for example `repeat_region` in Circular), and Reset all strokes returns
  each Legend row to the stroke it was drawn with. Before, these gave the features no
  stroke there, because the default stroke was read from the underlay, which is drawn
  without one, and Reset all strokes gave every Legend row that default stroke, GC
  rows included.
- Specific color rules (web app): a rule added or edited while an automatic
  rerender is updating the diagram is applied once the rerendered diagram is
  ready. Before, it could be dropped without a message when the diagram finished
  updating during the edit.
- Color names (CLI and Python API): `seashell` now resolves to `#FFF5EE`, the CSS
  color, so the CLI draws it as the Web app does. Before, it resolved to `#2E8B57`
  (seagreen). `rebeccapurple` (`#663399`) is now accepted. The other 146 CSS color
  names were already correct.
- Color names (web app): the web app resolves a color name with the same CSS table
  as Python, also where the browser offers no canvas, so a named Legend color or
  stroke in a Session (for example `gray`, the stroke an SVG had before an edit)
  loads as its hex value instead of being dropped (OV-160). An unknown name is
  still dropped.
- Sessions (CLI and Python API): a Session that a CLI run wrote with `--session_output` or
  `--save_session` no longer replays an empty label text as the label "nan". A label
  override row with an empty label text hides the label when the written Session is
  replayed, as it does in the first run, and the row keeps its empty text. The same
  reader now keeps an empty cell or the text "NA" in a feature visibility, label
  whitelist, qualifier priority, color, or annotation table as written; before, these
  cells were read as missing values.
- Text measurement: label, Legend, and title widths now apply the bundled fonts'
  GPOS pair kerning, as browsers draw the text. Before, only the older `kern` table
  was read. The two hold the same pairs for Latin, Greek, and Cyrillic text, so those
  diagrams do not move. Hebrew text in Liberation Sans (the default font, also used
  for `Arial`, `Helvetica`, and `sans-serif`) is now kerned, so its labels and Legend
  rows can shift by up to about 1 px per kerned letter pair at 14 pt.
- Text measurement: repeated spaces, tabs, and line breaks in labels and captions, and
  characters the bundled fonts lack (for example Japanese), are now measured as Chromium
  draws them: collapsed into one space, and one em wide instead of nothing.
- SVG output: the Legend part (`legendReflow`) of the `data-gbdraw-composition`
  metadata also records the font file, font size, DPI, and wrap width the Legend was
  laid out with, and the diagram bounds it was placed against, so that the web app
  can lay out an edited Legend as Python does. Drawn content is unchanged.
- CLI: `gbdraw render --session FILE` renders the drawings of a saved Session:
  every drawing with a committed render by default, or the ones named with
  `--drawing ID` (an ID or a name, repeatable). A drawing without a committed
  render is skipped with a notice, and naming one is an error.
  `--list_drawings` prints each drawing's ID, mode, name, and whether it has a
  committed render. One drawing keeps its output names; several write
  `<prefix>_<ID>`, and every diagram path and the Session output are checked
  before the first file is written. `--save_session` or `--session_output`
  writes the Session again with only the rendered drawings replaced; the
  Session is brought to the current version and validated before any render,
  so a Session that cannot be written again leaves no diagram behind. For a
  Session 31–39, which saved no feature catalog, the command prints a warning
  that names each saved Result it drops; the re-saved Session holds the Results
  of the new render.
  `gbdraw circular --session` and `gbdraw linear --session` take `--drawing ID`
  and, without it, render the Session's only drawing of their mode; a drawing of
  the other mode is an error that lists the drawings.
- Python API: a Session's diagrams are drawings. `SessionDocument.drawings`
  lists them as `SessionDrawing` values (ID, name, mode, and whether a committed
  render exists); `SessionDocument.drawing()` selects one by ID or name, and
  `SessionDocument.active_drawing_id` names the one the Web app opens. A
  selection that is missing, unknown, or ambiguous raises
  `SessionDrawingSelectionError` with the list of drawings. Pass `drawing=` to
  `session_to_request()` and `render_session()`. The new
  `render_session_drawings()` renders several drawings together: several
  drawings write `<base>_<id>` outputs, every output path is checked before the
  first write, and each embedded resource is parsed once for all drawings.
  `build_session_document()` and `save_session_document()` take
  `drawings=[...]` (typed requests or `SessionDrawingSpec` values) and
  `active_drawing=`; the undocumented `adjunct=` argument is replaced by
  `SessionDrawingSpec(state=...)`. The new `upgrade_session_document()` returns
  a Session 31–44 in the current version without rendering it, as a
  `SessionUpgrade` with the `document` and its `warnings`. Sessions 31–39
  saved no feature catalog, so the upgrade drops their Results; each such
  drawing gets a warning that names every dropped Result, also logged.
  Session 46 holds at most one Circular and one Linear drawing, named by their
  mode.
- Python API: `gbdraw.api.derive_region_drawing()` derives a new drawing of
  selected regions (`RegionSelection`) from a materialized Session. It keeps the
  drawing's look for those records, refuses regions that cross the origin, and
  returns a size you set to Auto only when its Auto value changes at the new
  length and mode (`adaptation.reset` lists each one).
- Legend colors (web app): Generate no longer fails with "The generated result could
  not be accepted" after a color or stroke is set on a row added in the Legend editor,
  in Linear and Circular. The added row is drawn with its color and stroke after
  Generate, a rename, and Session Save and Load. At Generate, an added row no longer
  takes the edited stroke of the first Legend row.
- Legend editor (web app): a row added in the Legend editor now has the same stroke
  when it is added and after Generate: the stroke of the first Legend row as the
  diagram drew it. Before, the row took that row's edited stroke, and the Block
  Stroke settings before Generate, until the next Generate drew the drawn stroke.
- Legend editor (web app): Generate draws the Legend with a row added in the Legend
  editor where the add placed it. Before, Generate moved the added row away from the
  other rows, in Linear past the right edge of the canvas, where its caption was cut.
- Legend editor (web app): a stroke set on a Legend row stays on the row's features
  after Generate, in Linear and Circular, and Undo, Redo, and Reset of that stroke
  change the features as Generate draws them. A feature with a stroke of its own from
  the feature popup keeps it. Before, Generate drew the row's features without the
  row's stroke, and the row's stroke replaced a feature's own stroke until Generate.
- Legend editor (web app): after a Legend row is removed, Generate draws the Legend
  and the canvas as the removal left them, in Linear and Circular: the other rows
  close the gap. Before, Generate kept the gap and the earlier canvas size. A batch
  Result shown after the removal is laid out the same way.
- Legend editor (web app): removing a row, adding a row, changing a row's stroke
  width, and Reset Stroke are each one Undo step, as the stroke color is. Undo and
  Redo of a removed or added row return the Legend, its strokes, and the canvas as
  they were. Before, these edits made no Undo step when they were not made through
  a control, adding a row never made one, Undo of a removal left the rows out of
  place, and an added row that Undo returned took the first row's stroke edit.
- Legend editor (web app): rows added or removed in the Legend editor keep their
  layout through an automatic rerender (for example after a Feature Visibility
  rule), as Generate draws them. Before, the rerender moved an added row, in Linear
  far to the right of the other rows, and kept the earlier canvas size.
- Legend editor (web app): after a row without features, such as GC content, is
  renamed in the Legend editor, Generate draws the Legend and the canvas as the
  rename left them, and the rows keep their spacing. Before, Generate moved only the
  renamed row, which left uneven gaps, and kept the earlier canvas width.
- Legend colors (web app): after **Apply to all** in the feature popup colors every
  feature of a Legend row, removing the rules it wrote (**Clear All** in SPECIFIC
  RULES, or deleting the last of them) also removes the Legend color it stored for
  the row, in the same Undo step. The row then takes the palette color its features
  return to, live and after Generate. Before, Generate drew the row in the popup's
  color, and the saved Session kept it. A Legend color set on a row in the Legend
  editor is kept.
- Legend editor (web app): opening or closing a row's **Stroke options** no longer
  makes an Undo step, and Sessions no longer save which rows show their stroke
  options. Loading a Session closes them.
- Legend editor (web app): rows removed in the Legend editor are listed under
  **Deleted items**, each with **Restore**, and **Restore all** returns every one. A
  restored row returns at once where Generate draws it; each click is one Undo step,
  and Sessions save the rows still removed. Before, only Undo returned a removed row.
- Legend editor (web app): a renamed feature row (for example `tRNA`) keeps its place
  in the Legend, live and after Generate, in Linear and Circular. Before, the rename
  moved the row to the end of the Legend.
- Legend editor (web app): each row's name field, its move, **Stroke options**, and
  remove buttons, and the stroke color, stroke width, and reset controls inside
  **Stroke options** now have accessible names that include the row's caption (for
  example "Move CDS up", "Remove CDS", and "Reset stroke of CDS to default"), so
  assistive technology tells the rows apart. Before, the name field and the stroke
  width had no name of their own and every row's buttons shared one name.
- Sessions (web app): a Session saved before Merge was limited to two rows of one
  feature type (for example one that merged `GC skew (+)` into `GC skew (-)`) loads
  as saved and keeps its merged row; a new Merge between such rows is no longer
  offered.
- Result names (web app): a live edit that redraws a loaded Session's Result (such
  as a Legend color, or Undo or Redo of a color step) keeps the Result's saved name.
  Before, the redraw renamed it after the **Output prefix** (a Gallery Session's
  `HmmtDNA_basic_circular` became `out.svg`), so **SVG**, **PNG** and **PDF**
  downloaded `out.*` and **Save Session** saved the new name. Generate still names
  its Results after the **Output prefix**.
- Legend names (web app): Generate no longer fails with "The generated result
  could not be accepted" after a **Depth** row renamed in the Legend is hidden by
  **Show Depth**, in Linear and Circular. The row is not drawn; as for a renamed GC
  row that is switched off, it returns under its series caption when Depth is
  shown again.
- Legend names (web app): a Legend rename of a row that only track data names
  (a Depth series, an annotation set legend label) is now retired together with
  that data, with the color and stroke stored under the new name. Generate no
  longer fails with "The generated result could not be accepted" after the Depth
  file or track of a renamed row is removed, or the file is replaced by one with
  another label. A replacement that keeps the row's caption keeps the rename, and
  Undo restores the data, the name, and the styles.
- Depth tracks (web app): a Depth series that no record gives a Depth TSV now
  fails Generate with "The depth input or settings are invalid. Depth series N.
  Supply the required value." and the **Depth** actions. The same failure showed
  "The operation failed without recognized diagnostic information".
- Feature strokes (web app): a stroke set on one feature in the feature popup is no
  longer lost when Generate draws a Result without that feature, such as a Generate
  in the other mode (Circular or Linear). The stroke stays in the draft and in a
  saved Session, and the next Generate that draws the feature draws it again.
- Web: Load of a Session 44 or older whose Feature visibility edit holds a value such
  as `constructor` (a hand edit; no writer saves one) drops that edit and reports it
  with the other dropped edits, as the CLI does. Before, Load failed with an invalid
  draft. A Session 41-44 Feature placement draft keyed `__proto__` is now rejected by
  Load like any other key that does not name its feature, as the CLI rejects it.
  Before, Load lost the row without a message.
- CLI: `gbdraw circular|linear --session <Session 40-44> --session_output out.json`
  (and `--save_session`) moves an annotation's `hash=` target to the feature's
  source identity where the Session's saved feature catalog makes the figure
  certain, as the Web app does on Load. Before, the rewritten Session kept the
  `hash=` target in `config.annotationSets`, so after a later crop or reverse
  complement the annotation no longer matched its feature. The replay logs how many
  targets moved.
- CLI: `gbdraw circular|linear --session <Session 40-44> --session_output out.json`
  (and `--save_session`) keeps the Feature visibility, Label visibility, and label
  text edits that the Web app saved by rendered feature ID. Before, the rewritten
  Session dropped them: the diagram kept their effect through the request's tables,
  but the Web app no longer listed them as edits. The replay now moves each edit to
  `features.featureOverrides` through the Session's saved feature catalog, as the
  Web app does on Load, and logs how many edits it dropped or now applies to fewer
  features.
- Labels (web app): **Load Label TSV** with a table that applies to no label of the
  displayed Results no longer removes the existing label edits. Before, it cleared the
  bulk and per-feature label edits, recorded a "Load label edits" Undo step, and
  reported "Applied to 0 label(s)." It now changes nothing, records no step, and says
  that no row matched and that the existing label edits were kept.
- Sessions: a Session whose resource declares `checksum` (the SHA-256 digest of its
  bytes, as `sha256:<hex>` or bare hex) now loads in Python and on the command line,
  as it already did in the web app. `load_session_document()` and
  `gbdraw circular|linear --session` failed with "has unknown field(s): checksum",
  also for such a Session after the web app saved it again. A resource whose bytes
  do not match its `checksum` is rejected when the Session is loaded.
- Diagram modes (web app): Circular and Linear each keep their own Result. A
  Generate replaces only the Result of its own mode, so a Linear Generate no
  longer replaces the Circular Result, and switching back shows the Circular
  Result with its moved Legend, title, and live edits. A mode without a Result
  shows "No Circular Result yet" or "No Linear Result yet". A moved Legend or
  title on one mode's Result no longer makes the other mode's Generate fail.
  The mode buttons wait while Generate, a label update, or an Undo or Redo
  runs, and an Undo or Redo of a switch waits for Generate or a label update.
- Diagram modes (web app): Circular and Linear each keep their own settings and
  edits: form and advanced settings, palette and color rules, filters, Legend
  edits, feature and label edits, Feature placements, annotation sets, Depth
  series, custom track stacks, and layout positions. A change made in one mode,
  and its Generate, no longer reach the other mode, and a mode switch changes no
  setting. **Reset Settings** resets both modes, and one **Undo** restores both.
  The input files, LOSAT run settings, and the rich feature popup stay shared.
  This fixes Show Depth turned off by a switch or by the other mode's Depth file
  removal (OV-82, OV-101); a Linear Generate that failed after two Circular
  Depth series; a Circular Generate that reset the Linear Label rendering
  (OV-159); a Legend style or rename on a row the other mode does not draw that
  failed Generate (OV-80); a Legend rename and its styles lost while Generate
  hides the row (OV-120); the other mode's Legend reordered or repainted on a
  switch (OV-142); an annotation's target record picked in one mode lost after
  a switch to the other mode (OV-105); an unmanaged config override of one mode failing the
  other mode's Generate and Save (OV-106); and the Depth panels writing the
  Depth series when they draw (OV-109). A source replacement now also removes
  the strokes of features the new source does not have (OV-84) and the Legend
  renames of rows the new Result does not draw.
- Sessions: **Save Session** writes the settings and the Result of each mode.
  Session 46 keeps each mode's settings and edits in `modes.circular` and
  `modes.linear` and gains an optional `otherModeResult` field for the Result
  set of the mode that is not at the top level; its request uses the same
  resource table. Loading a Session 27–44 gives each mode its settings from the
  one saved draft; the development-only Session 45 is rejected. **Load Session**
  shows the saved mode when it has a Result, otherwise the mode that has one.
  A Result from a Session saved by an older gbdraw that waits in the other mode
  still needs one Generate before **Save Session**: the message names its mode,
  and its **Generate** button switches to that mode and generates there.
  `gbdraw circular --session` and `gbdraw linear --session` render the Result
  set of their own mode; a re-save with `--save_session` or `--session_output`
  replaces only that set, where it is, and keeps the other mode's set and the
  saved `ui.mode`.
- CLI: `gbdraw circular|linear --session <file> --session_output out.json` (and
  `--save_session`) keeps the resource IDs and file names of the Session's unchanged
  inputs. Before, the rewritten Session renamed them to positional IDs such as
  `record-1-genbank`, and table files lost the names they were uploaded with.
- Depth tracks: Depth TSV positions now follow a crop and a reverse complement,
  in Linear and Circular, on the command line, in the Python API, and in the web
  app. A crop such as `--region chr:601-800` drew the TSV rows at positions 1-200
  in place of the rows at 601-800, and a reverse-complemented record drew its
  coverage mirrored. Positions are source coordinates of the named record; a crop
  keeps only the positions inside it, and a reverse complement flips them. The
  automatic Depth maximum now comes from the drawn positions only.
- CLI: `gbdraw circular|linear --session <Session 41-44> --session_output out.json`
  (and `--save_session`) no longer fails with "Feature placement drafts require a
  circular or linear scope." for a Session whose Feature placement drafts were saved by
  the Web app. The diagram was written but the Session was not. The replay now gives
  each draft a mode as the Web app does on Load: a lane placement the mode of its side,
  a Main placement both modes. A Session whose drafts cannot be written now fails before
  any diagram is rendered.
- Legend colors (web app): a Legend color or stroke on a row that only track data
  names (an annotation set legend label, a Depth series) is retired together with
  that data, so Generate no longer fails with "The generated result could not be
  accepted" after the set, the Depth file, or the Depth track is removed. Styles
  follow the caption: replacing the data keeps the style of a caption the new data
  still names. Undo restores the data and the style.
- Legend colors (web app): Generate no longer fails with "The generated result could
  not be accepted" after a Legend color or stroke is set on a **Depth** row and
  **Show Depth** is then turned off, in Linear and Circular. The row is not drawn, and
  the stored Legend style applies again when Depth is shown.
- Depth tracks (web app): Undo of a removed **Depth TSV** file or Depth track
  now shows the Depth track with its tick text (for example `0.13x`) on the
  current Result, as Generate draws it. The track came back without its ticks
  until Generate.
- Legend colors and names (web app): Generate no longer fails with "The generated
  result could not be accepted" after a Legend color, stroke, or name is set on a row
  whose features are then all hidden (Feature visibility Off, a Feature Visibility
  rule, or a color rule that recaptions them). The row is not drawn, and the stored
  Legend edit applies again when the row returns. A Legend edit for a name that no
  feature can produce still fails Generate. The Python engine now reports, for each
  Result, the Legend rows it drew and the rows its draft removed.
- Feature Visibility rules (web app): adding, editing, moving, or deleting a
  rule in **Features → Feature Visibility** now updates the current Result at
  once, as Generate draws it, like the other visibility edits. The change used
  to appear only after a label rerender, Undo, a Result switch, or Generate.
  A rule whose value regex Generate rejects stays in the list, leaves the
  Result as it is, and shows the error Generate reports, with the same line.
- Features list and Search features (web app): a hidden feature now stays
  listed, so it can be shown again. The list holds the displayed Result's
  features of the selected types, every feature it draws, and every feature
  with its own **Feature visibility**. Each row has a **Visibility** checkbox
  that shows whether the feature is drawn and sets its **Feature visibility**
  **On** or **Off**; **Edit** and **Open** open the popup of a hidden feature,
  which says that it is hidden. A feature used to drop out of the list and
  search after Generate once it was hidden, and its popup could not be opened.
  Showing a feature that the Result does not draw renders it again live, also
  with **Auto Reflow** off. A GFF3 input loads the type of a feature whose
  **Feature visibility** a row turns **Off** as well as **On**, and a `source`
  feature with its own **Feature visibility** stays in the feature catalog.
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
- Label whitelist, Qualifier priority, and Default colors file imports (web
  app): a row with too many or too few tab-separated columns is now rejected
  with a table error that names the row and the required column count. The
  import used to skip a short row and ignore extra columns, so a file that the
  CLI rejects was accepted with rows dropped.
- Default colors, Specific colors, Qualifier priority, Label whitelist or
  blacklist, Label overrides, and Feature visibility tables: a row with more
  columns than the table has is now an error that names the file and line
  ("Malformed line ... expected N columns"). Before, a Label whitelist,
  Qualifier priority, or Default colors file whose first row had an extra
  column was read shifted: `CDS`, `product`, `two`, `words` became feature type
  `product`, qualifier `two`, keyword `words`, and no error was raised. Remove
  the extra cells, or the tab inside a value, from such files.
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
  preview, after Generate, and through Save and Load. Session version 46 stores
  them as each mode's `features.featureOverrides`; Session 44 and older Web Sessions move
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
[Retired inputs](./docs/REFERENCE/session-and-request-compatibility.md#retired-inputs).

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
- Layout (CLI and Python API): Circular track rows next to a row with a **Radius**
  follow two rules. Two adjacent rows keep the larger of their facing gaps, so an
  explicit **Outer gap** after a pinned row no longer fails with "Circular track slot
  order cannot be honored" (GX-04). A row with a Radius and an Auto width compresses
  like Auto, centred on its Radius, so typing a row's Auto radius back gives the Auto
  layout and a Radius just above it renders (GX-05, GX-17). A row with an explicit
  width keeps it.
- Layout (CLI and Python API): a custom stack that repeats the default stack, which
  the web app sends when **Use custom stack** is turned on and nothing is typed,
  draws the default figure. Before, the GC rows packed under the ticks instead of
  keeping their preset radii (MG1655). Tracked Sessions and reference SVGs do not
  change (GX-19).
- Custom Track Slots (web app): changing a row's renderer keeps only the parameters
  the new renderer accepts, so Generate no longer stops with "retains '…' from
  another renderer" after Ticks to Dinucleotide content or Depth to another renderer
  (TK-04). An invalid Circular Depth track index (`1.5`, `-1`) stays in the field with
  the message "Enter a whole number of 0 or more." and the row keeps its last valid
  index (TK-07).
- Custom Track Slots (web app): **Hide GC Skew** and **Hide GC Content** reach only the
  rows that use the drawing's **Dinucleotide** setting, so an AT skew row stays enabled
  (TK-09). Reset and Reset to preset build one Depth row per loaded Depth series, also
  while Show Depth is off (TK-10). **Duplicate** is disabled while a Session operation
  runs (GX-01).
- Custom Track Slots (web app): a Depth series without a Depth TSV fails Generate with
  "Attach a Depth TSV to this series, or remove the series." (TK-06). The help tips
  for a Depth series name say that typing renames the series everywhere and what a
  blank gives (TK-08). **Width** and **Radius** reject `0x10`, `0b11`, and `0o7`, and
  a bad value gives one alert that names the field (TK-12). A stack without a Features
  row generates when the record has no underlay features, as in the CLI (TK-13).
- Custom Track Slots (web app): before the first Generate the Auto notes read
  "≈ 37 px (estimate)", and the ticks row is estimated instead of "0 px" (TK-15).
  After Generate, the **Radius** note of a ticks row shows the tick anchor, so typing
  it back keeps the ticks where they were (GX-18).
- Generate (web app): Generate with no input says so at once instead of after the
  diagram Worker starts (5.7 s warm, 26-42 s cold) (UI-07). An input error band
  (`NO_RECORDS`, `FASTA_REQUIRED`, `INPUT_REQUIRED`) ends when the replaced or
  completed input is read; other errors stay (UJ-08). Canceling a Generate when
  nothing changed since the shown Result says "It matches the current settings."
  (UJ-10).
- Input files (web app): replacing or removing a Circular GenBank file resets the
  single-record selector, crop, reverse complement, Record label, and Subtitle in one
  History step, in Circular and Linear, so Generate no longer fails with
  `RECORD_SELECTION`. A one-record file replacing a one-record file keeps the crop and
  titles (CI-03).
- Run info (web app): the Source recipe of a Circular multi-record canvas keeps the
  order set with Up and Down, so the CLI draws the records in the order the web app
  does (CI-02).
- Labels (web app): a label drawn only because of Label visibility **On** or **Off**
  returns to the Default state after Undo, the popup's Default, Reset all label text,
  and Label TSV import, as Generate draws it (UJ-01, GX-21).
- Region Annotations (web app): the panel notice names an annotation or set id the
  editor changed (for example `region_1` to `region_1_2`). Importing a table without
  rows asks before it removes every set, and Import errors appear in the panel notice
  instead of a browser alert (FL-14).
- Feature captions (web app): a feature with no label, product, gene, locus tag or note
  is named `tRNA at 577..647` (1-based, split locations joined) in the popup title,
  the label text default, the default Legend name, and the Interactive SVG popup.
  It read `576..647`. Stored Legend names keep their text (GX-10).
- Load Session (web app): **Load Session** asks "Replace the current work?" before it
  replaces work changed since the last Save or Load (UJ-09). On gbdraw.app the empty
  state offers **Load an example**, which lists every Gallery example and loads the
  chosen one's Session, fitted to the Preview; the local `gbdraw gui` ships no Gallery and shows no such
  button (UJ-06).
- Preview (web app): the Preview feature search belongs to its mode and a Session load
  clears it (observation). After a Generate, the Linear Auto notice says "Auto
  shows/hides these fields" when the shown Result already does (observation).
- Popups (web app): Escape closes the top layer only, a popup before the Editor
  drawer, and focus returns to the drawer's **Edit** button. A popup keeps clear of
  the footer when opened low or dragged down, and releasing a drag past its limit no
  longer closes it (GX-22).
- Legend Name Scope (web app): **Cancel** closes the dialog at once and records no Undo
  step, also right after a Session load (GX-20).
- Comparison rings (CLI and Python API): a ring without `--conservation_labels` is labelled with
  its comparison file name without the last extension, as in the web app (`NC_002333.2.fna`
  draws `NC_002333.2`). This covers LOSAT rings and precomputed `--conservation_blast` rings.
  Sessions now store every ring label, and a saved Session keeps the labels it was drawn with.
  GenBank and DDBJ rings still take the DEFINITION, then the organism.
- Scale interval (CLI, Python API, web app): the CLI and the Python option objects reject a scale
  interval of 0 or less with an error that names `--scale_interval` or `objects.scale.interval`,
  and the web app's Scale Interval field starts at 1. Before, the value was accepted and drew the
  automatic interval. A Session or web request that holds such a value still draws the automatic
  interval.
- Tick labels (CLI, Python API, web app): a manual interval that is not a whole unit prints the
  decimals it needs. A 500 bp Circular interval read "0 kbp, 1 kbp, 1 kbp, 2 kbp", and a 250 bp
  Linear interval wrote 2250 as "2.2 kbp". Automatic intervals are unchanged (OV-201).
- Sessions (CLI and Python API): the LOSAT runtime record of a Session stores only the
  executable's name for a runtime outside the package (explicit, conda, managed, or PATH), not its
  absolute path. The command line in `cliInvocation` is still stored as typed.
- Sessions: replaying a Session 27–30 with an unlabelled precomputed ring labels it with the
  original file name instead of the temporary `arg<n>-<name>` (OV-202), and the web app loads an
  older CLI Session without ring labels with the full file name, as the CLI draws it (OV-203).

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
