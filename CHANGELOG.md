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

## [0.14.0](./docs/RELEASE_NOTES_0.14.0.md)

- Added circular-record display-start rotation and manual feature lane placement.
- Improved Run Info / Exact replay, Save Session resource preservation, preview
  navigation, and Generate processing-stage feedback.
- Includes the beta's package-root Python API and current session/request
  compatibility, plus isolated wheel/sdist installation and local GUI packaging fixes.
- See the [full release notes](./docs/RELEASE_NOTES_0.14.0.md) for migration,
  compatibility, and installation availability. Publication dates are recorded
  in [GitHub Releases](https://github.com/satoshikawato/gbdraw/releases).

- Fonts (CLI, Python, web app): diagram text can use any font, including a
  font installed on the computer or loaded from a file, and gbdraw measures
  labels, the legend, and titles with that font instead of Liberation Sans
  (#1046, #1058).
  - CLI: `--font_family` takes a family name or a CSS font-family list.
    `--font_file PATH` (repeatable; `.ttf`, `.otf`, `.ttc`, `.otc`) adds a font
    file, also with `--session` and `gbdraw render --session`. In Python, set
    `objects.text.font_family` through `config_overrides`. When the first
    family of the list is found nowhere, the text is drawn in the next family
    found, else the nearest Liberation font, with a warning.
  - SVG and Interactive SVG embed the characters they use of a font that is
    not bundled with gbdraw, as `@font-face` data in
    `<defs id="gbdraw-embedded-fonts">`. A font whose license forbids
    embedding, or allows only bitmaps, is named in a warning and not embedded.
    `--no_font_embedding`, `embed_fonts=False`, and the web app's **Embed font
    in exported SVG** keep only the name in SVG and Interactive SVG.
  - PNG draws the chosen font whatever the embedding setting and the font's
    license. PDF, EPS, and PS never contain a font whose license forbids
    embedding, allows only bitmaps, or forbids subsetting; that text is drawn
    in the nearest Liberation font, and a warning (CLI, Python) or notice (web
    app) names both. From the command line and Python, PNG, PDF, EPS, and PS
    draw a `--font_file` font that is not installed where cairo uses
    fontconfig, as on Linux.
  - Web app: a new **Font** section offers Liberation Sans, Serif, and Mono,
    **Load font file…**, and **Choose installed font…** (Chrome and Edge).
    SVG and Interactive SVG exports embed the font, and PNG export draws it.
    PDF export embeds the characters it draws, converting OpenType CFF
    outlines to TrueType. In PDF export, a font that cannot be embedded is
    drawn in the nearest Liberation font, and a notice names both.
  - Sessions keep the font name, not the file; saving warns about it. Loading a
    Session whose font is not available asks to load the font file, use an
    installed font, or load with the nearest font, which keeps the name for
    later.
  - Fixes: OpenType CFF fonts are measured by their glyph outlines, and a font
    without a Windows Unicode cmap no longer stops measurement (OV-680..683;
    PD-OI-092, PD-OI-093, PD-OI-098).
- Circular labels (CLI, Python, web app): after the automatic repair of small
  label defects, the labels around the repaired defects make room (up to 12
  neighbours on each side; moves within 24 px of text and 12 px of elbow), so
  their leader lines are shorter and straighter. The step never adds a defect
  or a leader line through another label's text. A repaired map takes about 1
  to 8 s longer to draw and its canvas can grow by up to about 24 px. Maps that
  fit without a repair are byte-identical.
- The label fit warning names the search bound the caller passed (the web app
  passes 0 s) instead of always saying 5 s.
- Legend names (CLI, Python, web app): legend rows that share a name follow
  one rule, and no legend name carries a color ID. Rows with the same name and
  the same swatch (fill and line) are one row. When two or more rows with the
  same name and different swatches are drawn, every one is numbered in legend
  order (`Packaging (1)`, `Packaging (2)`); a name drawn once has no number.
  The rule replaces the `name [#hex]` captions of one caption in several
  colors and of a rule named like a generated row (now `tRNA (1)` and
  `tRNA (2)`), and it covers track and region-annotation rows with one
  `legend_label`.
  - Changed drawing: two depth tracks labelled `Coverage`, which 0.13.0 drew
    as `Coverage` and `Coverage (2)`, now draw `Coverage (1)` and
    `Coverage (2)` when their colors differ and one `Coverage` row when they
    are equal. Two GC content tracks always have the palette's GC content
    color, so 0.13.0's `GC content` and `GC content (2)` are now one row. A 0.13.0
    Session still loads, and its legend row edits stay on the rows of the same
    tracks; edits saved for a second row that now shares the first row are not
    used.
  - gbdraw no longer rewrites a rule's caption. The color table, the rule
    list, the popup, the downloaded color-table TSV, and the Session keep it
    as written, and the web app notice "Updated N specific-color caption(s) to
    distinguish their colors." is gone. The numbers follow the legend order:
    they change when a row is added, removed, or moved, and a name left with
    one row loses its number.
  - Compatibility: when a color table file or a saved Session is read, a
    caption that ends in its own rule's color ID, such as
    `Packaging [#ffe700]` on a rule colored `#ffe700`, loses that ID, also
    when you typed that name yourself.
  - Web app: the Legend tab, the popups, and the dialogs show the drawn name.
    The feature popup's **Legend name** field no longer offers suggestions;
    the list under it is the only place that offers an existing group
    (OV-583, OV-584, OV-585, OV-587; PD-OI-042, PD-OI-095).
  - Web app: **Keep both** in the **Legend Name Conflict** dialog keeps the
    name you typed, and every row of that name is numbered in legend order.
    Before, the renamed row alone got a numbered name and the other row kept
    its name. The button shows the name the renamed row is drawn with
    (**Keep both: name it "CDS (2)"**), and a sentence names the other row's
    new name when it changes (`The existing "CDS" becomes "CDS (1)".`).
    Typing a numbered name the legend shows, such as `CDS (2)`, opens the
    dialog for that row, and **Keep both** keeps the text as written. Rows of
    one name and one swatch are one row, also in the live preview, so
    **Keep both** with an equal swatch gives one row; a color or line edit of
    that row changes every row it stands for. Features renamed off a
    row get their new row right after the row of the name (PD-OI-061).
  - Sessions: the Web app now saves Session 47. A Legend row's color and line
    edits are kept under the row's key instead of the name it shows, so a
    renamed row keeps its edits and two rows of one name keep their own
    colors through Generate, **Save Session**, and **Load Session**. A
    Session 46 is read by moving each edit saved under a renamed row's shown
    name to that row; the CLI and Python re-save and upgrade it the same way.
  - Fixed (web app): a color rule that would draw a Legend item that already
    exists in another color was refused with "The operation failed without
    recognized diagnostic information." and a **Retry** that failed again. It
    now says "A color rule would draw a Legend item that already exists in
    another color. Change the rule's Legend name or color, or rename that
    Legend item." (OV-644).
  - Fixed (web app): the color scope dialog hid or offered the **Label** and
    **Source label** scopes by comparing the label text with the Legend row's
    key. It now offers a label scope only when it selects other features than
    the **Legend item** scope (OV-648).
- Region annotations (CLI, Python, Session): an overlay annotation slot can
  span several tracks. With `to_slot=<slot id>` next to `anchor_slot`, the
  overlay runs from one slot's band to the other's in Circular and Linear, for
  example one tint per region from the GC track to the feature track
  (`regions:annotations@set_id=regions,side=overlay,anchor_slot=features,to_slot=gc_content,layer=underlay,z=-1`).
  `to_slot` follows the anchor rules, and the span adds no space between the
  slots; in Linear, a span that shows labels keeps room for them at its edge.
  The web app draws and saves a loaded `to_slot`; its slot editor has no
  control for it yet. The CLI reference now documents the overlay slot
  parameters (`anchor_slot`, `layer`, `cover_anchor`, `to_slot`). Fixed: a
  Python slot mapping that spelled the key with capitals (`Anchor_Slot`) was
  accepted but drawn on the axis instead of on its anchor slot (OV-475).
- Linear tints that cover the feature track (a `highlight` annotation row
  without custom slots, a feature drawn as an underlay such as
  `repeat_region`, and an overlay slot with `cover_anchor=true`) now cover
  the feature lanes with `--track_layout above` and `below`. They were drawn
  9.5 to 18.5 px nearer the axis than the lanes, and their own height pushed
  the lanes away from the axis and left empty space. The tint now adds no
  space, its labels keep room beyond the lanes as for `to_slot`, and a region
  without a label keeps none (OV-474, OV-476, OV-842). On a record whose only
  features are underlays, a `to_slot` span from the feature track now covers
  their nominal lanes too. Fixed: a feature-bound ribbon of a feature moved to
  another lane stopped at the outermost lane of a `to_slot` span instead of
  at its own lane.
- Linear span labels (a `highlight` annotation row without custom slots, an
  overlay slot with `to_slot`, or one with `cover_anchor=true` on the feature
  track) no longer lie on the next slot or on the axis line. The room for the
  labels (9 px with the default label style) is now kept between slots too,
  not only at the edge of the record: the next GC, skew, or feature-lane slot,
  or the axis line, moves out by up to that much on the side with label text.
  With the default `--track_layout below`, a labelled `highlight` row no
  longer puts its label on the axis line; the lanes move away from the axis
  only by the part of the room the gap did not already give (2 px for the
  tobacco chloroplast and AP027078). Figures without a labelled span are
  unchanged (OV-843).
- Linear annotation labels (a labelled row, `to_slot` span, or annotation
  slot, with any mark) keep `label_offset` (4 px by default) between their
  painted letters and the painted edge of their mark, its stroke and a
  bracket's end caps included, and the same gap beyond the letters to the
  axis line, the next slot, or the record edge; a `to_slot` or
  `cover_anchor` span's label stays outside every slot the span covers,
  whatever its mark and its own `h` (OV-1156). `highlight`, `band`, and
  `bracket` labels no longer touch or overlap their mark; `line` labels,
  which had 4.7 px, move 0.7 px closer. The label room grows to match;
  records grow only where nothing else already gave that room, by 0 to 7 px
  per record in the tested figures (OV-844). A label that crossed a feature
  label's leader line or overlapped a feature label now moves along its
  region, its centre over one of the region's marks (never into the gap of a
  region across the origin) and its letters on the bar (OV-1157), to the
  nearest clear place; labels without such
  a conflict stay where they were, for example the tobacco chloroplast `IRb`
  and `IRa` labels with `--track_layout above` (OV-845). Where no place is
  clear, the label goes to the place with the fewest conflicts and is drawn
  there without a fit warning (OV-1155). The label room now also follows a
  slot's `style_override` font size and offset (OV-848).
- Web app feature popup: the header shows the gene name, the product on up to
  two lines, and the record and location with thousands separators. **Edit**
  groups **Appearance** (fill, stroke, the label text with its visibility and
  one **Apply Label**, feature visibility, legend name) before **Similarity
  group** and **Layout**; the always-shown hints, the label key, and the source
  label line are gone, and a missing label is explained in one line.
  **Details** and **Qualifiers** list each field beside its value with one copy
  button per tab. The record rotation buttons are now **Apply when I Generate**
  and **Apply now**, and a record that cannot rotate shows its reason under a
  disabled toggle.
- Web app popups: feature and match popups open inside the window beside the
  clicked item instead of clipped at the bottom right (OV-371), and a popup no
  longer covers the header buttons when opened or dragged (OV-384). A popup
  opened with **Edit** in the Editor's Features list opens beside the Editor
  when the window has room there (OV-503).
- Web app on a phone: the docked Editor uses most of the space between the
  canvas and the Generate bar (about 400 px instead of 200 px at 390x844).
- Web app: a click on a Legend row or a record definition in the preview opens
  a popup while Layout edit is off. The Legend item popup shows the feature
  count and the rule behind the color, and edits the color and legend name or
  removes the row. The record popup lists the accession, length, reverse
  complement, display start, and region, in Linear also the organism line and the
  alignment, and edits the definition text and Reverse complement, which apply
  on Generate.
- Web app: a click on the Circular scale ticks, the Linear scale, a GC content
  or GC skew (dinucleotide) track, or the plot title in the preview shows that
  element's settings section in a popup beside the figure while Layout edit is
  off: **Axis & Scale**, **Dinucleotide content/skew**, or **Titles & Record
  Labels**. The pointer on one of these elements shows the hint **Click to open
  its settings**. The popup holds the section itself, so the settings panel does
  not list that section until the popup closes; an edit there is the same
  setting as in the panel and applies at the same time. This popup and a
  feature, match, Legend row, or record definition popup replace each other,
  so the preview shows one popup at a time. Every popup of the preview now moves back inside the window when the
  window is resized, instead of staying partly outside it (OV-626).
- Web app: with **Layout edit** on, a Legend row shows the grab cursor, and the
  grabbing cursor while you drag it (it showed the pointing hand). With
  **Layout edit** off a row keeps the pointing hand and opens its popup. In a
  Circular Multi-Record Canvas with **Layout edit** on, a record definition
  shows the grab cursor where it showed the pointing hand (OV-604). On a
  Linear figure with **Layout edit** on, a record definition no longer shows a
  grab cursor where nothing can be dragged (OV-625).
- Web app: the Editor's Features list uses one compact row per feature that
  expands to its label text field; the label and feature-edit TSV actions move
  into one **Import / export** menu. The Legend tab reorders rows by a drag
  grip, also with the Up and Down arrow keys, and its help moves into a help tip.
- Interactive SVG: the similarity-group section of a match popup shows the
  group name and gene count, the members, and the identifiers under
  **Technical details**, as in the web app.
- Web app: confirmations use the app's own dialog with named buttons instead of
  the browser's OK/Cancel box: **Reset Settings** lists what it resets, and a
  risky color rule pattern, a repeated or large Session download, a large
  reproducibility or raw LOSAT TSV download, and hiding a GC track of a custom
  stack each name their action. Reset Settings also says that it turns the
  Linear comparisons off. Load Session, Save Session, and table imports and exports report
  in a dismissible notice above the Result instead of a browser alert; Session
  warnings show as a list.
- Web app: the feature-edit choice dialogs (color and stroke scope, label
  visibility, label text scope, Legend name, palette colors, Feature placement
  reset) show one primary choice, group the shared scopes, and name each option
  by its outcome; a Legend name conflict shows the name **Keep both** gives.
- Web app: feature search fits one row, shows **Qualifier key** only for
  qualifier search, labels the regex toggle **Regular expression** (also in a
  downloaded Interactive SVG, whose invalid-pattern message now names it), and
  shows a match count only after a search. **Reset Layout** is unavailable until
  something moved, and **Canvas padding** fields show their unit (px).
- Web app fixes: the Editor's Legend tab showed the default skew colors on the
  rows of a skew track with its own colors, and a live palette change repainted
  those rows in the preview legend; both now keep the track's colors (OV-484).
  The Legend rows of a second skew track with the same label now follow a live
  palette change, and a palette change no longer applies a track slot color
  that still waits for Generate (OV-488).
  A Session whose request names a Circular track slot's renderer by the alias
  `content`, `gc_content`, `conservation` or `circular_comparison` (a Python
  API request) now loads it as that track: a spec string was read as a skew
  track (OV-489), and a slot object, also with `skew` or `gc_skew`, stopped
  the Load (OV-491). A Linear content or skew slot without a dinucleotide (a
  CLI or Python API Session) now loads, and shows in the Source recipe, with
  the Session's dinucleotide, not GC (OV-490). In a drawing with a dinucleotide
  other than GC, the empty dinucleotide field of a Custom Track Slots row and
  the Circular preset summary said GC; both now name the drawing's dinucleotide
  (OV-497).
- Raised the short/long size boundary of the Auto sizes: Linear diagrams now
  use 100 kb and Circular diagrams 75 kb; both used 50 kb. A record of 50 kb
  has too few features for the long styles, which are sized for whole genomes.
  Circular uses 75 kb, not 100 kb, because the short label size overprints
  outer labels on denser records near 100 kb. Only Auto sizes change; set a
  size explicitly to keep the previous look. What changes:
  - A Linear diagram whose longest record is 50,000 to 99,999 bp now gets the
    short styles (larger labels, features, strokes, and legend).
  - A Linear diagram with records of both size classes now sizes its labels
    like its features, by the longest record in the label scope
    (`--show_labels`): with `--show_labels all`, a record under 100 kb beside
    one of 100 kb or longer gets the long label size (5 pt). Before, one
    record under 50 kb in the scope made every label 24 pt on 20 px features
    (OV-482). With the scope limited to a short record, as with
    `--show_labels first`, the labels keep the short size; a label filter
    such as a whitelist does not change the scope. A Session saved earlier
    is drawn by this rule too and keeps its stored boundary.
  - A Circular record of 50,000 to 74,999 bp drawn alone or beside shorter
    records now gets the short styles.
  - A Circular Multi-Record Canvas with records of both size classes now
    takes every Auto size of every record, and of the legend, from the long
    class (OV-492). A record of 50,000 to 74,999 bp beside one of 75 kb or
    longer is drawn as before. A record under 50 kb beside a long one
    changes: its label leader lines are 0.5 px wide like the long ones
    (OV-483; they were 1.5 px beside the 8 pt long label font), and its
    labels sit at the long class's distances. When that record is the first
    one, the legend changes too: 16 px text and 20 px swatches, where the
    first record's class gave 20 px and 24 px. Before, such a canvas copied a
    fixed list of long sizes to its short records and left their other sizes
    short. A Session saved earlier is drawn by this rule too and keeps its
    stored boundary.
  - A Session that stores its whole configuration, such as one written by
    `--session_output`, keeps the boundary it was saved with. In the Web app,
    the Auto values beside the size fields and the Circular track preview
    follow that boundary, as Generate does.
  - The aminoglycoside BGC Gallery figure is unchanged. Its tutorial now
    types Head Length Ratio 0.64 in step 10, **Finish the figure**, and its
    intermediate previews show the short Auto sizes.
  - In the Web app, the Auto values beside the size fields and the Circular
    track preview now read the records Generate draws. The Linear **Label Font
    Size** value follows the records in the label scope (**Show Labels**), so
    with **First Record Only** on a short first record it reads 24 px, not
    5 px (OV-493). The Circular values follow the records that are drawn: a
    Multi-Record Canvas as one diagram, separate diagrams one record each,
    and neither a record turned off nor one that is not selected (OV-494).
    They read **Auto (by genome size)** when separate diagrams take different
    sizes. A record cropped to a region is sized by the length of the region,
    as Generate sizes it: a 306,008 bp record cropped to 1..30,000 reads
    14 px for **Label Font Size**, not 8 px (OV-512). The sliding window and
    **Scale interval** follow each Circular record's own length, as Generate
    does; where the drawn records take different values they read **Auto (by
    genome size)** and `by record length (auto)` (OV-512). The track preview
    estimates the diagram that draws the record of the shown Result, or the
    first diagram while no Result is shown.
- Web app: the Feature Visibility rule action **Exclude from matching** is now
  **Leave out of protein search**, also in the selected-feature bar's
  visibility list, where it read **Suppress**. The color rule field **Legend
  caption** and the selected-feature bar's name field are now **Legend name**,
  the words the feature popup uses, and the selected-feature bar's color
  button and the Reset Fill Color dialog say legend name too. Saved Sessions
  and rule files are unchanged. A Feature Visibility rule shows **Action** on
  its own line under **Value regex**, so the closed list shows the whole
  action.
- Web app: a Circular record counts as cropped only where Generate crops it:
  one record offered alone, with a region Generate can draw (OV-1093). Its
  **Display start** was unavailable, as for a cropped record, when the region
  fields were emptied, or when separate diagrams or a Multi-Record Canvas draw
  the record whole. With a region left in the fields while separate diagrams
  or a Multi-Record Canvas are drawn, the record's reverse-complement override
  and **Display start** now count in the record list, as they already did in
  the drawing.
- Fixed (web app): on a phone, or with a short Preview, the selected-feature
  bar's color, visibility and stroke rows ran under the preview controls and
  could not be reached. The bar now stays inside the Preview and scrolls
  inside itself, and it no longer runs under the Editor toggle at 1024 px. Its
  three edits wrap to new rows when the Preview is narrow (OV-1098).
- Fixed (web app): loading a Session that stores a full resolved
  configuration, such as one written by the CLI `--session_output`, no longer
  turns the packaged font list into an explicit **Font**. **Font** stays at the
  default, Run info lists no `--font_family`, and **Save Session** writes no
  font override; the drawing is unchanged. The Gallery Sessions are saved again
  with the current writer (OV-1014).
- Web Comparison card (Linear): **Settings** lists the mode, the result filters
  and the appearance first, then **Search limits** and **Performance** (formerly
  **Runtime and reproducibility**). **Selected pairs** shows each pair as one row
  with its source, file, problem and **Save Raw LOSAT TSV**; the separate **Raw
  LOSAT results** list in **Advanced comparison and layout** is gone. Pair-plan
  problems appear once, in a **Fix before Generate** list with an action each.
  The Collinear **Evidence scope** setting is now **Search between**. An open pair
  **File name** field takes its own row under **Save Raw LOSAT TSV** instead of
  the narrow column beside it.
- Web errors appear once, at their cause: a failed Generate no longer repeats the
  "last successful result" notice, a Linear comparison error offers **Show in
  Comparison**, an invalid colour-rule pattern shows its cause under the field,
  and a failed Session load shows **Session not loaded** with the reason.
- **BLAST TSV missing** lists one row per pair with its own **Choose file…**,
  then **Generate** or **Skip these comparisons and Generate**.
- The first Generate shows each runtime preparation step ("step 2 of 5, ...").
- Fixed: a Session load error raised in Circular stayed on screen after switching
  to Linear (OV-372).
- Fixed (web app): a Session or configuration file whose legend position is not
  one of its diagram mode's positions (a Circular corner such as Upper Left in
  Linear, or any other text) loaded, and **Generate** then failed with "Input
  validation failed". Such a position is now read as unset, so the legend takes
  the position an unset one gets; capitals and outer spaces are still accepted
  (OV-652). The side a loaded Result's Legend was laid out with is also read
  through its mode's positions, and a side outside them reads as the mode's
  default: a Linear Result saved by 0.13.0 with a corner position keeps its
  figure, and its first Legend edit places the Legend below the diagram; such
  a side with capitals or spaces no longer stops the load (OV-671, OV-672).
- Fixed (web app): a Session that 0.13.0 saved after **Generate** and a switch
  to the other diagram mode opened its Result in the shown mode, listed its
  feature information as not restored, and moved the Legend and the title at
  its first edit. When the Session records the mode of its last Generate, it
  now opens in that mode, with the Legend side and title side the Result was
  drawn with; the shown mode keeps its settings, without a Result. A Session
  0.13.0 saved without that mode (for example after a Linear LOSAT comparison
  run) still opens in the shown mode (OV-677, OV-969).
- Changed (CLI `--session`, Python API): the canonical request reader refuses
  values the request format does not define: a comparison table path at any
  request schema; `null` or another kind of value where a field holds text, a
  flag or a name; and, at request schema 9, anything outside
  `describe_canonical_request()`, other spellings of the comparison mode
  settings, and other spellings of `conservationReference`,
  `pairwiseMatchStyle` and `generatedProteinComparison.mode`. A fractional
  annotation `lane`, `start` or `end` still reads as its integer part. The
  encoder writes a track slot's `renderer` and `side` by their names, so a
  slot built with an alias such as `gc_content` reads back as
  `dinucleotide_content`. A Session whose stored request holds a value that
  earlier versions read and this reader refuses still loads in the web app,
  and its Generate shows "An input value is invalid". A record display
  selector with non-ASCII digits is refused, as the web app already did
  (OV-324, OV-528, OV-596, OV-611, OV-631, OV-632, OV-651, OV-653).
- Changed (web app): a Session or configuration file that Load cannot fully
  restore now opens one **Load this file?** dialog before anything changes. It
  lists the settings, Legend positions of the layout preferences or a saved
  Result, Legend entries and edits it re-keys, table rows, Depth rows that
  name a Depth series their drawing does not have, feature edits and feature
  information Load drops, replaces with a default, or reads differently (up
  to 12 lines): a setting with its place in the file and the value used, a
  feature edit or feature information as a count or a statement. A setting
  that a Session of version 44 or older stores once for both diagram modes is
  one row. It says that Save Session afterwards writes the current format. Values that had no
  effect on the figure (a pending palette Instant Preview drops, edits of
  similarity groups the Result does not have) are dropped without a row. A
  Depth series without a Depth file loads as saved, and so does a Depth row
  the app keeps when a Depth file is removed (a row with its own size or
  position) (OV-755). An unknown text value of `track_type`, `scale_style`,
  `labels_mode` or `show_labels_linear` in a Session or configuration file is
  still kept as stored, because the request format lists no values for them. No file that
  loaded before is refused: a feature shape or Linear Accession or Length /
  Coordinates visibility that is not text keeps the drawing's value, and an
  unknown shape for a feature type the defaults lack is dropped, each with a
  row (a Session of version 44 or later with such a visibility, refused
  before, now loads the same way). In a configuration
  file, a color rule or label whitelist entry that is not an object, which
  made the Load fail, is now dropped and listed, and a label filter mode or
  feature overlap tolerance the app does not accept, which made Save Session
  fail after the Load, is now listed and set to its default. A configuration
  file loads over the current drawing: a settings group (`form`, `adv`,
  `losat` or `circularConservation`) that is not an object, a value of
  another kind inside one (except the values kept as stored, above), and
  `rules` or `whitelist` that is not a list (such as `"form": 5`) are not
  used, and the drawing keeps its value, which
  the row names or says the drawing keeps. **Load** loads the file;
  **Cancel** leaves the settings, Results, input files, History and screen as
  they were. Dropped feature edits, narrowed Feature visibility edits, `hash=`
  values that could not be renamed, joined or dropped table rows and Legend
  entries with a color the app does not accept, which notices reported after
  the load, are listed in this dialog instead. The notice of `hash=`
  annotation targets moved to the feature they name, and the notice of a
  loaded configuration file, still follow the load. A file that loads without
  a loss, such as every Gallery Session, opens no dialog. A Session older than
  version 27 (gbdraw 0.11.0 or older, or gbdraw.app before 0.12.0) is still
  refused, and the error now says so (OV-710, OV-712, OV-719).
- Changed (web app): Load repairs the saved diagram's settings it cannot use
  instead of refusing the file, and lists each repair in the **Load this
  file?** dialog: an unknown setting, an invalid value, an unreadable list
  entry, a part of the settings, or a reference to a file or record the
  Session does not hold. A Session whose saved settings cannot be read at all
  loads as settings only (the saved Results, inputs, feature list and caches
  are left out) when it is a Session of version 42, 44, 46 or 47 saved by the
  web app; a Session of version 42, 44 or 47 written by the CLI or the Python
  API, and a Session of version 31 to 41, is refused, and so is a Session of
  version 31 to 41 whose request names a record file it does not hold. A
  `null` request reads as a missing one in a file that holds saved Results or
  a feature list. Eight search-scope
  spellings only the web app accepted are no longer read: a generated
  comparison that holds one is dropped with a row, and a Session's comparison
  settings that hold one read as the default, with a row. A file that holds
  `top1` or `legacy` now loads. See [Saved diagram settings that cannot be read](docs/REFERENCE/session-and-request-compatibility.md#saved-diagram-settings-that-cannot-be-read).
- Fixed (web app): since late August 2026, gbdraw.app had refused every Session
  31–32 whose request has schema 1, and some Sessions 32–33 whose request has
  schema 2, among them most older Gallery Sessions with schema 2 (OV-790). They
  load again. Load reads the request in the current format, with no
  **Load this file?** dialog when nothing is lost. Their saved comparisons,
  protein comparison settings and similarity alignment stay editable (OV-796),
  and a live edit works before the first **Generate Diagram**. **Save Session**
  writes request schema 9; right after Load of a Session with a saved Result
  it first asks for one Generate, as for other Sessions saved before Session 40.
- Changed (web app): older or hand-edited Sessions and configuration files
  that Load used to refuse for their content now load, and the **Load this
  file?** dialog lists what Load drops or replaces: an unknown or invalid
  setting, track slot, Legend or stroke edit, an unreadable Result, feature
  list, stored input file or saved diagram of the other mode, and a key that
  Session version never held. A comparison sequence file or Depth file that
  cannot be read keeps its place empty, and a Circular comparison series whose
  BLAST table cannot be read is left out with its sequence files, with a row
  that names it; the other series keep their own files and names. A Linear
  comparison whose BLAST TSV cannot be read loads without that file, as after
  removing it from its pair card. Run notices, recomputable caches, empty
  values and values with no effect on the figure give no row. Load still
  refuses a file that is not a gbdraw Session or configuration file, an
  unsupported or too old version, a file too large or with an unsafe name, and
  a file from which it restores nothing. **Save Session** refuses what it
  refused before.
- Fixed (web app): the shape set for a feature type added to the feature list
  was dropped when its Session was loaded again; Load now keeps it.
- Fixed (web app): a cleared **Feature overlap tolerance (bp)** made Save
  Session fail and Generate send an invalid value; it now draws and saves as
  its default, 0. A blank number in a configuration file now reads as the
  field's default instead of keeping the current drawing's value, and a blank
  **Column Gap Ratio**, **Row Gap Ratio** or LOSAT **Merge conflicts** no
  longer reads as 0. A cleared **Record gap (px)** under **Record Layout**
  (Linear) now draws and saves as its default, 24 px (it drew 0). Save
  Session before the first Generate no longer fails on an enabled Depth row
  whose Depth file was removed: the saved render request leaves the row out,
  and the Session's settings keep it. A Session 44 whose font was
  a preserved override and whose `adv` is not an object now lists the dropped
  `adv` in the **Load this file?** dialog.
- Fixed: the records of one multi-record File were named differently in
  **Record Layout**, **Selected pairs** and the **BLAST TSV missing** dialog
  (one by its definition, the next by its accession); each now shows the File's
  definition plus its accession, and the accession stays visible in narrow rows
  (OV-436). The Align panel keeps its own form, which always adds the accession.
- Web app sidebar: the input, Colors, Region Annotations, Titles, Layout and
  Custom Track Slots cards are more compact. Each setting that applies on
  Generate keeps a short **Applies on Generate** marker. The Colors card uses
  plain names (**Default colors file**, **Specific rules file**), shows the
  colors this figure draws and the colors you changed first and the rest under
  **More colors**, and opens each specific rule from a one-line summary.
  Annotation rows and track slots keep their less used options in a disclosure.
  Auto window and step values show the value for the current genome once its
  records are read. Text with `<i>` markup shows an italic preview. Removing a
  Linear File asks **Empty this slot** or **Remove file**.
- Web app settings panel: a row of icon buttons above the settings lists the
  sections of the current mode. A button opens its section, brings it to the
  top, and moves the keyboard focus to its heading; the section title shows on
  hover and on keyboard focus, and the section at the top is marked. Linear
  mode has entries for **Comparison** and **Advanced comparison and layout**.
  Dragging the panel's resize handle no longer selects page text (OV-534).
- Web app settings panel: **Collapse settings to icons**, at the right end of the
  section row, shows the sections as a narrow column of icons with short names
  and one section beside it; **Close section** leaves the column alone, so the
  diagram gets the room, and **Expand settings** returns to the full panel. The
  browser remembers the choice for the site address; it is not part of a Session.
  Phone-width windows always show the full panel. **Edit pattern** on a color
  rule error now opens the **Colors** section when it is closed (OV-538).
- Web app Linear File cards: a card folds to its header line (File N, the file
  name, the move buttons); click the line, or press Enter or Space on it, to open
  the card. After **Load Session** or **Load example**, the cards that hold their
  files start folded; a card you add, or one with a file still missing, stays
  open. A folded line still shows a file that could not be read, a missing file
  of a GFF3 + FASTA pair, and Gencode required. With five files, the file list is
  218 px high instead of 1,778 px. The fold is not saved in a Session, and Undo
  does not change it. A file name wider than the upload box now ends with an
  ellipsis and keeps its icon, instead of being cut at both ends (OV-600).
- Web app (Linear): a name typed on an empty File card is kept when its file is
  chosen. For a file with several records it becomes the File-level name
  (OV-486).
- Web app (Linear): a File with one record shows that record's **Organism /
  strain** and **Subtitle / title** on the File card instead of **Names for all
  records in this file**. A File with two or more records, or a one-record File
  that holds a File-level name, is unchanged. A name typed on the card of a
  one-record File is now the record's own: **Reset Settings** clears it and
  **Undo** restores it, and another one-record file chosen for the card keeps
  it.
- Web app (Linear): in **Record options**, the **Using file default** badge and
  the label of its field each stay on one line (OV-485).
- Web app: before the records are inspected, Region Annotations rows name their
  record instead of calling it unavailable (OV-373), and Linear record selects
  say **Not inspected** instead of **Loading records...** (OV-379). After
  Inspect, **Display start** shows the current start instead of an empty field
  (OV-374).
- SVG feature IDs: a feature of a cropped or reverse-complemented record
  (`--region`, `--reverse_complement`, **Reverse complement**, or a direction
  chosen in **Review alignment options…**) is now drawn with its source-record
  hash, as on an untransformed record. `data-gbdraw-feature-id`, the element IDs
  built from it, and `data-gbdraw-stable-feature-id` change for such records and
  now match the feature catalog in plain, interactive, and web-app SVGs, so
  `gbdraw.api.enrich_svg(..., feature_catalog=...)` accepts these diagrams and
  the web app's interactive export writes the same stable ID as the CLI
  (OV-401, OV-412). A `hash=` selector (color, label, and visibility rules,
  annotation targets and tables, Feature Edits) now names that source-record
  hash on every record; before, rules matched the hash after the crop or
  reverse complement. Loading or upgrading a Session saved before 0.14.0 renames
  the **This feature only** color rules it saved for such records to the
  source-record hash. Replaying a Session 31–44 with `--session` renames the
  same rules in its color table and the label text edits in its label override
  table. The Web app and the CLI report each value they cannot rename (see the
  Session compatibility reference for the other saved `hash=` values, which
  load as saved). `hash=` values that
  you wrote yourself for a feature of a cropped or reverse-complemented record,
  in a TSV file, a CLI option, or a table typed into the Web app, need that
  feature's new hash: the one its newly drawn SVG ID carries.
- Sessions: a **This feature only** color that 0.13.0 saved on a multi-record
  Linear diagram colors its feature again after Load and Generate (OV-416).
- Gallery: the majanivirus example shows Litopenaeus vannamei majanivirus
  (LC738872.1) reverse-complemented with its record's **Reverse complement**
  setting, so its links to the neighboring genomes no longer cross. The rows stay
  centered, and the Gallery command adds `--reverse_complement 1` for that record.
- Records (web app): choose which records of a file to draw. Linear record cards
  gain **Draw this record**, and a File with several records gains **Choose
  records…**, a list that searches and sorts by file order, record ID, or length
  and turns records on or off in bulk. **Choose records…** also appears beside
  **Record** in Circular, where **Multi-Record Canvas** and separate diagrams
  draw only the records that are on. A record that is off is left out of
  **Generate Diagram**, comparisons, row layout, and the Legend, keeps its
  settings (also in a saved Session), and returns in its place when turned on;
  **Delete settings** removes its card settings, feature edits and Feature
  placements, annotations that target it, and comparison pairs that name it.
  LOSAT keeps the whole subject file as its E-value database. After an upload, a
  file with more than 20 records opens the list, and turning off the last drawn
  record of a file asks whether to remove the file. Before, Linear drew every
  record of a file and Multi-Record Canvas drew every record. Session version
  46 gains the optional `recordsOff` field.
- Sessions (CLI): a CLI re-save of a Web Session no longer drops **File
  definition**, **File subtitle**, or the inferred definition from Linear cards
  (OV-402).
- Similarity alignment (web app): the **Select alignment anchors** dialog is now
  **Review alignment**. One sentence that names the reference gene and an
  **Apply required** label replace the three introductory paragraphs; the
  record counts are on the section headings and in the panel's description
  for screen readers. **Arrow direction** offers **Keep as is**, **All right →**, **All
  left ←**, and **Custom**, with one line that says which records reverse. Each
  record is a card with its organism, accession, and current arrow; a record
  that Align reverses shows the old arrow struck through, the new arrow, and a
  **reversed** chip. The reference card says how far its record's left edge
  moves. A record with one usable candidate shows a compact card with a
  **Skip** checkbox; records that still need a choice come first under **Choose
  a gene**, then **Several candidates**, **One candidate**, and **No usable
  candidate**, and the footer names how many records block **Apply**. The
  reference card and the one-candidate cards give the gene's coordinates and
  strand. The panel is
  docked in the Result, beside the figure or below it, at every width and
  without dragging, instead of floating over
  it: the figure is fitted to the room that remains, the Editor and the feature
  popup step aside, and all three return when the panel closes. Clicks and
  selections on the figure leave the hidden popup alone; **Open** in the
  feature search row makes its feature the popup's, shown when the panel
  closes. A one-click **Align…** takes nothing over. Each card draws
  its record as it will lie after **Apply**, and on the figure a reference
  line, a ring on each candidate, and an arrow replace the numbered badges; a
  click on a ring chooses that gene. The
  **Recommended** chip stays on the recommended candidate after you choose
  another. Behaviour is unchanged: choices stay local until **Apply**, which
  validates, regenerates, and commits in one step.
- Similarity alignment (web app): while **Review alignment** is open, the
  figure is fitted again whenever its room changes: the window is resized, the
  settings panel is dragged wider or narrower or switched between the full
  panel and the icons, a section beside the icons opens or closes, the review
  panel moves between the side of the figure and below it, or a notice above
  the Result appears or goes. Before, the figure kept the fit of the moment the
  panel opened, so it was cut or left room unused. A zoom or position you set
  by hand during the review stays while the room stays and is replaced by the
  fit at the next such change; **Fit** works as before, and the view from
  before the review still returns when the panel closes (OV-536). The settings
  panel now follows the pointer from the first move of a width drag; it jumped
  10 px wider first (OV-620).
- Similarity alignment (web app): changing a feature colour, **Show**/**Off**,
  or a **Legend** or **Plot Title** setting while **Review alignment** is open,
  or undoing such a change, no longer ends the review's rings and fitting: the
  rings are drawn again, the figure is still fitted when its room changes, and
  the view from before the review returns when the panel closes (OV-1082,
  OV-1120). Apply keeps a Legend or Plot Title setting changed during the
  review. Switching between Circular and Linear and back
  with the panel open no longer loses the view from before the review (OV-622).
  Apply now keeps a colour changed before or during the review; before, it drew
  the old colour or failed with "The operation failed without recognized
  diagnostic information" (OV-1083).
- Similarity alignment (web app): the **Reset alignment** dialog says "Undo
  the latest Align." and offers **Positions only** and **Positions and
  directions (N records)**. The reason a scope is unavailable sits under it.
  The per-record list, one line each (current → restored direction), appears
  only when the directions scope is selected, and the primary button names the
  scope (**Reset positions** or **Reset positions and directions**).
  Behaviour is unchanged.
- Similarity groups (web app): in the Editor's **Similarity groups** tab the
  active alignment is one line ("Aligned to neoU in BGC0000709" with **Reset
  alignment…**) that opens to show the reference and records. Group rows have
  two lines (name, "5 genes · 5 records"), renamed groups carry an **Edited**
  badge, and the selected group's detail and **Align…** buttons scroll in one
  region under the list. The exact-reference options name the record and the
  gene instead of the record key and feature hash. **Reset Rename** is now
  **Reset name** beside the group's ID; it still restores the name and the
  description. As in the popups, one **Copy** menu holds the group's FASTA
  copy and download, each member has a **FASTA** menu instead of four
  buttons, **Description** shows only a description you wrote (the computed
  one is under **Name suggestions**), and the list scrolls to the group that
  a new alignment selects.
- Similarity groups (web app): the disabled **Similarity groups** tab, for
  example in Circular, says in its tooltip that Similarity groups appear after
  a Linear comparison with Similarity groups turned on. The group's
  **Display name** and **Description** fields are now named by their labels
  for assistive technology (OV-427).
- Similarity groups (web app): member tables in the feature popup, ribbon
  popup, and Similarity Groups drawer name each record in plain text, as the
  diagram heading does. Before, a definition with markup such as
  `<i>Streptomyces</i>` showed the tags, and a File default definition was
  ignored (OV-369).
- Collinear blocks (web app): a collinear ribbon pairs its genes by group and
  shows each group's members again. Before, a ribbon without a collinear group
  scope looked its groups up only among scoped collinear groups and showed only
  the group IDs (OV-370).
- Collinear blocks (web app): the block popup leads with one line: both spans
  with their lengths in kb, the orientation, the identity rounded to a whole
  percent, and the number of gene pairs, with **Copy** for that line and one
  **FASTA** menu for the spans. **N gene pairs** opens one table of query gene,
  subject gene, and group, replacing the separate group and Query/Subject
  tables. Block ID, kind, color mode, score, anchors, the exact identity, the
  aligned length in aa, and the note on block envelopes are under **Details**.
- Similarity groups (web app): a group ribbon popup is titled with the group
  name, shows the group ID as a badge and one **N genes in M records** line, and
  lists the members with one **Copy** menu for the group and a **FASTA** menu
  per member. RBH seeds, group paths, and a computed description are under
  **Technical details**; a description you wrote stays under the name. A
  member's role and confidence are shown only when they are not anchor and
  high.
- Feature popup (web app): the **Similarity group** section is one line (group
  ID, name, **N genes in M records**) with **Show in editor** and **Align…**;
  **Review alignment options…** is a link under **Align…**. **Highlight**,
  **Clear highlight**, and the group's FASTA copy and download are under
  **More**, and the member table opens from **N genes**. The **Current
  member** line and the consensus sentence are gone; the current member's row
  is marked in the table.
- Pairwise match popup (web app): a LOSATP pairwise ribbon whose genes are in
  no similarity group no longer shows a **Similarity group** section holding
  the query's product as its name (OV-420).
- Similarity alignment (web app): while an alignment is active, including one
  loaded from a Session, the Similarity Groups drawer selects its group and
  preselects its exact reference, so **Align…** and **Review alignment
  options…** work from the drawer at once. Before, the drawer showed the
  active alignment's reference but kept both buttons disabled until you chose
  the same reference again.
  A new alignment, an Undo, or a Session load replaces an earlier drawer
  group and reference choice with the new plan's (OV-381).
- Similarity alignment (web app): an Undo that restores a cleared alignment
  also removes its "Alignment cleared" notice. Before, the notice stayed while
  the alignment was active again (OV-382).
- Similarity alignment (web app): choosing **Upload BLAST TSV** for a pair
  before choosing its file keeps the active alignment, because the pair has no
  new data yet. Before, the alignment was cleared at once with "Alignment
  cleared: comparison configuration changed." Choosing a file, another
  program, or another LOSATP mode still clears it (OV-380).
- Similarity alignment (web app): a Result filter or grouping setting that
  moves the alignment reference or a target out of its Similarity Group (for
  example a higher **Minimum identity**) keeps the alignment, and every later
  **Generate Diagram** draws it while its anchor features exist. Before, the
  Generate after such a change was blocked with "Alignment needs repair". When
  an anchor feature itself is gone, the repair notice now says whether the
  reference or an aligned target feature is gone instead of "The clicked
  feature is not a unique member of this Similarity Group." (OV-341, OV-342).
- Feature popup (web app): while the first popup after a page or Session load
  waits for Python to match the color rules, the toolbar shows **Preparing
  feature details…**, and **Save Session** and **Load Session** stay available.
  Before, the toolbar said "Applying an edit. Retry after the edit finishes.",
  both buttons were disabled, and the click seemed to only highlight the
  feature's similarity group (OV-377, OV-378).
- History (web app): Undo and Redo of a Linear setting keep the active
  similarity alignment and record positions. Before, undoing any Linear
  setting, for example the first edit after loading an aligned Session,
  silently dropped the alignment, and the next Generate drew the records
  unaligned (OV-383).
- Gallery: the aminoglycoside BGC example draws its arrows with **Shaft Width
  Ratio** 0.6 and aligns its records on similarity group og_6 (neoU of
  BGC0000709) instead of og_1. Its tutorial adds the **Features › Arrow Geometry
  › Shaft Width Ratio** step and walks through **Review alignment options…**,
  where BGC0000708 has two og_6 members and keeps the recommended livU.
- Gallery tutorials: the BGC, Hepatoplasmataceae, and majanivirus tutorials now
  say that the LOSAT runtime controls are under **Comparison › Settings ›
  Runtime and reproducibility**, not under **Advanced comparison and layout**
  (OV-366).
- Gallery tutorials: the BGC and *Vibrio harveyi* group tutorials now teach in
  stages. The reader adds one file at a time and fills its **File defaults**,
  clicks **Generate Diagram** with the web app's defaults, then changes only
  what the preview shows is wrong and generates again. Each intermediate
  preview is shown. LOSAT runs once, and later Generates reuse its results.
- Gallery: the aminoglycoside BGC example's Session is rebuilt in the web app
  from the five MIBiG files. BGC0000713 is the original record with **Reverse
  complement** on, so its popup coordinates are those of the source file. The
  figure now draws what its command gives: the orange palette's blue identity
  ribbons and CDS default #dddddd.
- Linear ruler (Python API and web app): with **Ruler on Axis** active, the
  ruler's ticks and labels follow the axis color unless **Scale Line Color**
  or **Ruler Label Color** is set, as the CLI already drew them; before, they
  were black (OV-440). When **Ruler on Axis** is ticked but does not apply
  (**Bar** scale style or **Middle** track layout), the axis keeps the plain
  Linear color instead of the ruler's darker one, as in the CLI (OV-442).
- Gallery tooling (maintainer): `tools/capture_gallery_tutorial_screenshots.py`
  captures an intermediate tutorial state by replay. A capture's
  `replayThrough` loads the example Session, clicks **Reset Settings**, then
  performs the `replay` actions of the tutorial steps up to the named one.
- Web: In Colors › SPECIFIC RULES, color rule changes made in quick succession
  are all kept, each as its own Undo step. This covers **Add**, color and
  caption edits, Move, Remove, a popup or Legend recolor, a rules TSV import
  and a preset. Before, a change made while an earlier one was still being
  applied could be silently dropped. A pattern edit that is still being
  applied when another change is made stays in its field as an unapplied
  draft, as before (OV-445, OV-443, OV-444, OV-448, #1006).
- Web: **Load Label TSV** now says when it was not applied because something
  changed while the table was loading: a label or feature visibility edit
  (which is kept), a replaced diagram, or changed colors, color rules, Legend,
  or inputs. Load the table again (OV-447, OV-449, #1008).
- Web: an undoable action started while an earlier one is still being applied
  (for example a record display edit while **Load Label TSV** is still loading
  its table) now applies at once and is its own Undo step under its own name,
  and what the earlier action writes afterwards is an Undo step under that
  action's name. Before, the two were one step named after the earlier action,
  and what the earlier action wrote after the later one was lost on Undo then
  Redo. The steps follow the order of the writes. As before, a change made
  meanwhile in a plain control (a text, number, or color field, a checkbox, a
  radio button, or a select) or by a Legend or diagram drag goes into the
  earlier action's step, and so does any edit made while a **Specific Table
  (-t)** upload, **Reset Settings**, or a color rule change that adds or
  removes a Legend row is being applied (OV-460, OV-461).
- Changed (web app): Undo and Redo of a color change that adds or removes a
  Legend row no longer draw the diagram again in Python once the automatic
  redraw of that change has been shown while the change was the last step:
  its History step keeps the redrawn diagram. An Undo or Redo made while that
  redraw is still running, or after a redraw that never finished, still draws
  again, so the diagram and its Legend always show the restored colors
  (OV-731).
- Fixed (CLI, Python): an integer written with non-ASCII digits (for example
  `#١` or `#²` as a record selector, or `-w ٣` as a window size) is refused
  with the reader's usual error. This covers record selectors, regions, record
  rows, table cells, Session slot indexes, track slot `z=`, GenBank
  `/codon_start` and `/transl_table`, GFF3 `phase`, and every integer option
  of the `gbdraw` commands. Before, some of them read the value as the ASCII
  number and others ended in `invalid literal for int()`. An integer written
  with underscores (`1_0`), or with more than 4300 digits, is refused too;
  signs and surrounding spaces are still accepted (OV-674, OV-676, OV-960,
  OV-961).
- Fixed (CLI, Python): the same rule now covers decimal numbers. A number
  written with non-ASCII digits or with underscores is refused with the
  reader's usual error, where Python's `float()` read it as the ASCII number:
  for example `--gc_content_radius ٠.٥`, `--evalue 1_0`, `--track_axis_gap ١`,
  a track slot width `٢px`, a config or Session depth bound, or an annotation
  table `stroke_width` of `1_0` (the table's `start` and `end` integers too).
  This covers every decimal option of the `gbdraw` commands, the decimal
  readers of Sessions and tables, and every number, integers included, of a
  config file or of the `config` of a Session or a Python API call
  (`radius: "٣٩٠"` drew a radius of 390). Everything else `float()`
  accepted is read as before: `0.5`, `+.5`, `5.`, `1e-5`, ` 2 `, and `inf` or
  `nan` where a reader allowed them (OV-962, OV-1006).
- Fixed (web app): **Load Session** of a Session 31 to 33 or 39 whose saved
  conservation series has an empty (`null`) row stopped with an unknown error.
  That row is now read like a row without a label (OV-963).
- Fixed (web app): a color change made while the automatic redraw of an
  earlier color change is still running is one Undo step again. Before, the
  Legend rows that redraw drew (or, after a Legend or title move, the moved
  position under the redrawn diagram) were recorded as an extra step under the
  same name, whose Undo changed nothing visible (OV-872, OV-873).
- Fixed (web app): with LOSATP Collinear and **Search between** set to **All
  records** on a CLI Session's record grid that places a record above one
  listed before it, **Save Raw LOSAT TSV** on one comparison card saved the
  search of another pair, and the other pair's card had no **Save Raw LOSAT
  TSV** button. The reproducibility files of two such pairs were also numbered
  alike, so one overwrote the other. Each card now saves the search of its own
  two records under its own file name, and each pair has its own
  reproducibility file (OV-980, OV-985).
- Fixed (web app): a Session whose raw LOSAT cache entry has a file name that
  is not text (for example `{}`) loads with the entry's default name,
  `losat_pair_<n>.tsv`. Before, **Save Raw LOSAT TSV** named the file
  `[object Object]` (OV-981).
- Fixed (CLI, Python): a Session whose raw LOSAT cache entry has a file name
  that is not text (for example `7` or `["a.tsv"]`) is read as having no name,
  as the web app reads it, instead of writing that value as the name. A Session
  saved again names such a shown entry `losat_pair_<n>.tsv`, the name the web
  app gives it, and writes a hidden (`display: false`) entry without a name, as
  the web app does; before, the CLI kept a hidden entry's name (OV-986,
  OV-1007, OV-1013).
- Fixed (web app, CLI, Python): in an interactive SVG, a feature or orthogroup
  member whose start is missing (`null`) no longer shows as starting at 1. The
  web app wrote such a start as 0, so the feature without location parts read
  `1..20`. Both writers now read a feature's or member's start and end as they
  read a location part's: a number or ASCII number text, truncated toward zero;
  blank text, a boolean, or a list is missing (OV-1009, OV-1011).
- Changed (web app, interactive SVG): the web app and the interactive SVG
  viewer read the start or end of a location part, a feature, or a Similarity
  Group member by the same rule as the two writers. This concerns hand-edited
  Session or SVG data only; gbdraw writes whole numbers. A fraction is
  truncated toward zero: a part from `1.5` to `20` shows `2..20` (was
  `2.5..20`), its sequence is extracted, and Align accepts such a member
  (it refused it with "invalid display facts"). Hex or binary text (`0x10`),
  text with any non-ASCII character (a no-break space included), and whole
  numbers beyond 2^53 - 1 are missing: the location, length, FASTA header
  coordinates, and sequence are left out for them and Align refuses such a
  member, where `0x10` read as 16 before. The member rows of a Similarity
  group in the feature popup and the Similarity drawer show the location as
  the member popup does: a missing start no longer shows as `1..20`, and a
  start of `1.5` shows `2..20` (was `2.5..20`) (OV-1015). Rotating a record to
  a feature reads the feature's location parts by this rule too, so a part
  bound stored as text (`"20"`) is no longer refused as outside the record
  (OV-1016). The Interactive SVG viewer shows a member with one missing
  bound the same way as the app (OV-1017).
- Fixed (web app): loading a Session older than version 46 that recolored a
  Legend row named like a built-in JavaScript member (`constructor`,
  `toString`) without a saved original color no longer records that member as
  the swatch's original color. Before, Reset set the swatch fill to
  `function Object() { [native code] }` (OV-982).
- Fixed (web app): Load drops a raw LOSAT cache entry whose `searchContext` is
  not text (for example a list holding the 64-character hash), as it drops
  other entries it cannot read. Before, the entry was kept, never reused, and
  written again by Save Session (OV-983).
- Fixed (web app): when a Session is loaded, a LOSATP Collinear comparison card
  whose search ran in the other direction than the card (a CLI Session's record
  grid that places a record above one listed before it) keeps its **Save Raw
  LOSAT TSV** button, both for the CLI Session itself and after Save Session
  and Load. Before, the button was missing until the next Generate (OV-987).
- Fixed (web app): a Session with an unsafe key (`__proto__`, `constructor`,
  `prototype`) anywhere in the file, for example under `resources`, is refused
  with "An input value is invalid" and the advice to load a supported Session
  file. Before, it was refused as a failure without recognized diagnostic
  information (OV-952).
- Fixed (web app): when a Session saved by 0.13.0 after a mode switch stores a
  Legend side outside its Result's mode, the Load dialog lists it as the side
  of a saved Result ("Legend of a saved Result: ..."), as it does for the same
  value in a Session saved without a mode switch. Before, it was listed as a
  Legend setting ("Legend: the position ...") (OV-964).
- Fixed (web app, interactive SVG): a feature location part whose start or end
  is missing (null or blank) is read as missing, not as 0. The interactive SVG
  viewer and the web app's feature details show no location or length for
  such a part, the FASTA header names no coordinates for it, and neither takes
  its sequence from the record start. Before, the location began at 1 and the
  sequence at the record start (OV-968).
- Fixed (web app, interactive SVG): a Similarity Group member whose start or
  end is missing (null or blank) shows no coordinate for it in the group
  popups of the web app and the interactive SVG, gets no `1-20(+)` location as
  its display ID, and makes Align, and Generate with a saved alignment of
  that Similarity Group, refuse the whole group with "Similarity Group member
  N has invalid display facts." (at Generate, after "Alignment needs
  repair: "), as an undefined start already did; this concerns hand-edited
  Session data only. Before, the missing value was read as 0 (OV-1010).
- Web: when **Load Label TSV** is applied while the Label Not Shown dialog is
  open, the dialog is cancelled and the edited label text is not written.
  Before, the dialog disappeared without an answer, and Undo, Redo, Save and
  Load could stay unavailable until the next undoable edit or Generate
  (OV-462).
- Gallery publication (maintainer tooling): `tools/refresh_gallery_sessions.py`
  publishes a Session that the web app saved. It rebuilds the request from the
  request's own color, whitelist, and priority tables rather than from the
  files bound to those inputs (OV-367). It stores Default colors as the CLI
  replay resolves them (OV-375) and one-record GenBank sources as `exactly_one`
  (OV-376). It also names the example's GenBank files in `examples.json` when
  the Session has no CLI invocation (OV-397).
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
- Editing (web app): an edit to a Feature visibility rule made while a label rerender is running is drawn as **Generate Diagram** draws it; before, the shown Result could keep the older rule until the next Generate (#987).
- Editing (web app): a color edit no longer adds a Legend row live for a color
  rule that no feature drawn on the shown Result matches, such as a rule for a
  feature of the other Result of a batch or a feature that a Feature visibility
  rule hides; before, the row showed until the next **Generate Diagram**
  (OV-344).
- Web app: **Load example** is now also in the header, after **Load Session**, so a
  second Gallery example can be loaded once a Result is shown. It opens the same chooser
  as the **Load an example** button in the empty Preview, which stays. At 769 to about
  830 px the header no longer runs past the viewport, and its button labels no longer
  wrap onto two lines up to about 1300 px (OV-368, #986).
- Sessions (web app): Session Load keeps Specific-colors rules in the shape Generate
  stores them, and keeps a Linear row without a Depth file as saved (`null`, not
  `[null]`). Before, a Session saved after Generate and then loaded did not give back
  the same editor state. What the app draws does not change (OV-334, OV-337, #984).
- Sessions (web app): the Session Load notice quotes a dropped Legend entry's caption
  with double quotes, so a caption such as `5' UTR` reads clearly (#983).
- Editing (web app): every editor edit (stroke, fill, palette, color rule, feature
  visibility, Legend rename, recolor, delete, Restore, and sort) now reaches the
  displayed Result through the compile and Result executor that Generate uses, so the
  live Result shows what the next Generate draws. This fixes Reset Stroke, Reset all
  strokes, and Undo of a stroke for connectors and every size class; a batch Result
  shown again after an edit was removed; Undo, Redo, and a batch Result display after a
  Legend color, color rule, or palette step; and Feature Edits TSV loads and popup label
  edits with Label Rendering = Embedded Only (OV-129, OV-144, OV-150, OV-200, OV-293). A Legend whole-row rename,
  merge within one feature type, recolor, delete, Restore, and sort show live without
  the automatic rerender when Python would draw the same rows. Deleted Legend rows stay
  in the Result, hidden; exports strip them (#975).
- Editing (web app): **Reset Settings** returns every Legend rename, including track
  rows such as `GC content`, to the default caption, and Undo of Reset restores the
  renames. A Legend rename onto a deleted row's caption opens the Legend Name Conflict
  dialog. When both rows draw one feature type it offers Restore the deleted row and merge
  into it, Keep current color and add a suffix, or Cancel; rows of different types get
  the suffix or Cancel. A Legend row stroke follows the color Python gives the row's
  features, so it stays after a Legend row color edit in either order (OV-288). Rule
  rows are addressed by Python's rows (OV-292, OV-294), and a palette row's Legend color
  is read from the color the row shows now, not from the last Generate (OV-282)
  (#975).
- Removed (web app): **Add legend item** is retired. No control called it, and it was
  reachable only from the browser console. Sessions with editor-added Legend rows still
  load and show them, and the popup reset dialog's "This only + add legend entry" is
  unchanged (#975).
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
  rule captioned `rRNA` on tRNA features, keeps its own color row instead of
  being dropped; the two rows are numbered `rRNA (1)` and `rRNA (2)` (N-06).

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
- A legend row that a specific-color rule draws under the name of a generated
  row in another color (the two are numbered, such as `rRNA (1)` and
  `rRNA (2)`) now edits that rule in the Legend editor, and a live rule edit
  keeps that row instead of recoloring or failing on the generated row (N-06).

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
- Placement (web app, CLI): the Feature placement table can name one of two identical
  features (same record, type, and coordinates). `hash=` also compares the complete
  biological feature ID such as `fe6f54094~1`, and the Run Info source recipe writes
  `feature_selector=hash=<biologicalFeatureId>`, so the CLI accepts the recipe. A value
  without `~` keeps its meaning (#766).
- Generate (web app): a Generate that fails before it replaces the displayed Result now
  ends with the classified **Generation Error**, the previous Result stays, and Undo and
  **Save Session** work again. Before, the page stayed on "Generating Diagram..." until
  reloaded. The History rollback of **Rotate record to feature** uses the same restore
  (#785).
- Labels (web app): the Label Not Shown dialog and the feature popup note state the
  actual reason a label is absent (hidden, Underlay, Embedded Only, Show Labels None or
  First, the whitelist, or the blacklist) instead of one sentence about Show Labels and
  the label filter (#799).
- Placement (web app): every control that changes the direction of the feature slot
  (**Use custom stack**, the custom stack resets, adding, deleting, disabling, moving,
  duplicating, or re-typing a stack row, and the Circular **Separate Strands**) opens
  the "Reset Feature placements?" dialog when a lane placement would become undrawable.
  **Reset N placements to Auto** is one Undo step, and **Cancel change** keeps the
  value. Before, they changed the slot silently and the next Generate stopped on
  `FEATURE_PLACEMENT` (#805).
- Annotations (web app): an annotation made from a selected feature is kept in the mode
  it was selected in. Before, a Circular selection could be drawn on the Linear diagram,
  and the reverse, when both used the same record key. Annotation sets stay shared by
  both modes (OV-21, #806).
- Linear (web app, CLI): a multi-record Linear diagram with a cropped record that has
  two features at the same coordinates no longer stops at result admission or at
  `--session_output` with "Feature metadata identity does not agree with rendered SVG
  ID" (#808).
- Sessions (web app): loading a Session 40 to 44 moves an annotation target
  `hash=<hash>` that a selection made onto its feature, where the figure cannot change,
  so it stays on its feature after a crop or reverse complement. The Load notice reports
  how many targets moved (#811).
- Sessions (web app): loading a Session before 40 that has feature edits but no saved
  feature catalog now fails with the `WORKER_INIT` diagnostic when the diagram runtime
  cannot start, and the previous Session stays. Before, Load reported the edits as
  dropped, and the next Save lost them (OV-39, #833).
- Sessions (web app): a Session 27 to 30 saved with Custom Track Slots off, such as the
  0.13.0 Gallery Session `BGC0000708-BGC0000713`, loads. Another retired slot field
  still fails Load, and the message now names the field and the track row (OV-38, #834).
- Sessions (web app): a Session 31 to 39 whose Default colors, Label whitelist, or
  Qualifier priority table has a row with extra cells now loads. Extra cells join the
  last column with one space, a row that lacks a required column is dropped, and the
  Load notice names each changed row. A remaining table error names the table and says
  it comes from the Session (OV-40, #835).
- Live edit (web app): a Generate that succeeds removes the "Live edit failed: ..." note
  above the Result. A failed, canceled, or superseded Generate keeps it (OV-36, #838).
- Legend (web app): in a Circular batch, a second **This feature only** fill or stroke
  on a feature colored at the last Generate no longer fails Generate when another Result
  does not draw that row. A missing Legend row that the plan requires now reports
  `RESULT_INVALID` instead of the generic error (OV-45, #845).
- Sessions (web app): after loading a Circular Session, the Record control reads
  "NC_012920.1 (not inspected)" until the source records are inspected, and says "not
  found" or "ambiguous" only after an inspection (#852).
- Legend (web app): a live visibility or specific-color rule edit that changes what the
  Legend derives from (drawn feature types, their first-drawn order, or the rules the
  drawn features use) asks for the automatic rerender, with Auto Reflow on or off, so
  the Legend of every Result matches Generate (OV-44, #857).
- Legend (web app): a Legend-only color or stroke on a row that only one Result of a
  batch draws no longer fails Generate. Result admission requires the row in at least
  one Result (OV-46, #871).
- Generate (web app): an automatic rerender never starts while Generate runs. Before, it
  could take the newer generation token, and Generate ended with no Result and no
  message. A rerender requested meanwhile runs once after Generate (OV-48, #879).
- Legend (web app): recoloring a Legend row, or its rule in the Specific color rules
  table, after coloring that row's features in the feature popup keeps the new color at
  the next Generate. Before, the Legend swatch came back to the popup's color (#887).
- Errors (web app): a preserved setting that belongs to the other diagram mode, and a
  stale Reset alignment, now fail Generate with messages that name the setting or the
  evidence and the ways out, instead of generic text (OV-130, #925).
- Performance (web app): the popup feature color change prepares only the clicked
  feature's rules, so it no longer hangs or crashes on Sessions with 6,000 or more
  features (#949).
- Performance (web app, Python): popup "same Legend item" color choices and stroke edits
  on large Sessions do work linear in features plus rules, and Generate validates LOSATP
  manifests once. On a 6,304-feature Session the color choice went from over 60 s to
  12.9 s. SVG bytes are unchanged (#952).
- Search (web app): the active match line names its record, such as `Record NC_022349.1
  · Matched Qualifier gene: dnaA`. The sidebar buttons **Use selected feature 5′ end**
  and **Use selected feature midpoint** are removed; start a record at a feature with
  **Rotate record using this feature** in the feature popup (PD-OI-087, #954).
- Performance (web app, Python): large Sessions pay repeated copies, validations, and
  per-feature Worker messages less often on Load, Generate, rule commits, and `gbdraw
  <mode> --session`. SVG and Session JSON bytes are unchanged (#956).
- Gallery: the Vibrio collinearity example starts every chromosome at its replication
  initiator, the other nine examples have titles that name the result, `docs/GALLERY.md`
  mirrors the Web Gallery, and the label settings stay visible while a feature label is
  On (OV-222, #957).
- Web app: the popup color scope, Legend name, and **Reset fill color** dialogs open at
  once and show **Applying an edit…** until the choice is applied. Opening them no
  longer redraws the Specific color rules rows. gbdraw checks every user color once when
  it reads it, and an unreadable color stops with a gbdraw error (PD-OI-088, #958).
- Gallery: the Vibrio Gallery Session stores its finished comparison result, so the
  first Generate after Load takes about 20 to 29 s instead of 37 to 39 s. The figure is
  unchanged (#961).
- Colors (web app): a Default colors value that differs from the selected palette is a
  user color, as `-d` over `-p` on the command line. A palette switch asks before it discards user
  colors (**Keep my N colors**, **Use the palette's colors**, or **Cancel**); Default colors
  **Reset** asks first (**Reset to the palette's colors** or **Cancel**), and **Apply to all "<type>"** on a palette Legend
  row sets the type's default color instead of one hash rule per feature (PD-OI-089,
  #963).
- Colors (web app): choosing **Apply to all** in the popup's Color Change Scope dialog
  while an automatic rerender runs now applies after the rerender, with the dialog
  showing **Applying an edit…** (OV-280, #969).
- Reset (web app): after **Reset Settings** the displayed Result shows the reset Legend
  strokes, Legend row colors, deleted rows, feature fills, hidden features, and label
  edits at once. Load then Generate of the Vibrio Gallery Session draws the published
  figure with V. alginolyticus chromosome II starting at `display_start 430232`.
  **Save** after **Load** of a CLI-written Circular Session that gives two records of
  one file a display start no longer fails, and the Linear comparison **Result filters**
  say they apply on Generate (OV-325, OV-332, #977).
- Performance (web app, CLI, Python): each render request serializes the SVG once
  instead of twice, the web app checks the records it keeps between Generates feature
  by feature once instead of twice, and a Linear plan projects its comparison rows once
  instead of twice. SVG, PNG, and PDF bytes are unchanged. A Linear request with an
  invalid comparison binding and also an invalid annotation selector or alignment plan
  now reports the annotation or alignment error first (#1045).
- Performance (web app): a Generate after a small change reuses the parts of the render
  whose inputs did not change (color, visibility and label-filter tables, GC, skew and depth
  data, the legend table, the interactive context, the feature catalog, and Linear comparison
  rows, protein extraction and feature layers). SVG, PNG, and PDF bytes are unchanged on
  every surface. `gbdraw.api.apply_config_overrides` now shares `labels.filtering` with its
  input when no override touches it. A collection region with a file selector no longer
  reuses the crop of another file inside one cache (OV-700, #1055).
- Performance (web app): a Generate after a small change also reuses the legend and plot
  title contents, the composition plan, the contents of every Linear group, the feature
  catalog item of a Linear Result and the protein comparison analysis. SVG, PNG, and PDF
  bytes are unchanged. The internal record annotations `gbdraw_record_label` and
  `gbdraw_record_subtitle` are no longer written; labels come from the request.
  `gbdraw.io.comparisons.load_comparisons` is removed (#1059).
- Performance (web app): a Generate after a small change on a Circular Result also reuses
  the record features, labels and group contents, the Multi-Record merge and the feature
  catalog item. A title change on a large Multi-Record Session drops from about 1.4 s to
  about 0.1 s warm (native Python). SVG, PNG, and PDF bytes are unchanged (#1070).
- Performance (web app): a Generate whose feature catalog was not built again keeps the
  catalog the page already holds: the catalog is not sent, decoded, or admitted again.
  Python returns it beside the metadata, not inside it, and the Worker render protocol is
  now 3, so a cached older wheel is reported as a version mismatch. The figure, the
  catalog, and Saved Sessions are unchanged.
- Performance (web app): a Generate sends and sanitizes only the SVG groups whose content
  changed. Python writes each SVG document as ordered pieces with a sha256 revision for
  each group's content; the page keeps each group's sanitized text on its Result, and the
  next Generate names those revisions so Python leaves their text out. The joined Result
  equals DOMPurify of the whole document. Python (CLI and Python API too) now raises
  `ExportError` for an SVG group that holds a comment, CDATA or a processing instruction;
  gbdraw writes none. The Worker render protocol is now 4. The figure, SVG, PNG and PDF
  bytes, and Saved Sessions are unchanged.
- Performance (web app): after a Generate that keeps the feature catalog, the page no
  longer derives three values again whose inputs did not change: the feature set the
  color and visibility rules match against, the record and feature names of the
  Similarity alignment drawer, and the Features list. Nothing visible changes.
- Performance (web app): after "Preparing render inputs and session..." and "Preparing
  preview..." a Generate no longer waits two animation frames; it yields one task, so these
  two statuses may not be painted or announced on a fast Generate. The first status still
  paints before the work starts. A Cancel clicked during either step is still honored there
  in Chromium: no Worker restart, and the canceled Result is never shown (PD-OI-101).
- Performance (web app, CLI, Python API): a Generate spends less time in Python. On
  every surface the request contract check is compiled once per type and the overlay
  Legend search stops at the first conflicting obstacle. In the web app, a linear
  Generate no longer builds the derived LOSAT cache entries that no reply reads (their
  input checks still run), the Circular definition's GC percent and the placed labels
  are keyed once and reused, and the cache no longer checks render-stage entries against
  a constant. Figures, SVG bytes, replies and Saved Sessions are unchanged.

- Web app, Custom Track Slots: each row keeps its Move up and Move down arrow
  buttons, and one **Row actions** menu holds its other actions with text
  labels: the moves across the Axis, **Track on Axis** where the row can sit on
  the Axis, **Duplicate**, and **Remove** (Circular) or **Delete** (Linear), in
  red. The menu closes after an action, the arrow keys move through it, and
  Escape closes it. After a move, focus stays on the row's **Row actions**;
  after a removal it moves to the next row's, else the previous row's, else
  **Add track**.
- Web app: the Session title is asked in the app's own dialog instead of the
  browser's prompt: **Save Session** without a title opens **Save Session**,
  and the title in the header opens **Rename Session**. The field starts with
  the title selected; Enter saves or renames, Escape cancels, and focus returns
  to the button that opened it.
- Web app (OV-1097): a Region **End** past the record end is no longer refused.
  Generate draws to the record end, as the CLI and Python already do, and a
  note under the Region fields (Circular, and Linear once the record is read)
  names the record length. A **Start** past the record end stays an error.
- Web app, Custom Track Slots (OV-1099): **Duplicate** in Linear names the copy
  as Circular does, with the renderer's default ID and the first unused number
  (`gc_skew_2`, then `gc_skew_3` for a copy of either), instead of appending
  `_copy`. A copy is named from its renderer's default ID, also when the source
  row has an ID of its own. **Duplicate set** in Region Annotations follows the
  same rule (OV-1139): the copy takes the ID **Add set** gives (`annotations`,
  `annotations_2`, ...) instead of `<id>_copy`. Saved slot and set IDs are
  unchanged. A Linear copy of an unedited Depth row is a managed Depth row, as
  in Circular: removing the Depth files removes it with its source row.
- Web app (OV-1137): the **FASTA** menu of each similarity group member, in the
  feature popup, the pairwise match popup, and the **Similarity groups** drawer,
  has an accessible name that names the member, such as "FASTA of BDT61648.1",
  the same text as its tooltip.

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
