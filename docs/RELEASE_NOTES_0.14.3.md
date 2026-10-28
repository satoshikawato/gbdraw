[Documentation home](./DOCS.md) | [Installation](./INSTALL.md) | [Changelog](../CHANGELOG.md) | [0.14.2 release notes](./RELEASE_NOTES_0.14.2.md) | [0.14.1 release notes](./RELEASE_NOTES_0.14.1.md)

# gbdraw 0.14.3 release notes

gbdraw `0.14.3` is a patch release of [0.14.2](./RELEASE_NOTES_0.14.2.md). In
the Web app, per-feature edits keep applying when the view of their record
changes, and redraws after small changes do less work. It also makes Circular
label placement much faster on large label sets, stops Circular labels that
cannot fit before spending minutes on them, makes the command line fail
clearly when a requested image format is unavailable, and renames the palette
`hearbal_garden` to `herbal_garden`. A new Circular **Show Definition**
setting hides the center definition. gbdraw 0.14.3 saves Sessions as session
version 49 to store that setting (the canonical `renderRequest` schema stays
9): Sessions saved by 0.14.2 open as before, and gbdraw 0.14.0 to 0.14.2 cannot
open a Session saved by 0.14.3.
[GitHub Releases](https://github.com/satoshikawato/gbdraw/releases) gives the
publication date, [PyPI](https://pypi.org/project/gbdraw/) lists the published
packages, and [Installation](./INSTALL.md) explains each install route.

## Fixed

- **Per-feature edits follow their record in the Web app.** A "This feature
  only" color, label text and visibility, Feature visibility, stroke, Feature
  placements, and annotations of selected features keep applying when the
  record is drawn another way: a Circular record drawn alone or on the
  Multi-record canvas, and a Linear File card that draws all of its records or
  one of them. Before, they stopped applying after a **Generate Diagram** that
  drew the record under another key. Undo of that Generate puts them back. The
  `hash` value of such a color in **Specific rules** names the record key
  of the current view, so it changes with the view. The fix does not change
  Session files; a Session saved by 0.14.2 or earlier moves its edits at its
  first Generate.
- Label text and **Label visibility** typed in the Web feature popup and not
  yet applied stay when the labels are placed again while the popup is open.
- A Legend row's stroke in the Web app stays on the Result when a color rule
  edit gives the row another key (for example a caption edit of another rule),
  as **Generate Diagram** draws it. Before, the row showed the default line
  until the next Generate.
- A Legend row's stroke in the Web app no longer comes back on the features
  of the row's color after its row has left the Legend (for example the CDS
  row that becomes "other proteins" when a color rule takes a CDS feature under
  another caption). Before, a later Legend delete stroked those features again
  until the next **Generate Diagram**.
- After a palette switch or a Default colors change in the Web app, a Legend
  row whose name a color rule also uses takes its new color while the
  rule colors no drawn feature (for example **Use the "CDS" Legend color** on
  a CDS feature that a visibility rule hides), as **Generate Diagram** draws
  it. Before, the row kept the rule's color until the next Generate.
- A Legend row that a popup Legend name or a Legend rename adds in the Web
  app is drawn with the feature line, as **Generate Diagram** draws it, also
  after a Legend sort or move has put a GC content, GC skew, or Depth row
  first. Before, it had no line until the next Generate.
- A color rule edit in the Web app, made while a batch Result is shown on
  which the rule colors no feature, keeps that Result's Legend row of the
  rule's name and color (for example the CDS row after **Use the "CDS"
  Legend color** on a CDS of the other Result). The same holds on one Result when the rule
  colors only features that Feature visibility hides. Before, the row
  disappeared until the next **Generate Diagram**.
- A popup Legend name in the Web app on the only feature of a color
  rule's Legend row keeps the row in its place at **Generate Diagram**, as
  the live Legend shows it: when the new rule only gives that rule's
  condition a new name, it takes that rule's place in the color rule table.
  Before, Generate could draw the row ahead of rows it followed live.
- The Web app's header from 769 to about 970 px stays inside the window and
  keeps the Session title in view. Before, it ran past the right edge by up to
  32 px in French (and up to 8 px in English, Spanish, Portuguese, and
  Indonesian), cutting off the translation feedback link, and below about
  810 px it narrowed the Session title to nothing. Where the title still
  narrows, the closed Language menu now takes the width of the shown
  language's name (in browsers that support CSS `field-sizing`, such as
  Chrome and Edge); where the title would be narrower than 5rem, the header
  takes the two-row layout it has in a window 768 px wide or narrower. 0.14.2
  had both problems.
- The packaged palette `hearbal_garden` is spelled `herbal_garden`. The old
  name still works in `-p/--palette`, `FeatureOptions`, `ColorOptions`,
  Sessions, and requests, and draws the same colors; lists show only
  `herbal_garden`, and new requests name it.
- With no LOSAT or BLAST+ runtime, the error points to `gbdraw setup-losat` and
  no longer says that macOS and Windows have no LOSAT binary.

## Changed behavior

- **Circular labels that cannot fit.** Horizontal labels that need more than 1.5
  times the room around the map, when no single fix makes room, stop at once
  with the ways to make room, instead of after minutes of placing and searching
  (AP027078 `--labels out` at 8 pt: about 1.4 s instead of 28 s; MG1655
  `--labels out`: about 8 s instead of about 37 min). Labels that need 0.8
  times the room or more warn before they are placed (`LabelsMayNotFitWarning`);
  the Web app shows the warning under **Labels may not fit** while Generate
  runs, and after a stop the **Labels do not fit** dialog lists the ways to make
  room. Python raises `LabelsDoNotFitError` with no rows and a `capacity`
  (`LabelCapacityRecord`, exported from `gbdraw.api`) for each record that was
  not placed, and warns with `LabelsMayNotFitWarning`.
  See [Labels that do not fit](./REFERENCE/palettes-feature-rules-labels-shapes-and-tracks.md#labels-that-need-much-more-room-than-the-map-has).
- **Missing CairoSVG.** `gbdraw circular` and `gbdraw linear` exit with status 1,
  write no file, and print one error that names the formats and the install
  command when `-f` asks for PNG, PDF, EPS, or PS and CairoSVG is not
  installed. Before, they wrote only the SVG with a warning and exited 0.
- **LOSAT and BLAST+.** The LOSAT development build tracked in the source tree
  is removed. In a source checkout, run `gbdraw setup-losat` or put `losat` on
  `PATH`. The automatic order is the explicit option, conda, the managed
  install, `losat` on `PATH`, then NCBI BLAST+ on `PATH`. Each search logs one
  `INFO: Running ...` line that names the backend, its source, the executable,
  and its version. New Session runtime records no longer carry the `cli` field;
  older records still load.
- **Building from source.** Building gbdraw needs setuptools 77 or newer. The
  package lists Python 3.13.
- **Web header.** In the one-row header (769 px and wider), where every button
  already shows only its icon and the Session title would still be narrowed
  (for example from 769 to 1000 px in most languages), the Session block hides
  its "Session" caption and the toolbar gaps close to 2 px, so the Session
  title is wider (French at 900 px: its full 128 px, with the Language menu
  change under Fixed, instead of 28 px in 0.14.2).

## Added

- **Show Definition (Circular).** A **Show Definition** check in **Record
  Labels** (default on), `gbdraw circular --hide_definition`, and the config
  key `objects.definition.circular.show`. Off draws no center definition
  (name, strain, replicon, accession, length, GC) for any record, and the
  inside tracks keep no room for it; the plot title stays. A **Center Reserved
  Radius** you set still applies. Sessions store the setting from version 49;
  an older Session opens with the definition shown.
- The Web app leads to the Interactive Gallery's tutorials: a **Gallery &
  tutorials** header button, a **Learn with tutorials** card in the empty
  Preview on gbdraw.app, a **Tutorial** link on each example in **Load
  example**, and a link in the notice of a loaded example.
- `SECURITY.md` names the supported versions and how to report a vulnerability.
- The README and the Quickstart give two commands that draw a Circular map and a
  Linear comparison from the inputs in `gbdraw/web/tutorial-data/`. The README
  commands find that folder in the installed package, so they work after a
  `pip install` and in a source checkout.

## Faster

- **Circular labels.** Label placement and the label-fit search on large label
  sets do less repeated work and hold less memory, with the same SVG. Radial
  labels of MG1655 take about 8 s instead of 166 s, and of NC_002745 about 4 s
  instead of 25 s.
- **Web app redraws.** A Generate or live redraw that only moves or restyles
  what is drawn keeps the feature catalog the page holds instead of rebuilding
  it. On a Multi-record Circular canvas, a change to some records copies and
  re-reads only those records' parts. Circular labels no longer redraw when a
  definition font or the center radius changes, and Linear comparison ribbons
  no longer redraw when a definition or record label changes. The Worker
  parses each Depth file once and keeps it while **Show Depth** is off, so
  showing Depth again does not parse the file again. A Circular Species change
  that moves no track keeps the placed labels, and keeps the record definitions
  and the plot title where they do not draw the species (for example a typed
  plot title or the record-summary profile), and a Linear change that moves rows without changing a comparison keeps its
  ribbons. Turning the GC content or GC skew track on or off no longer redraws
  the feature contents, the Linear records, or the Circular ticks, and a Depth
  Legend label change keeps the drawn Depth paths. In a warm Generate of the
  vibrio Gallery Session, the prepared-input check takes 155 ms of CPU instead
  of 331 ms in CPython. A change of Separate strands, Resolve overlaps, or
  Label rendering rebuilds fewer render stages, and a warm Generate of a Linear
  similarity-groups or collinear figure with a reversed record is faster.

The drawn figures are the same as before, except where listed under
[Changed behavior](#changed-behavior).

## Installation and compatibility

Installation from PyPI and Bioconda is unchanged; building from source needs
setuptools 77 or newer, and a source checkout needs `gbdraw setup-losat` for
LOSAT (see [Changed behavior](#changed-behavior)). See
[Installation](./INSTALL.md). gbdraw 0.14.3
supports the same Python versions as 0.14.2 (3.10, 3.11, and 3.12; 3.13 is
listed and tested on `dev`). A Session saved by 0.14.2 opens in 0.14.3. gbdraw
0.14.3 saves Sessions as version 49, which gbdraw 0.14.0 to 0.14.2 refuse;
open them with 0.14.3 or later. See
[Session and request compatibility](./REFERENCE/session-and-request-compatibility.md).

## Known issues

- Edits made on a Linear File card that draws every record of its file, in a
  Session saved by gbdraw 0.13.0, do not follow its records after the first
  **Generate Diagram** splits the card into one card per record.

Carried from 0.14.1:

- A Session saved by 0.14.0 with a crop and ring LOSAT still holds unconverted
  rows. Its ring draws as in 0.14.0 on a forward crop, and at another
  unconverted place on a reverse-complemented crop, until the next **Generate
  Diagram**, which reuses the raw LOSAT cache and converts the rows.

Carried from 0.14.0:

- After a palette switch, or the **Legend item** choice on a palette Legend
  row, a color
  rule in the row's earlier color can show as one Legend row until the next
  **Generate Diagram**, which draws two.

See the [0.14.0 known limitations](./RELEASE_NOTES_0.14.0.md#known-limitations),
which still apply.

[Documentation home](./DOCS.md) | [0.14.2 release notes](./RELEASE_NOTES_0.14.2.md)
