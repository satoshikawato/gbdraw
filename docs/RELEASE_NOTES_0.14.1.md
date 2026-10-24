[Documentation home](./DOCS.md) | [Installation](./INSTALL.md) | [Changelog](../CHANGELOG.md) | [0.14.0 release notes](./RELEASE_NOTES_0.14.0.md)

# gbdraw 0.14.1 release notes

gbdraw `0.14.1` is a patch release of [0.14.0](./RELEASE_NOTES_0.14.0.md). It
fixes layout and Circular similarity ring problems, makes the web app faster
when you generate again, and gives each drawing in the web app its own
Undo/Redo History. Session files are unchanged: gbdraw 0.14.1 saves Sessions as
session version 48, with canonical `renderRequest` schema 9, and Sessions saved
by 0.14.0 open as before. [GitHub Releases](https://github.com/satoshikawato/gbdraw/releases)
gives the publication date, [PyPI](https://pypi.org/project/gbdraw/) lists the
published packages, and [Installation](./INSTALL.md) explains each install route.

## Fixed

### Circular similarity rings

- The command line and the Python API search the whole reference file for
  Circular similarity rings (`--losat losatn`, `--losat tlosatx`), as the web app
  does: every record with a sequence, uncropped, also records that are not drawn,
  and every FASTA entry in GFF3 mode. Ring E-values (and so the E-value filter)
  and LOSAT raw-cache keys change for runs that draw some records of a file, one
  record of it, or a crop. They now equal the web app's, so a Session saved from
  the command line and opened in the web app reuses its ring LOSAT results.
- Every ring table (an uploaded table, `--conservation_table`, the ring tables of
  a Session) is read in the search frame of the drawn record: the record after
  its crop, on its source strand, 1-based. On a reverse-complemented record the
  hits are drawn at their place; they were drawn mirrored. A row outside the
  drawn record is skipped with a warning.
- In the web app, ring LOSAT hits on a cropped record are drawn at their place.
  A Session saved by 0.14.0 with a crop and ring LOSAT keeps unconverted rows
  until the next **Generate Diagram**; see [Known issues](#known-issues).
- LOSAT, the Circular conservation database, and restored match sequences read
  the same records as the drawing when a GenBank file holds a record without a
  sequence, a file starts with a UTF-8 byte order mark, or a GFF3 `#N` selector
  names a record.

### Layout

- A Circular legend at a corner (`upper_left`, `upper_right`, `lower_left`,
  `lower_right`) on a Multi-Record Canvas, which the web app uses by default,
  was drawn over the records. It now clears them and the canvas grows.
- A top or bottom plot title is drawn outside a corner legend that was moved
  above or below the diagram; it was drawn through the legend rows.
- A one-gradient legend block (for example "Pairwise match identity") is as wide
  as its title, with the bar centred under it, in Linear legends and in top or
  bottom Circular legends. A title wider than the bar was cut off.
- A Linear plot title's box is the drawn width and height of its text, italic
  parts included. Linear figures with a plot title can be about 10 px shorter.
- A Circular corner legend that cannot sit beside the circle no longer moves
  feature labels.

- Adding a plot title to a single-record Circular figure no longer moves the GC
  content, GC skew, and tick rings. The layout reserves the center definition
  the figure draws without the title, so the figure keeps its untitled layout
  and the tick labels no longer run together. In 0.14.0 the shorter definition
  let the rings jump to the preset's lanes. Figures without a plot title are
  unchanged.
- A Circular custom track stack that adds an inside row the preset does not
  have, such as an AT skew row or an inside annotation slot, packs its rows
  from the feature ring toward the center again, as in 0.13.0.
  Its Auto widths follow the full definition the record draws alone. In 0.14.0
  such a stack moved to the preset's lanes whenever the center definition left
  room, so the record summary each record draws on a Multi-Record Canvas could
  move the GC content, skew, and tick rings and run the tick labels together.
  A record now draws these rings in the same place alone and on a Multi-Record
  Canvas. Figures from such a stack change. The default stack, a custom stack
  with only the preset's rows (comparison rings and depth series count as
  the preset's), and the Gallery figures are unchanged.
- A record drawn from another start (`--display_start_coordinate`, a records
  table start, or **Rotate record using this feature**) shows its length in the
  definition, for example `16,569 bp`, in Circular and Linear. It showed the
  wrapped source spans, such as `[3,307..16,569], [1..3,306] bp`.

### Web app

- A legend or plot title that you moved stays inside the canvas: the canvas
  grows to hold it, and **Reset Layout** returns the automatic canvas.
- A Circular **Inner gap** or **Outer gap** typed as `1px` is stored as `1`, so
  Save and a fresh Load agree.
- The compact Result Preview layout and the width of **Reset zoom** apply in
  every UI language; they applied only in English.
- The Run info tab, the feature popup's label note, the Label Not Shown dialog,
  and the Choose records dialog show the UI language in every language other than English.
  The Run info text that Generate builds stays English.
- Panning the Result from blank space and releasing over a feature or a
  pairwise match no longer opens that popup.
- A Legend offset that nothing names is released with its Result, and a reflow
  lands in its own Result's drawing.

## Documentation

- The HmmtDNA_ATskew, tobacco-chloroplast, and Vnig_TUMSAT-TG-2018 Gallery
  tutorials start from the Web defaults and generate after each change, with a
  preview of every stage.
- The Vibrio harveyi group collinear tutorial says that the Linear ruler keeps
  each record's own coordinates.

## Added

- gbdraw warns when the center definition text of a Circular record touches
  what an inner track draws (GC content or GC skew, a feature block, or the
  axis circle), and offers a smaller definition font that fits. The command line
  prints a `WARNING:` with the fitting `--definition_font_size` and exits 0;
  Python raises `LabelLeadersCrossWarning`; the web app's label-fit notice
  offers it under **Choose a fix**. The figure is drawn exactly as before.

## Faster web app

The web app no longer repeats work on inputs that have not changed:

- **Generate Diagram** encodes a newly added file in milliseconds instead of
  about 0.5 s per 12 MB. **Save Session** no longer re-reads or re-hashes the
  input files, and returning to Linear no longer reads every File card again.
- **Align** and the feature-edit helpers reuse the inputs the last Generate
  parsed. On the hepatoplasmataceae Gallery Session the Align helper takes about
  0.5 s instead of about 10 s.
- Circular similarity rings read each ring file and the reference sequence once
  per file, not on every Generate.
- A blastp threshold change filters the kept hits instead of parsing every raw
  protein table, and the protein extraction reads the GenBank records the Worker
  already holds.
- The Generate after a failed one (for example after the Label fit choice)
  reuses the sources the failed Generate parsed: about 2.9 s instead of 4.0 s on
  a 5.5 Mb GenBank file.
- LOSAT hashes each record's FASTA once, and the first blastp Generate after an
  upload sends and parses each new file once: about 4.4 to 5.1 s instead of 7.6
  to 8.0 s for six cards of a 17.5 MB GenBank file.
- Switching between Circular and Linear no longer re-sends or re-parses the
  inputs; over four Generates of the `vibrio-harveyi-group-collinear` Gallery
  Session the bytes sent to the Worker fall by 42 %.

The drawn figures are the same as before.

## Changed behavior

### Undo and Redo per drawing

- A mode switch (Circular or Linear) shows the other drawing, writes nothing,
  and records no step. Undo and Redo walk the shown drawing's steps.
- Each drawing keeps 30 Undo steps, and a project stream keeps 30 more, under
  one shared memory budget. 0.14.0 kept 30 in all.
- Reset settings, and edits to the Session title, Auto Reflow, PNG DPI,
  Instant Preview, the LOSAT execution, or the input drafts, are project steps
  that either drawing can undo. Reset settings is one step: one Undo restores
  both drawings. The dialog is unchanged.
- Reset settings also gives the multi-record layout a new drawing's defaults,
  derives the Circular comparison rings again from the kept files, and resets a
  hidden drawing's applied palette and alignment Reset receipt with it.
- Instant Preview **On** applies every drawing's queued palette in one step.
- A drawing switch waits while a color rule table edit is being matched.
- In **Canvas padding**, an empty or partial entry keeps the side until you leave
  the field, then becomes 0, and Save writes a number.
- The five auxiliary file cards (Priority File (TSV), Whitelist File, Blacklist
  File, Default colors file, Specific rules file) are shared by both drawings as
  in 0.14.0. **Remove** takes the file's rows out of both drawings and one Undo
  brings them back. Because such a change is a project step, an Undo pressed
  while the other drawing is shown also undoes it.

### Circular rings from the command line

Ring E-values and LOSAT raw-cache keys of command-line and Python runs change as
described under [Fixed](#circular-similarity-rings). The ring rows stay in the
search frame of the drawn record, so `--losat_output_dir`, `conservation.tsv`,
and Sessions hold the same frame as before.

## Installation and compatibility

Installation is unchanged; see [Installation](./INSTALL.md). gbdraw 0.14.1
supports the same Python versions as 0.14.0 (3.10, 3.11, and 3.12). A Session
saved by 0.14.0 opens in 0.14.1. See
[Session and request compatibility](./REFERENCE/session-and-request-compatibility.md).

## Known issues

New in this release:

- A Session saved by 0.14.0 with a crop and ring LOSAT still holds unconverted
  rows. Its ring draws as in 0.14.0 on a forward crop, and at another
  unconverted place on a reverse-complemented crop, until the next **Generate
  Diagram**, which reuses the raw LOSAT cache and converts the rows.

Carried from 0.14.0:

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
  on cropped or reverse-complemented records until you click **Generate Diagram**.
  Colors are unaffected. An SVG downloaded before **Generate Diagram** carries the
  old IDs.
- A label whitelist row that the former **Whitelist only** action saved before 0.14.0,
  for a feature without a unique qualifier, names the feature by its old drawn ID.
  On a cropped or reverse-complemented record it no longer matches that feature,
  and its label is no longer whitelisted. To show the label, set **Label visibility**
  to **On** in the feature popup, or add a matching row under **Labels** →
  **LABEL FILTERING** → **Whitelist**.

See the [0.14.0 known limitations](./RELEASE_NOTES_0.14.0.md#known-limitations),
which still apply.

[Documentation home](./DOCS.md) | [0.14.0 release notes](./RELEASE_NOTES_0.14.0.md)
