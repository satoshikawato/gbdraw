[Documentation home](./DOCS.md) | [Installation](./INSTALL.md) | [Changelog](../CHANGELOG.md) | [0.14.1 release notes](./RELEASE_NOTES_0.14.1.md) | [0.14.0 release notes](./RELEASE_NOTES_0.14.0.md)

# gbdraw 0.14.2 release notes

gbdraw `0.14.2` is a patch release of [0.14.1](./RELEASE_NOTES_0.14.1.md). It
fixes the "This feature only" fill, draws the Auto organism name and the
Circular centre definition by one rule, lays out labels and the Legend at the
size their text is drawn, and makes the web app, the command line, and SVG files
faster or smaller. Session files are unchanged: gbdraw 0.14.2 saves Sessions as
session version 48, with canonical `renderRequest` schema 9, and Sessions saved
by 0.14.1 open as before. [GitHub Releases](https://github.com/satoshikawato/gbdraw/releases)
gives the publication date, [PyPI](https://pypi.org/project/gbdraw/) lists the
published packages, and [Installation](./INSTALL.md) explains each install route.

## Fixed

- A "This feature only" color colors only the clicked feature of the clicked
  File card and record, live and after Generate, also when two File cards hold
  the same file or two records share an ID. The rule's `hash` value now names
  the feature instance (`<feature id>@<record key>`). Resetting one copy leaves
  the others as they are.
- Two identical `source` features of a record, one shown by its own Feature
  visibility, have the same feature IDs in the Web feature list as in the
  drawing and in a feature identity table.
- The Legend is laid out at the width of its painted captions, so a caption in
  a monospace or other wide face is no longer cut off at the canvas edge. A
  corner Legend that fits beside the records no longer grows the canvas.
- A Bakta `.` organism no longer draws an italic "."; a Circular record with no
  organism and no typed title draws no plot title; a typed plot title with
  nothing to draw no longer fails.
- Group color edits keep the colors of the loaded features, and a palette row
  restored after a palette switch comes back in the color the Legend shows.
- The Editor's "N features ready" shows the feature count of the mode on screen.
- A Result saved or exported with the web app in a language other than English
  no longer keeps the preview's button attributes on comparison bands.
- A registered font with a corrupt table fails every text measurement the same
  way, instead of measuring without that table after the first failure.

## Changed behavior

- **Auto organism name.** The Circular centre, the Circular plot title, and the
  Web Linear Auto definition draw line 1 as the scientific name, italicized by
  the nomenclature code of the record's GenBank lineage, and line 2 as the
  strain, drawn once even when `/organism` already holds it. Saved Linear Auto
  definitions keep their text until the file is read again.
- **Circular centre definition.** Its lines are never wrapped. The inside
  tracks first move outward at their widths, and narrow only when every
  narrowed track stays at or above 30 px (`--track_narrowing_min_px`);
  `--no_track_narrowing` turns narrowing off. A width you set is never
  narrowed for the definition. The `DEFINITION_RESERVED` failure is gone. The
  Web **Labels do not fit** dialog offers **Edit Species / Strain…**, and its
  Cancel reads **Keep as drawn** over a drawn figure.
- **Text boxes.** Feature labels, tick labels, the Linear ruler and scale bar,
  the Circular centre definition, and the Legend are laid out at the box their
  text is drawn in. Canvases change by a few pixels, and stacked Circular
  labels sit slightly closer. Saved Results keep their drawing until Generate
  redraws them.
- **Circular radius.** Before Generate, a Circular track row's **Radius** says
  "Auto" and its **Width** estimate is the Auto size the Python layout starts
  from. The ring warning above **Generate** and in the **Pairwise Comparisons**
  card is gone; Generate still reports a ring stack that does not fit.
- **Group colors.** A group color edit writes one exact rule (feature type plus
  product, gene, or label text) where the renderer confirms it colors only the
  edited features. Coloring the 2761 "hypothetical protein" CDS of NC_010162
  writes 1 rule instead of 2761. The rule also colors matching features of
  files added later.
- **SVG files.** Every path number has at most 3 decimals, and Depth, GC
  content, and GC skew paths drop points that do not move the outline. Files
  are smaller (Circular AP027078 with Depth: 5.6 MB to 0.48 MB). A comparison
  match carries its ID only in `data-gbdraw-match-id`; scripts that read
  `data-gbdraw-pairwise-match-id` from new files must read the new name. Older
  files that carry only the old name still load and still go through
  `enrich_svg`.

## Faster

- **Web app.** The first Generate starts faster (parallel wheel downloads,
  Pyodide bytecode, and a prefetch of the boot files after you choose a file).
  Cancel keeps the diagram Worker on a cross-origin-isolated page such as
  gbdraw.app, so the next Generate does not boot Pyodide again. Loading a large
  Session is about twice as fast, and a live edit no longer serializes the whole
  SVG. A returning visitor reloads gbdraw.app without downloading Pyodide, the
  wheels, fonts, icons, or LOSAT wasm again.
- **Command line and Python.** Replaying a Session, rendering, and Similarity
  group inference do less repeated work and hold less memory, with the same
  SVG and the same Sessions; GC content and GC skew of a 52 Mb Circular diagram peak at 576 MiB RSS instead of 750 MiB. The installed `gbdraw` command sets
  `OPENBLAS_NUM_THREADS=1` before NumPy loads, which cuts the CPU time of a
  one-shot render without changing its wall time.

The drawn figures are the same as before, except where listed under
[Changed behavior](#changed-behavior).

## Documentation

- The aminoglycoside BGC Gallery tutorial shows the antiSMASH `gene_kind`
  qualifier in a step of its own before the color step.

## Installation and compatibility

Installation is unchanged; see [Installation](./INSTALL.md). gbdraw 0.14.2
supports the same Python versions as 0.14.1 (3.10, 3.11, and 3.12). A Session
saved by 0.14.1 opens in 0.14.2. See
[Session and request compatibility](./REFERENCE/session-and-request-compatibility.md).
The browser wheel grows from 3.7 MB to 6.5 MB because it carries Pyodide
bytecode, so the PyPI wheel and the sdist each grow by about 2.8 MB.

## Known issues

New in this release:

- A "This feature only" color belongs to its record key, like the other per-feature edits, so it does not apply while the record is drawn under another key (OV-1949).

Carried from 0.14.1:

- A Session saved by 0.14.0 with a crop and ring LOSAT still holds unconverted
  rows. Its ring draws as in 0.14.0 on a forward crop, and at another
  unconverted place on a reverse-complemented crop, until the next **Generate
  Diagram**, which reuses the raw LOSAT cache and converts the rows.

Carried from 0.14.0:

- After a palette switch, or **Apply to all** on a palette Legend row, a color
  rule in the row's earlier color can show as one Legend row until the next
  **Generate Diagram**, which draws two.

See the [0.14.0 known limitations](./RELEASE_NOTES_0.14.0.md#known-limitations),
which still apply.

[Documentation home](./DOCS.md) | [0.14.1 release notes](./RELEASE_NOTES_0.14.1.md)
