[Documentation home](./DOCS.md) | [Installation](./INSTALL.md) | [Changelog](../CHANGELOG.md) | [0.14.3 release notes](./RELEASE_NOTES_0.14.3.md) | [0.14.2 release notes](./RELEASE_NOTES_0.14.2.md)

# gbdraw 0.14.4 release notes

gbdraw `0.14.4` is a patch release of [0.14.3](./RELEASE_NOTES_0.14.3.md). It
fixes LOSATP comparisons of circular records that have a CDS across the
origin, which were drawn over the whole record. It also draws partial feature
ends open (no arrowhead, no edge stroke), adds similarity-ring color schemes to
the command line and the Python API, and applies more Web app settings to the
shown Result as soon as you commit a value. LOSATP Collinear blocks are now
built as the longest chains of anchors and have one set of defaults on the
command line, in the Python API, and in the Web app. Sessions are still saved as session
version 49 (the canonical `renderRequest` schema stays 9): Sessions saved by
0.14.3 open as before.
[GitHub Releases](https://github.com/satoshikawato/gbdraw/releases) gives the
publication date, [PyPI](https://pypi.org/project/gbdraw/) lists the published
packages, and [Installation](./INSTALL.md) explains each install route.

## Fixed

- **LOSATP links and Collinear blocks across the origin.** A CDS whose location
  crosses a circular record's origin, such as
  `complement(join(3717939..3718359,1..656))`, was read as an interval over the
  whole record. Every LOSATP Pairwise link and Collinear block that used it was
  drawn over the whole record, in every version since 0.11.0. 0.14.4 draws
  them over the CDS's own arc, as one piece on each side of the origin, both
  pieces carrying the same block and protein metadata; a match popup on one
  piece shows that piece's coordinates, and its FASTA download uses that piece's
  span. The blocks that are called do not change; only their drawn coordinates
  do. **Redraw** figures that show LOSATP comparisons of circular records with a
  CDS across the origin, which is typical of complete bacterial genomes from
  RefSeq. Records without such a CDS draw the same bytes as before.
- **Collinear blocks as longest chains.** A run of collinear genes stays one
  block where it used to split into several blocks or single links: where a
  query gene has a second anchor, where a single link continues the run, and
  where a compatible cluster sits within **Max unit gap** but not next to it.
  Chaining keeps the same anchors, and every
  anchor is still drawn; an anchor that no chain takes is a single-link block.
  Collinear figures change, in the command line, the Python API, and the Web app.
  **Redraw** the figures you keep. Saved derived LOSATP results are rebuilt once
  on their next use. `--collinear_max_conflicts_in_merge_gap` keeps its meaning for `0` and
  `1`: it limits the other anchors strictly inside the box between two
  neighbouring anchors of a chain.
- **Collinear speed at a large Max unit gap.** Block calling no longer slows down
  as **Max unit gap** (`--collinear_max_unit_gap`) grows. 2000 sparse anchors took
  seconds at a gap of 1000 and take tens of milliseconds at any gap.
- **Similarity groups.** A Linear layout with two or more records in one row
  (`--multi_record_position`, `LinearMultiRecordOptions(multi_record_positions=...)`,
  or a request with grid rows) no longer fails after the whole search with
  "Ordered BLAST/protein comparison inputs are ambiguous when a Linear row
  contains multiple records". The links join every record pair in adjacent rows,
  as Collinear blocks already do, and a Session saved on such a layout stores the
  drawn pairs and replays them. **Max hits per protein for groups**
  (`--losatp_member_max_hits`) no longer counts a protein's
  hit to itself. The bitscore normalization of a within-record search no longer
  takes its scale from those self hits; a within-record table with fewer than 12
  other hits uses the square-root-of-lengths fallback that a small cross-record
  table already used. Group memberships can change when a member limit is set and
  when a record's within-record search has in-paralog hits, most in small genomes,
  phages, plasmids, and organelles. Layouts with one record per row, and runs
  without a member limit and without in-paralog hits, are unchanged.
- In the Web app, **Generate Diagram** after the drawn records of a Linear Session
  whose records were aligned or moved changed (a File removed, replaced, or added,
  or a record turned off) failed with an unnamed error (`UNKNOWN`) instead of
  drawing; since 0.14.3 at least. It now draws, each kept record at the place the
  Session gave it, and an added record at its default place. The Similarity
  alignment plan is still cleared by such a change; **Align** again to recompute.
- **Similarity alignment center.** A Similarity alignment anchored on such a CDS
  placed the record's middle at the shared center. It now places the CDS
  itself, in the command line, the Python API, and the Web app alike.
- `import gbdraw` works on Python 3.14 and later. Python 3.14 removed
  `importlib.abc.Traversable`, which gbdraw's configuration loader imported;
  it now reads `Traversable` from `importlib.resources.abc`.
- The Circular legend lists only the similarity rings a slot draws; a ring of a
  disabled `sequence_conservation` slot no longer appears.
- In the Web app, a Legend name that joins another Legend row of the same color
  keeps that row where it was after **Generate Diagram** as well, unless moving
  the renamed rules would recolor another feature. A feature given a new Legend
  name apart from the rest of its row shows its new row at once where Generate
  draws it, instead of last.
- In the Web app, after a Legend name edit, a feature color change, or a feature
  color reset that adds, joins, splits off, or moves a color rule's Legend row,
  the Legend shows the rows at once in the order **Generate Diagram** draws
  them, including a row of another feature type and an old row the edit moves,
  instead of first showing them elsewhere and then jumping. Rows placed this way
  do not count as a hand-edited Legend order.
- In the Web app, on a Result loaded from a Session saved without a render
  request (Session versions 27 to 30), settings no longer say **Live edit**
  while nothing applies until Generate.
- In the Web app, a color picked without a pointer first (a picker opened from
  the keyboard) is one **Undo** step and applies as its setting does. Before,
  Undo skipped it.
- In the Web app, changing a feature's color for its color rule, when another
  rule of the same Legend name and color also draws that row, shows the row
  where **Generate Diagram** draws it at once. Since 0.14.0 the row stayed in
  its old place until the next Generate.
- In the Web app, a color rule edit is applied when another rule's Legend row
  shows a Legend color of its own, for example a feature type's Legend color
  that a **Legend name** given in the feature popup kept on the new row. 0.14.3
  refused every such edit, and a Legend item color change showed no message.
- In the Web app, **Save Session** refuses, before it downloads anything, a
  Session that would expand beyond the 512 MiB that **Load Session** opens, and
  a Session whose gzip file is over the 200 MiB that **Load Session** opens. The
  message names the Session's size and its largest inputs.

## Changed behavior

- **Circular legend for similarity rings.** The legend lists the drawn rings by
  name from `outer` to `inner` above one `identity` band, with a swatch per
  ring when the ring colors differ. It replaces the `Conservation identity` row
  and the per-ring gradient rows.
- **Live settings in the Web app.** These settings rerender the shown Result as
  soon as a value is committed, and the Result equals a **Generate Diagram**
  from the same draft: **Species**, **Default font size**, the Linear definition
  row **Size**, **Weight**, and **Color**; record names (Linear **Organism /
  strain**, Circular **Record label** and **Subtitle**, and the record popup's
  name fields); the **Dinucleotide content/skew** **Window**; Region Annotation
  **Fill color** and **Line color**; Comparison **Match height**; **Center
  Reserved Radius**; and a Linear custom-stack row's **Height**. **Strain**,
  **Keep Full Definition**, and the Linear subtitle still apply on Generate.
  Where a value cannot apply yet (an annotation item with a Legend label, an
  item added since the last Generate, a features or Depth row, a record whose
  File was replaced), its note says **Applies on Generate**.
  **Show GC Content**, **Show GC Skew**, **Show Depth**, **Hide GC Content**,
  **Hide GC Skew**, the on/off switch and removal of a GC content, GC skew, or
  Depth row of the custom stack, and **Use custom stack** redraw the shown
  Result at once (a track that other stack rows reach asks first, see
  below), as Generate would, also for a track the Result does not draw
  yet; the track comes back as the last Generate built it. Undo of a turned-off
  switch that also turned off an Overlay row, and of **Hide GC Content** or
  **Hide GC Skew** in a Circular custom stack or with Circular track widths or
  radii, applies at once and returns the rows it turned off; Redo hides them
  again, and a Save after the Undo holds them. An Overlay row turned back on by
  hand while its track is back returns at once with the values of the last
  Generate, and a row that left when **Use custom stack** changed returns when
  it changes back. These wait for Generate: a track the Result did not draw when its
  Session was loaded, a Depth file change, a row added since the last
  Generate, **Hide GC Content** or **Hide GC Skew** turned off by hand in a
  Circular custom stack or with Circular track widths or radii, and **Use
  custom stack** after a Session load or while Region Annotation sets wait for
  Generate.
- **Collinear defaults on every surface.** Similarity-group inference is off by
  default on the command line and in the Python API too
  (`--collinear_infer_orthogroups off`), and the Web Collinear **Max target
  seqs** and **Max hits per protein for groups** are unbounded by default, as on
  the command line. A Session, request, or recorded command that names Collinear
  without the inference setting still replays with inference on. The anchor mode
  `one_to_one` and the merge orientation `either` are gone: a saved `either` reads
  as `strand`, which it drew as; a saved `one_to_one` reads as `rbh`, and the
  Web app lists that change at Load (the command line and Python log a warning).
  Without similarity groups, Collinear searches only the displayed record pairs
  on every surface, whatever **Search between** or `--collinear_search_scope`
  holds; the drawn blocks are the same. **Search between** shows only while
  **Build similarity groups with within-record searches** is on, and **Max hits
  per protein for groups** only while it is on or **Anchor mode** is All hits.
  The block popup row reads **Anchors**, and block rows built without similarity
  groups no longer name a group scope. The Collinear option help text and the
  Collinear tutorial now describe what each setting changes.
- **Turning off a track that other stack rows reach.** Turning off a track that
  another custom stack row starts at or extends to (an Overlay row's anchor or
  **Extend to**) with the row's switch (Circular and Linear), or with **Hide GC
  Content** or **Hide GC Skew** while **Use custom stack** is on (Circular),
  asks first. The dialog names those rows and offers to turn the track and those
  rows off as one **Undo** step, or to cancel. Before, the track hid alone into
  a stack that **Generate Diagram** refuses.
- **Use custom stack after a Session load.** On a Circular Result drawn with a
  custom stack or track widths or radii, **Use custom stack** set the other way
  than the Session saved it applies on Generate, with every track it would hide
  or show. Before, the tracks it hides left at once and the stack stayed drawn.
- **Smaller Session files.** **Save Session** stores an input file that the
  Session uses more than once (for example one Depth file shown as two series)
  only once; only files with equal bytes, kind, and original file name are
  merged. Such a Session file is smaller. Sessions are still saved as version
  49 and the file format does not change.
- **Partial feature ends.** A feature end that its location marks as partial
  is drawn open in Linear and Circular diagrams: no arrowhead and no stroke on
  that edge, the way a record's display-start seam is already drawn. This covers
  ends written in the file (`<1..1635`, `48045..>48169`,
  `complement(<1..242)`) and ends cut by `--region` or a records-table
  `region`. A feature whose 5' end is partial keeps its arrowhead and has no
  stroke on its tail edge. Only the outer ends of a multi-part location count.
  The figures of records with such features change, in the command line, the
  Python API, and the Web app (**Generate Diagram** and Session replay draw
  them the same way). A region crop now marks every cut end as partial: a
  minus-strand feature cut by the region kept no `<`/`>` marker, and any crop
  dropped the markers written in the file. Because of this, a CDS that a crop
  made partial no longer shows its translation in the Web feature popup, and
  a minus-strand CDS cut at its 5' end by a crop no longer gets its first
  residue replaced with M in LOSATP. A block outline with one open end joins
  its corners at the closed end.

## Added

- **Every Collinear option on the command line.** `--collinear_anchor_mode
  {all,rbh}` and `--collinear_merge_orientation {strand,order}`, the two settings
  that only the Web app had, and the Python field
  `LinearComparisonOptions.collinearity_infer_orthogroups` (default `False`).
  `--collinear_unit_mode`, `--collinear_max_conflicts_in_merge_gap`, and
  `--collinear_max_paralog_links_per_orthogroup` are listed in
  `gbdraw linear --help`.
- **Similarity-ring color schemes.** `--conservation_color_scheme` (`same`, the
  default; `ten`; `per_five`), the config key
  `objects.conservation.color_scheme`, and `ComparisonRingOptions.color_scheme`
  color the rings without a color of their own by drawn position from the
  outermost ring. `--conservation_colors` accepts `auto` for such a ring, and a
  Python ring track with `color=None` no longer requires every ring to leave its
  color unset. The default draws the same one-color rings as before. The Web app
  setting follows in a later release.

## Removed

- The **Block e-value** row of the Collinear match popup, the
  `collinearity_block_evalue` column of the block table, and the
  `data-collinearity-block-evalue` SVG attribute. The value has been empty since
  0.14.0.

The drawn figures are the same as before, except where listed above.

## Installation and compatibility

Installation from PyPI and Bioconda is unchanged; see
[Installation](./INSTALL.md). gbdraw 0.14.4 supports Python 3.10 to 3.14
(0.14.3 does not import on 3.14). A Session saved by 0.14.3 opens in 0.14.4,
and 0.14.4 saves Sessions as version 49, as 0.14.3 does. See
[Session and request compatibility](./REFERENCE/session-and-request-compatibility.md).

## Known issues

- After **Reset Settings** on one diagram mode, the other mode's Result keeps its
  hidden labels and edited label text until the next **Generate Diagram**. Since
  0.14.3 at least.
- A Similarity group member popup and the Similarity alignment panel still show
  the coordinates of a CDS across the origin as `1..` the record length. The
  drawing is correct; only this text is wrong.
- Edits made on a Linear File card that draws every record of its file, in a
  Session saved by gbdraw 0.13.0, do not follow its records after the first
  **Generate Diagram** splits the card into one card per record.

- After you choose to turn off a track that other rows reach in the dialog that
  asks first, an Overlay row you turn back on by hand, while its track is back,
  returns in the Result at once, while the note beside the track stack still says
  **Applies on Generate**. The Result then equals the next Generate.

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

[Documentation home](./DOCS.md) | [0.14.3 release notes](./RELEASE_NOTES_0.14.3.md)
