[Documentation home](./DOCS.md) | [Installation](./INSTALL.md) | [Changelog](../CHANGELOG.md) | [0.14.3 release notes](./RELEASE_NOTES_0.14.3.md) | [0.14.2 release notes](./RELEASE_NOTES_0.14.2.md)

# gbdraw 0.14.4 release notes

gbdraw `0.14.4` is a patch release of [0.14.3](./RELEASE_NOTES_0.14.3.md). It
fixes LOSATP comparisons of circular records that have a CDS across the
origin, which were drawn over the whole record. It also adds similarity-ring
color schemes to the command line and the Python API, and applies more Web app
settings to the shown Result as soon as you commit a value. Sessions are still
saved as session version 49 (the canonical `renderRequest` schema stays 9):
Sessions saved by 0.14.3 open as before.
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
  draws it, instead of last, when every named color rule is of one feature
  type.
- In the Web app, a color picked without a pointer first (a picker opened from
  the keyboard) is one **Undo** step and applies as its setting does. Before,
  Undo skipped it.

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

## Added

- **Similarity-ring color schemes.** `--conservation_color_scheme` (`same`, the
  default; `ten`; `per_five`), the config key
  `objects.conservation.color_scheme`, and `ComparisonRingOptions.color_scheme`
  color the rings without a color of their own by drawn position from the
  outermost ring. `--conservation_colors` accepts `auto` for such a ring, and a
  Python ring track with `color=None` no longer requires every ring to leave its
  color unset. The default draws the same one-color rings as before. The Web app
  setting follows in a later release.

The drawn figures are the same as before, except where listed above.

## Installation and compatibility

Installation from PyPI and Bioconda is unchanged; see
[Installation](./INSTALL.md). gbdraw 0.14.4 supports Python 3.10 to 3.14
(0.14.3 does not import on 3.14). A Session saved by 0.14.3 opens in 0.14.4,
and 0.14.4 saves Sessions as version 49, as 0.14.3 does. See
[Session and request compatibility](./REFERENCE/session-and-request-compatibility.md).

## Known issues

- A Similarity group member popup and the Similarity alignment panel still show
  the coordinates of a CDS across the origin as `1..` the record length. The
  drawing is correct; only this text is wrong.
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
  row, a color rule in the row's earlier color can show as one Legend row until
  the next **Generate Diagram**, which draws two.

See the [0.14.0 known limitations](./RELEASE_NOTES_0.14.0.md#known-limitations),
which still apply.

[Documentation home](./DOCS.md) | [0.14.3 release notes](./RELEASE_NOTES_0.14.3.md)
