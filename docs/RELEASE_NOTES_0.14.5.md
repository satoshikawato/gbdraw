[Documentation home](./DOCS.md) | [Installation](./INSTALL.md) | [Changelog](../CHANGELOG.md) | [0.14.4 release notes](./RELEASE_NOTES_0.14.4.md) | [0.14.3 release notes](./RELEASE_NOTES_0.14.3.md)

# gbdraw 0.14.5 release notes

gbdraw `0.14.5` is a patch release of [0.14.4](./RELEASE_NOTES_0.14.4.md). The
Web app now draws a LOSATP link or Collinear block that uses a CDS across a
circular record's origin over the CDS's own arc, as the command line and the
Python API have done since 0.14.4. The release also adds `--wrap_length`, which
draws one long Linear record over several rows (command line and Python), four
built-in color tables, and reads more real-world GenBank, EMBL, and GFF3 files.
Sessions are now saved as session version 50 (0.14.4 writes 49; the canonical
`renderRequest` schema stays 9): **a Session saved by 0.14.5 opens only in
0.14.5 or later.** Sessions saved by 0.14.4 and earlier open as before.
[GitHub Releases](https://github.com/satoshikawato/gbdraw/releases) gives the
publication date, [PyPI](https://pypi.org/project/gbdraw/) lists the published
packages, and [Installation](./INSTALL.md) explains each install route.

<!-- PENDING: lead sentence of the summary above may need a mention of the LOSAT Wasm swap, OV-2369, OV-2370/OV-2380 and OV-2371 once they merge. -->

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
  - A Geneious LOCUS line with spaces in the name (`LOCUS       Contig 1 edited
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
<!-- PENDING: Fixed entry for OV-2369 (title and user-visible impact from its PR once merged). -->
<!-- PENDING: Fixed entry for OV-2370 / OV-2380 (title and user-visible impact from its PR once merged). -->
<!-- PENDING: Fixed entry for OV-2371, #1532: LOSATP Pairwise on a Linear row with two or more records no longer exits 1 and searches every adjacent-row record pair (command line and Python only); the ordered `-b/--blast` refusal now names `--comparisons_table`. Write it from the PR body if it has merged. -->

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
<!-- PENDING: Changed behavior for the LOSAT Wasm swap: engine commit (LOSAT dev dde737f18), and that a raw LOSAT search cached or saved by 0.14.4 runs again once. Write from its PR body. -->

## Added

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

<!-- PENDING: add or drop known issues for the LOSAT Wasm swap, OV-2369, OV-2370/OV-2380 and OV-2371 once they merge. -->

- A Similarity group member popup and the Similarity alignment panel still show
  the coordinates of a CDS across the origin as `1..` the record length. The
  drawing is correct; only this text is wrong.
- After **Reset Settings** on one diagram mode, the other mode's Result keeps its
  hidden labels and edited label text, and after a Session load the Legend rows of
  its removed color rules, until the next **Generate Diagram**. Since 0.14.3 at
  least.
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
