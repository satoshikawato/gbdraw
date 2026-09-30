[Documentation home](../DOCS.md) | [Tutorials](../TUTORIALS/README.md) | [Technical documentation](README.md) | [Command line](command-line.md) | [FAQ](../FAQ.md)

# Input formats and TSV schemas

## Sequence and annotation files

GenBank and GBFF files may contain one or more biological records. DDBJ records
are accepted when downloaded in GenBank flat-file format; native EMBL
flat-file syntax is not. gbdraw preserves the parsed record ID, description,
topology, sequence, feature locations, strand, and qualifiers. The reader does
not split one biological sequence into artificial records.

GFF3 must be paired with FASTA from the same biological source. GFF3 column 1
and the first token of the matching FASTA header must agree exactly.
Coordinates are 1-based and inclusive; strand is `+` or `-`; CDS phase is `0`,
`1`, or `2`. `ID` values should be unique, and `Parent` should preserve the
source annotation model.

CDS protein sequences use the GenBank `/translation` qualifier or the GFF3
`translation` attribute when present. Otherwise gbdraw translates the CDS
nucleotides with `transl_table` (default 1). The reading frame starts at
`codon_start` or, for GFF3, at the phase of the 5'-most CDS part. As in INSDC
`/translation`, the first residue is `M` when the 5' end is complete, the frame
starts at the first base, and the first codon is a start codon of that table;
for example, `GTG` and `TTG` become `M` in table 11. A 5' end is incomplete when
its position is fuzzy (`<` on the plus strand, `>` on the minus strand) or GFF3
has `start_range` (plus strand) or `end_range` (minus strand). Such a CDS, one
read from frame 2 or 3, and a `pseudo` or `pseudogene` CDS translate the first
codon literally. `transl_except` is not applied. Feature popups and Interactive SVG
metadata do not translate a `pseudo`, `pseudogene`, or fuzzy-location CDS.

When one source contains several records, select the intended record by ID or
index, or explicitly expand all records. Use Circular presentation for a
complete record whose biological topology is circular. Cropping a region or
splitting a sequence does not make it circular.

## Comparison and numeric tables

BLAST-compatible input is tab-separated UTF-8 text. gbdraw reads the first 12
outfmt 6 columns by position: `qseqid`, `sseqid`, `pident`, `length`,
`mismatch`, `gapopen`, `qstart`, `qend`, `sstart`, `send`, `evalue`, and
`bitscore`. Extra columns, such as those from `-outfmt "6 std qlen slen"`, are
ignored and reported in an INFO log. A line whose first non-blank character is
`#` is an outfmt 7 comment; a `#` or quote inside a field is part of the value.
An empty table or a table with only comment lines has no rows. The IDs must be
non-empty, `length`, `mismatch`, `gapopen`, and the four coordinates must be
integers, and the other values must be finite numbers. The command line, Python
API, web app, and Circular similarity rings apply the same rule. For a Linear
comparison, a row with fewer than 12 columns, a value of the wrong type, or a
missing or unreadable file stops the run; the error names the line of a
malformed row. A Circular similarity ring whose table is rejected is skipped
with a warning and keeps its ring position.

Query and subject direction must match the displayed endpoint mapping. For a
Linear comparison, a row whose `qseqid` or `sseqid` names the other endpoint or
another displayed record stops the run; in the web app, Generate reports a
comparison-endpoint error and keeps the previous Result. IDs are compared with
the record ID and name. A version suffix difference, such as `NC_000913` and `NC_000913.3`, is
accepted. IDs that match no displayed record keep the positional assignment;
the command line logs a warning. The SVG record-ID metadata always names the
endpoint records.

Depth input has `reference_name`, a 1-based positive `position`, and a
non-negative `depth`. Files are normally headerless. One header line is
accepted when the position or depth fields in the first row are nonnumeric.
Each file is one measured series for the named record; a missing series is not
equivalent to zero depth.

## Manifest tables

Manifest files are UTF-8 TSV with real tab characters. A UTF-8 BOM is accepted.
Unknown columns are rejected. Relative file paths resolve from the directory
containing the table.

| Table | Required columns | Optional columns |
|---|---|---|
| Records | One of `gbk`, or both `gff` and `fasta` | `record_label`, `record_subtitle`, `record_id`, `region`, `reverse_complement`, `topology`, `display_start`, `order`, `row`, `column` |
| Linear comparisons | `blast`, `query`, `subject` | None |
| Circular conservation | `blast` | `label`, `color`, `comparison_fasta` |
| Circular tracks | `id`, `renderer` | `side`, `r`, `w`, `inner_gap_px`, `outer_gap_px`, `z`, `params` |
| Annotations | `set_id`, `id`, `mark` | Target and presentation fields listed below |

One records-table row represents one displayed record. A table uses either
GenBank rows or GFF3/FASTA rows; it cannot mix the two forms. `record_id`
selects from a multi-record source. A row-scoped `region` contains coordinates
only and applies after selection. `reverse_complement` is a row-scoped boolean.

`order`, `row`, and `column` are positive integers. Explicit `order` values
sort before blank values; equal values retain table order. When placement is
present, every row needs a `row`. `column` controls left-to-right order, and
duplicate row/column cells are rejected. Use table placement instead of
repeated surface-specific position options.

Linear comparison `query` and `subject` values accept a displayed `#index` or
unique record ID. The endpoints must be different. A comparisons table cannot
be combined with the positional `--blast` form.

Circular track-table row order is slot order. A row with `side=axis` must use
the `features` renderer and establishes the track-axis boundary. Structural
values belong in their columns; `params` contains renderer settings such as
`nt`, `set_id`, or `legend_label`.

## Annotation table fields

Each annotation row defines exactly one coordinate target or one feature
target. `mark` accepts `line`, `bracket`, `band`, or `highlight`.

- A coordinate target uses `record`, `start`, `end`, `coordinate_space`,
  `wraps_origin`, and `out_of_bounds`.
- A feature target uses `record`, `feature_selector`, `envelope`, and
  `circular_path`.
- Presentation fields are `label`, `lane`, `legend_label`, `stroke`,
  `stroke_width`, `stroke_dasharray`, `line_cap`, `fill`, `fill_opacity`,
  `hatch_angle`, `hatch_spacing`, `hatch_color`, `hatch_width`, `hatch_cross`,
  `label_color`, `label_font_size`, `label_orientation`, `label_position`, and
  `label_offset`.

Coordinate targets are 1-based and inclusive. `wraps_origin=true` is
Circular-only; split an origin-crossing Linear range into two rows. A feature
selector uses qualifier expressions, with multiple conditions separated by
semicolons. Explicit `lane` values are zero-based. The renderer selects rows by
`set_id`, so one table can hold independently placed annotation sets.

In the Web app, **Region Annotations → Download TSV** saves the current editor
draft as `annotations.tsv`, including edits made since the last Generate.
It works offline and is disabled until the draft contains an annotation row.
The file can be loaded with **Import TSV**, Python's `read_annotation_table()`,
or the CLI's `--annotation_table` option.

The download preserves effective row targets and styles, including explicit
no-fill, unique record IDs, and one-based `#N` record bindings. A blank `fill`
cell in a styled row means no fill; omitting the column retains the Web import
default. TSV does not retain empty sets, metadata, or the distinction between
an inherited set style and a row override. Tabs and line breaks within cells
are replaced with spaces.

## Styling tables

Most styling tables are headerless TSV. Label-override and visibility readers
also accept their documented header row.

| Table | Columns in order |
|---|---|
| Default colors | `feature_type`, `color` |
| Specific colors | `feature_type`, `qualifier_key`, `value`, `color`, optional `caption` |
| Qualifier priority | `feature_type`, `priorities` |
| Label whitelist or blacklist | `feature_type`, `qualifier`, `keyword` |
| Label overrides | `record_id`, `feature_type`, `qualifier`, `value`, `label_text` |
| Feature visibility | `record_id`, `feature_type`, `qualifier`, `value`, `action` |

A Specific-colors `color` is `none` (no fill, any case), an SVG color name,
`#RGB`, or `#RRGGBB`. Other values, including hex colors with alpha, are
rejected with their line number. The Web app converts a color name to hex when
it reads the table. In this table, cells such as `None`, `NA`, and `null` are
values, not blanks.

`priorities` is a comma-separated qualifier list. Pattern `value` and `keyword`
fields use case-insensitive Python regular expressions. Specific-color and
Label-override patterns accept Python-only syntax such as `(?i)NADH`,
`(?P<enzyme>NADH)` and `NADH\Z`; Unicode matching follows Python semantics.
The Web app prepares Color rules and Label TSV patterns with that same Python
owner before committing live changes. Invalid syntax is rejected even when the
feature catalog is empty or contains no matches. A table-structure error or
runtime preparation failure is distinct from invalid regex syntax.

Feature Search and search inside downloaded Interactive SVG use
case-insensitive JavaScript regex instead. Python and JavaScript patterns are
not translated between these surfaces. Existing Color pattern fields can hold
an unapplied display draft; it is not part of the accepted TSV rule or saved
Session. See [Color and Label patterns](web-app.md#color-and-label-patterns) for
Retry, Revert and Save/Generate/Export behavior.

Selector qualifiers may also use the documented synthetic keys `location`,
`record_location`, and `hash` where that surface supports exact feature identity.
Visibility `action` is `show`, `off`,
or `exclude_matching`.

Table precedence and the meaning of those actions are documented in [Feature
presentation](palettes-feature-rules-labels-shapes-and-tracks.md#feature-presentation).

## Display and feature placement tables

Records-table `topology` accepts `auto`, `circular`, or `linear`; blank is auto.
`display_start` is an optional 1-based source base in `1..L`.
Blank means no override/no additional shift. Explicit start and crop cannot be
combined. Direct `--record_topology` and `--display_start_coordinate` flags target one
resolved record; use a records table for distinct per-record settings.

`--feature_placement_table` accepts UTF-8 TSV (including BOM):

| Column | Meaning |
|---|---|
| `record` | Optional unique record ID or displayed `#index`; omission must resolve to one record |
| `feature_selector` | Required exact selector, such as `protein_id=NP_054479.1`, `hash=<biologicalFeatureId>`, or a qualifier value |
| `placement` | Required `auto`, `main`, `outward`, `inward`, `above`, or `below` |
| `level` | Omit for Auto/Main; directional targets accept blank (lane 1) or `1` |

Selectors must match exactly one original-source feature. A gene qualifier may
match both a gene annotation and its CDS; a unique protein ID avoids that
ambiguity. Unknown/stale identities, duplicate resolved targets, extra columns,
wrong-mode sides and unsupported resolved layouts are errors. GFF duplicate
identities use complete original-source order, including features hidden by
loading or visibility rules. Changing visibility does not renumber them.

Auto removes an override. Source-known cropped-out, hidden or underlay features
retain dormant intent and reserve no foreground lane; restoring visibility or
the crop reactivates it. Persist exact record/biological-feature identities,
never an SVG fragment ID or a transient lane number.
