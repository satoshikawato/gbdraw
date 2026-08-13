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
source annotation model. A `translation` attribute is used when present.
Otherwise a valid CDS may be translated from sequence, strand, phase or
`codon_start`, and genetic code.

When one source contains several records, select the intended record by ID or
index, or explicitly expand all records. Use Circular presentation for a
complete record whose biological topology is circular. Cropping a region or
splitting a sequence does not make it circular.

## Comparison and numeric tables

BLAST-compatible input uses the 12 outfmt 6 columns, with optional outfmt 7
comment lines: `qseqid`, `sseqid`, `pident`, `length`, `mismatch`, `gapopen`,
`qstart`, `qend`, `sstart`, `send`, `evalue`, and `bitscore`. Query and subject
direction must match the displayed endpoint mapping.

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
| Records | One of `gbk`, or both `gff` and `fasta` | `record_label`, `record_subtitle`, `record_id`, `region`, `reverse_complement`, `order`, `row`, `column` |
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

`priorities` is a comma-separated qualifier list. `value` and `keyword` fields
are case-insensitive regular expressions. Selector qualifiers may also use the
documented synthetic keys `location`, `record_location`, and `hash` where that
surface supports exact feature identity. Canonical request schema 6 additionally
reserves `__gbdraw_instance_hash__`: its value is a case-sensitive literal with
the form `fi1_` followed by 26 lower-case base32 characters. The Web editor
generates this opaque value; do not derive or edit it by hand. It is stable only
for the canonical `recordKey` domain that produced it, so a downloaded exact-row
TSV must be replayed with the same record keys and a gbdraw version that supports
schema 6. The biological qualifier `instance_hash` remains an ordinary regex
qualifier. Schema 6 also reserves `__gbdraw_semantic_scope__` for Web-authored,
future-matching Feature-group rules. Its opaque `fs1:<kind>:<encoded-value>`
literal can represent a Feature type, underlying legend caption, rendered label,
source annotation label, or Similarity group. Keep this machine-generated row
unchanged; unlike exact instance rows, it intentionally applies to matching
Features introduced by regeneration. Under schema 5, the same spelling remains
an ordinary biological regex qualifier. Visibility `action` is `show`, `off`,
or `exclude_matching`.

Table precedence and the meaning of those actions are documented in [Feature
presentation](palettes-feature-rules-labels-shapes-and-tracks.md#feature-presentation).
