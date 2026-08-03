[Documentation home](../DOCS.md) | [Input how-to](../HOW_TO/GUI/use-genbank-and-gff3-fasta-inputs.md) | [Fixture provenance](tutorial-fixture-provenance.md)

# Input formats and TSV schemas

## Sequence and annotation files

GenBank, GBFF, EMBL, and DDBJ flat files may contain one or more complete biological records. gbdraw preserves record ID, description, topology, sequence, feature locations, strand, and qualifiers supported by the parser.

GFF3 must be paired with the FASTA sequence from the same biological source. GFF3 column 1 and the first token of the matching FASTA header must agree exactly. Coordinates are 1-based and inclusive; strand is `+` or `-`; CDS phase is `0`, `1`, or `2`. `ID` values should be unique, and `Parent` should preserve the source annotation model. A `translation` attribute is used when present; otherwise a valid CDS may be translated from sequence, strand, phase/codon start, and genetic code.

Tutorial fixtures do not manufacture contigs by cutting a complete sequence. The Lambda GFF3 + FASTA fixture represents the whole `NC_001416.1` record. Multi-record examples use genuine biological records or replicons.

## Comparison and numeric tables

BLAST-compatible comparison input uses outfmt 6 columns, with optional outfmt 7 comment lines: `qseqid`, `sseqid`, `pident`, `length`, `mismatch`, `gapopen`, `qstart`, `qend`, `sstart`, `send`, `evalue`, and `bitscore`.

Depth input contains `reference_name`, 1-based positive `position`, and non-negative `depth`. Files are normally headerless; one header line is accepted when the position or depth fields in the first row are nonnumeric.

## Manifest tables

Tables are UTF-8 TSV with real tab characters. A UTF-8 BOM is accepted. Unknown columns are rejected. Relative file paths resolve from the table file's directory.

| Table | Required columns | Optional columns |
|---|---|---|
| Records | One of `gbk`, or both `gff` and `fasta` | `record_label`, `record_subtitle`, `record_id`, `region`, `reverse_complement`, `order`, `row`, `column` |
| Linear comparisons | `blast`, `query`, `subject` | None |
| Circular conservation | `blast` | `label`, `color`, `comparison_fasta` |
| Circular tracks | `id`, `renderer` | `side`, `r`, `w`, `inner_gap_px`, `outer_gap_px`, `z`, `params` |
| Annotations | `set_id`, `id`, `mark` | Target and style fields listed below |

A records table cannot mix GenBank rows with GFF3/FASTA rows. `order`, `row`, and `column` are positive integers. A complete grid gives every row a grid row; duplicate row/column cells are rejected. A row-scoped `region` contains coordinates only and must not repeat a record selector.

Comparison `query` and `subject` accept a displayed `#index` or unique record ID. They must identify different displayed records. A Circular track row with `side=axis` must use the feature renderer and is the track-axis boundary.

## Annotation table fields

Each annotation row defines exactly one coordinate target or one feature target. `mark` accepts `line`, `bracket`, `band`, or `highlight`.

- Coordinate target: `record`, `start`, `end`, `coordinate_space`, `wraps_origin`, `out_of_bounds`.
- Feature target: `record`, `feature_selector`, `envelope`, `circular_path`.
- Presentation: `label`, `lane`, `legend_label`, `stroke`, `stroke_width`, `stroke_dasharray`, `line_cap`, `fill`, `fill_opacity`, `hatch_angle`, `hatch_spacing`, `hatch_color`, `hatch_width`, `hatch_cross`, `label_color`, `label_font_size`, `label_orientation`, `label_position`, `label_offset`.

Coordinate targets use 1-based inclusive positions. Linear annotations cannot wrap the origin. A feature selector uses supported qualifier expressions; several conditions are separated by semicolons. `lane` is zero-based when explicitly set.

## Styling tables

Default-color, specific-color, qualifier-priority, label-filter/override, and feature-visibility tables have separate schemas and precedence. Their token inventory is in [Palettes, feature rules, labels, shapes, and track renderers](palettes-feature-rules-labels-shapes-and-tracks.md). Use the public tutorial fixture copies rather than paths under `tests/test_inputs`.
