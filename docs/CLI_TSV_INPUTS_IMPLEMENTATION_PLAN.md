[Home](./DOCS.md) | [CLI Reference](./CLI_Reference.md)

# CLI TSV Inputs Implementation Plan

Date: 2026-07-07

This plan adds three table-based CLI entry points while keeping the existing
options stable:

- `--records_table`
- `--conservation_table`
- `--circular_track_table`

`--definition_line_style` is intentionally out of scope for this change.

## Goals

- Replace fragile parallel CLI lists with row-based TSV input where the values
  are naturally coupled.
- Keep short, simple commands short. Existing options such as `--gbk`,
  `--gff`, `--fasta`, `--conservation_blast`, and `--circular_track_slot`
  remain supported.
- Avoid overloading existing option names with two meanings. A BLAST TSV file
  and a metadata TSV file should not share one argument name.
- Reuse existing parser and renderer behavior where possible, especially for
  track slot specs and record ordering.
- Make session save/load and the web app preserve the new table-backed inputs.

## New Arguments

### `--conservation_table`

Chosen name: `--conservation_table`

Rejected alternative: `--conservation_blast`

Reason: BLAST outfmt 6/7 is already tabular. If `--conservation_blast` accepts
both raw BLAST output and a metadata table, the CLI has to guess from file
content or extension. A separate name keeps the interface obvious.

Mode: circular only.

Mutual exclusion:

- Cannot be combined with `--conservation_blast`.
- Cannot be combined with `--conservation_labels`.
- Cannot be combined with `--conservation_colors`.

Initial TSV schema:

```tsv
blast	label	color
ref_a.blast.tsv	Reference A	#E15759
ref_b.blast.tsv	Reference B	#4E79A7
```

Columns:

- `blast`: required. Path to a BLAST outfmt 6/7 file.
- `label`: optional. Conservation ring label.
- `color`: optional. SVG color name or hex color.

Path behavior:

- Relative paths resolve against the directory containing the table.
- Absolute paths are preserved.

Runtime mapping:

- `blast` becomes `conservation_blast_files`.
- `label` becomes `conservation_labels`.
- `color` becomes `conservation_colors`.
- Row order defines ring order.

Validation:

- Empty table: error.
- Missing `blast` column: error.
- Empty `blast` cell: error with row number.
- Unknown columns: error for the first implementation, so typos fail loudly.
- Existing conservation thresholds remain normal CLI options:
  `--conservation_reference`, `--evalue`, `--bitscore`, `--identity`,
  `--alignment_length`, `--conservation_ring_width`,
  `--conservation_ring_gap`.

### `--circular_track_table`

Chosen name: `--circular_track_table`

Rejected alternative: `--circular_track_slot`

Reason: `--circular_track_slot` currently means one inline slot spec. Reusing it
for a path would make the parser decide whether a string is a slot spec or a
file path. A separate argument avoids that ambiguity.

Mode: circular only.

Mutual exclusion:

- Cannot be combined with `--circular_track_slot`.
- Cannot be combined with `--circular_track_order`.
- Cannot be combined with `--circular_track_axis_index`.

Initial TSV schema:

```tsv
id	renderer	side	r	w	params
features	features	axis			
gc_content	dinucleotide_content	inside		0.1	nt=GC
gc_skew	dinucleotide_skew	inside		0.1	nt=GC
ticks	ticks	inside			tick_label_layout=label_in_tick_out
```

Columns:

- `id`: required. Slot ID.
- `renderer`: required. Slot renderer name.
- `side`: optional. One of `outside`, `axis`, `inside`.
- `r`: optional. Radius scalar, same meaning as inline slot `r`.
- `w`: optional. Width scalar, same meaning as inline slot `w`.
- `spacing`: optional compatibility alias, same meaning as inline slot
  `spacing`.
- `inner_gap_px`: optional.
- `outer_gap_px`: optional.
- `z`: optional.
- `params`: optional comma-separated renderer params, using the existing
  `key=value,key=value` parser.

Important design choice:

- The table should not expose `--circular_track_axis_index` directly.
- Instead, `side=outside|axis|inside` expresses user intent per row.
- `side=axis` is a table-only convenience value. It is not passed through as a
  `CircularTrackSlot.side` value, because the existing slot model accepts only
  `outside`, `inside`, and `overlay`.
- The loader computes the slot order and the axis boundary index from `side`.

Side mapping:

- Rows with `side=outside` are placed before the axis boundary.
- Rows with `side=axis` are placed on the axis boundary. Only one axis row is
  allowed initially, and it must be a `features` renderer row.
- The `side=axis` features row is emitted as a split feature slot:
  `side=overlay` with `lane_direction=split`, while the computed
  `circular_track_axis_index` is set to that row's slot index.
- Rows with `side=inside` are placed after the axis boundary.
- If `side` is omitted, default to `inside`, except the first `features`
  renderer row may default to `axis` only if no explicit axis row exists.
- Relative row order is preserved within each group: outside rows, the optional
  axis row, then inside rows.

Runtime mapping:

- Build canonical inline slot specs from rows.
- Reuse `parse_circular_track_slots()` and
  `normalize_circular_track_slots_with_axis()`.
- Pass the generated slots and generated axis index into the existing circular
  diagram path.

Validation:

- Empty table: error.
- Missing `id` or `renderer`: error.
- Duplicate `id`: handled by the existing slot parser or rejected before it.
- More than one `side=axis`: error.
- `side=axis` on a non-`features` renderer: error.
- A `side=axis` row with conflicting `side`, `lane_direction`, or `lanes`
  renderer params: error.
- Invalid `side`: error.
- Unknown columns: error.
- Ambiguous geometry combinations continue to be validated by the existing
  track slot parser.

### `--records_table`

Chosen name: `--records_table`

Rejected alternatives:

- `--records_list`: acceptable, but the file is a typed table rather than a
  simple list.
- `--gbk_list`: too narrow because GFF3 plus FASTA must be supported.

Modes: circular and linear.

First implementation rule:

- `--records_table` is an alternative input source.
- It cannot be combined with `--gbk`, `--gff`, or `--fasta`.
- When `--records_table` is used, row-coupled legacy options are also
  disallowed. Put those values in the table instead.
  - Linear: `--record_label`, `--record_subtitle`, `--record_id`,
    `--reverse_complement`, and `--region`.
  - Circular: `--multi_record_position`.

This keeps precedence simple. The table still provides the same capability as
`--gbk` plus repeated `--record_label`, but it does so as one coherent input
manifest.

GenBank TSV schema:

```tsv
gbk	record_label	record_subtitle	record_id	region	reverse_complement	order	row	column
a.gbk	Strain A		#1		0	1	1	1
b.gbk	Strain B		#1		0	2	1	2
c.gbk	Strain C		#1		0	3	2	1
```

GFF3 plus FASTA TSV schema:

```tsv
gff	fasta	record_label	record_subtitle	record_id	region	reverse_complement	order	row	column
a.gff3	a.fna	Strain A		chr1	1000-9000	0	1	1	1
b.gff3	b.fna	Strain B		chr1		0	2	1	2
```

Columns:

- `gbk`: GenBank path. Required for GenBank rows.
- `gff`: GFF3 path. Required for GFF3/FASTA rows.
- `fasta`: FASTA path. Required for GFF3/FASTA rows.
- `record_label`: optional. Linear definition top line.
- `record_subtitle`: optional. Linear definition subtitle line.
- `record_id`: optional. Record selector for the row.
- `region`: optional. Region crop spec for the row.
- `reverse_complement`: optional. Reverse-complement flag for the row.
- `order`: optional. Positive integer used to sort rows before loading.
- `row`: optional. Circular multi-record row.
- `column`: optional. Circular within-row ordering hint.

Row model:

- One TSV row represents one displayed record.
- If an input file contains exactly one record, `record_id` may be omitted.
- If an input file contains more than one record, the row must provide
  `record_id` so the loader selects exactly one displayed record.
- To display multiple records from the same input file, repeat the file path in
  multiple rows and give each row its own `record_id`.
- A selector that matches no records or multiple records is an error.
- `record_label`, `record_subtitle`, `region`, `reverse_complement`, `order`,
  `row`, and `column` apply to the row's selected displayed record only.

Input kind rule:

- A table must be either all GenBank rows or all GFF3/FASTA rows.
- Mixed `gbk` and `gff`/`fasta` rows are rejected initially.

Path behavior:

- Relative `gbk`, `gff`, and `fasta` paths resolve against the table directory.
- Absolute paths are preserved.

Linear runtime mapping:

- `gbk` rows become `args.gbk`.
- `gff` and `fasta` rows become `args.gff` and `args.fasta`.
- `record_label` becomes repeated `--record_label` values.
- `record_subtitle` becomes repeated `--record_subtitle` values.
- `record_id` becomes repeated `--record_id` values.
- `reverse_complement` becomes repeated `--reverse_complement` values.
- `region` becomes repeated `--region` values. Selectorless cells are row-scoped
  and become `#<loaded-index>:<region>` after row ordering and record selection.

Circular runtime mapping:

- `gbk` rows become `args.gbk`.
- `gff` and `fasta` rows become `args.gff` and `args.fasta`.
- `record_id` selects one displayed record per row. This requires extending the
  circular loader path to honor the same row-level selectors used by linear
  mode.
- `reverse_complement` applies to the row's selected displayed record.
- `region` applies to the row's selected displayed record. Selectorless cells
  are row-scoped and become `#<loaded-index>:<region>` after row ordering and
  record selection.
- `record_label` and `record_subtitle` remain linear display metadata in the
  first implementation. Circular mode accepts the columns in shared manifests
  but ignores their values.
- `row`, `column`, and `order` can generate `multi_record_positions` when
  `--multi_record_canvas` is set.

## How `row`, `column`, and `order` Absorb `--multi_record_position`

Current CLI:

```text
--multi_record_position "#1@1" --multi_record_position "#2@1" --multi_record_position "#3@2"
```

Equivalent table:

```tsv
gbk	record_label	row	column	order
a.gbk	A	1	1	1
b.gbk	B	1	2	2
c.gbk	C	2	1	3
```

Interpretation:

- `order` sorts the input rows before loading. If omitted, table row order is
  preserved.
- `row` selects the multi-record canvas row.
- `column` controls order within that row.
- The first implementation can translate this to existing tokens like
  `#1@1`, `#2@1`, `#3@2`.

Limit:

- Empty grid cells are not preserved at first. For example, columns `1` and `3`
  without column `2` are compressed into two records in the same row. This
  matches current `multi_record_positions` behavior.
- True fixed-column placement would require extending the circular multi-canvas
  layout data model beyond `selector@row`.

Validation:

- `order`, `row`, and `column` must be positive integers when present.
- `row` and `column` are meaningful only in circular mode.
- If any `row` is present, all rows should have a `row` value.
- If any `column` is present, all rows in rows with explicit placement should
  have a `column` value.
- Duplicate `(row, column)` pairs are errors.
- When `row`/`column` are present without `--multi_record_canvas`, warn and
  ignore, matching current `--multi_record_position` behavior.

## Parser Architecture

Add a small shared parser module, for example:

```text
gbdraw/io/cli_tables.py
```

Suggested dataclasses:

```python
ConservationTableRow
CircularTrackTableResult
RecordsTableResult
```

Suggested functions:

```python
read_conservation_table(path: str) -> ConservationTable
read_circular_track_table(path: str) -> CircularTrackTable
read_records_table(path: str) -> RecordsTable
```

Use the standard library `csv.DictReader(..., delimiter="\t")` unless an
existing pandas-based helper already makes this cleaner. The table reader should
centralize:

- UTF-8 text loading.
- Header validation.
- Blank line skipping.
- Row number reporting.
- Relative path resolution.
- Positive integer parsing.
- Boolean parsing for `reverse_complement`.
- Path dependency collection for session capture. The parser should expose the
  table file plus every referenced `gbk`, `gff`, `fasta`, and `blast` path so
  session code can embed the complete manifest.

Avoid adding a broad table framework. Three focused readers are easier to audit.

## Shared Input Validation

Update the common input validation so `--records_table` is a first-class input
source.

Rules:

- Exactly one of these input source families must be provided:
  - `--records_table`
  - `--gbk`
  - `--gff` plus `--fasta`
- `--records_table` cannot be combined with `--gbk`, `--gff`, or `--fasta`.
- Error wording should name all valid choices, for example:
  `Either --records_table, --gbk, or both --gff and --fasta must be provided.`
- The common validator should not require `--gbk` or `--gff`/`--fasta` when
  `--records_table` is present.

## CLI Integration Points

### `gbdraw/circular.py`

Argument additions:

- Add `--records_table`.
- Add `--conservation_table`.
- Add `--circular_track_table`.

Validation in `_get_args()`:

- `--records_table` cannot combine with `--gbk`, `--gff`, or `--fasta`.
- `--records_table` cannot combine with `--multi_record_position`; use the
  `row` and `column` table columns instead.
- `--conservation_table` cannot combine with `--conservation_blast`,
  `--conservation_labels`, or `--conservation_colors`.
- `--circular_track_table` cannot combine with `--circular_track_order`,
  `--circular_track_slot`, or `--circular_track_axis_index`.

Runtime in `run_circular_from_namespace()`:

- Read `--records_table` before loading records and populate local input lists,
  selectors, reverse flags, row-scoped regions, and placement metadata.
- Extend the circular `load_gbks()` / `load_gff_fasta()` path, or add a small
  table-specific loader wrapper, so each records-table row loads exactly one
  displayed record. Do not let one row silently expand into multiple rendered
  records.
- Read `--conservation_table` before assigning `conservation_blast_files`,
  `conservation_labels`, and `conservation_colors`.
- Read `--circular_track_table` before validating legacy geometry conflicts and
  populate `circular_track_slots_or_none` plus `circular_track_axis_index`.
- Generate `multi_record_positions` from records-table placement fields when
  applicable.

### `gbdraw/linear.py`

Argument additions:

- Add `--records_table`.

Validation in `_get_args()`:

- `--records_table` cannot combine with `--gbk`, `--gff`, or `--fasta`.
- `--records_table` cannot combine with `--record_label`,
  `--record_subtitle`, `--record_id`, `--reverse_complement`, or `--region`;
  use the table columns instead.

Runtime in `run_linear_from_namespace()`:

- Read `--records_table` before setting `file_count`.
- Populate `args.gbk` or `args.gff`/`args.fasta` equivalent local lists.
- Populate record labels, subtitles, record selectors, reverse flags, and
  regions from the table.
- Enforce the one-row-one-displayed-record rule after loading. If a row's input
  file would load more than one record without a selector, error and ask the
  user to add `record_id`.
- Keep existing explicit CLI handling unchanged when no table is supplied.

## Session and Web App Plan

### Session sidecars

Update session capture so table-backed runs can round-trip.

Recommended behavior:

- Store the original table files as file bindings.
- Also store every file referenced by those tables as file bindings:
  - `--records_table`: `gbk`, `gff`, and `fasta` paths.
  - `--conservation_table`: `blast` paths.
  - `--circular_track_table`: currently no file-reference columns.
- Preserve `cliInvocation.args` with `--records_table`,
  `--conservation_table`, and `--circular_track_table` when those were used.
- When restoring an older session, continue emitting legacy arguments.
- When restoring a new table-backed session, materialize referenced files first,
  write restored TSV copies whose path cells point to those materialized files,
  and pass the restored table paths to the new table arguments.
- Prefer paths relative to the restored TSV location when rewriting table cells,
  so the restored manifest remains readable and movable inside the temporary
  session directory.

Likely files:

- `gbdraw/session_io.py`
- `gbdraw/cli_utils/session.py`
- `tests/test_session_io.py`
- `tests/test_refresh_gallery_sessions.py`

### Web UI

The web app can continue using structured state internally. It does not need to
switch its UI model to TSV immediately.

Recommended first step:

- Add wheel capability detection for the new arguments.
- When the user uploads multiple conservation BLAST files with labels/colors,
  the web app may continue emitting legacy CLI args initially.
- For saved sessions that contain table-backed CLI invocations, materialize the
  table file, materialize every referenced input file, rewrite the restored
  table paths, and pass the restored table argument.

Recommended second step:

- Use generated in-memory TSVs in Pyodide for conservation series and circular
  track slots only if doing so simplifies argument generation.
- Avoid adding visible TSV import/export UI until the CLI behavior is stable.

Likely files:

- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/services/history-snapshot.js`
- `gbdraw/web/js/state.js`
- `tests/web/*.test.mjs`

## Documentation Plan

Update:

- `docs/CLI_Reference.md`
- `docs/TUTORIALS/2_Comparative_Genomics.md`
- `docs/TUTORIALS/3_Advanced_Customization.md`
- `docs/RECIPES.md`

Add examples for:

- Circular conservation table.
- Circular track table.
- Linear records table with labels.
- Circular records table with multi-record canvas placement.

Keep documentation wording clear that legacy options remain supported.

## Test Plan

### Unit tests

Add parser-focused tests:

- `tests/test_cli_tables.py`

Cover:

- Valid conservation table.
- Relative path resolution.
- Missing required columns.
- Empty required cells.
- Unknown columns.
- Invalid positive integers.
- Invalid booleans.
- Mixed `gbk` and `gff`/`fasta` rows.
- Duplicate circular track `side=axis`.
- Circular track `side=axis` on non-`features` renderer.
- Duplicate records-table `(row, column)`.
- Records-table path dependency collection for session capture.

### Circular CLI tests

Extend existing circular tests or add focused tests:

- `--conservation_table` forwards blast files, labels, and colors in row order.
- `--conservation_table` rejects legacy conservation list options.
- `--circular_track_table` forwards slots and computed axis index.
- `--circular_track_table` rejects slot/order/axis legacy combinations.
- `--records_table` loads GBK paths.
- `--records_table` applies `record_id`, `region`, and `reverse_complement`
  to the selected displayed record.
- `--records_table` rejects a row that would expand into multiple displayed
  records without `record_id`.
- `--records_table` with `row`/`column` generates expected
  `multi_record_positions`.
- `--records_table` rejects explicit `--multi_record_position`.

Likely files:

- `tests/test_circular_conservation.py`
- `tests/test_circular_track_slots.py`
- `tests/test_circular_multi_canvas.py`

### Linear CLI tests

Add tests for:

- `--records_table` supplies GBK inputs and `record_label`.
- `--records_table` supplies GFF3/FASTA inputs.
- `record_id`, `region`, and `reverse_complement` map to existing loader
  behavior.
- `--records_table` rejects explicit `--gbk`/`--gff`/`--fasta` combinations.
- `--records_table` rejects explicit row-coupled legacy options:
  `--record_label`, `--record_subtitle`, `--record_id`,
  `--reverse_complement`, and `--region`.
- `--records_table` rejects a row that would expand into multiple displayed
  records without `record_id`.

Likely files:

- `tests/test_linear_selectors.py`
- `tests/test_linear_definition_alignment.py`

### Session tests

Add tests for:

- Table file bindings are captured.
- Files referenced by table cells are captured.
- Session restore materializes referenced files, rewrites restored table paths,
  and emits the table argument.
- Older legacy sessions still restore legacy options.

Likely files:

- `tests/test_session_io.py`
- `tests/test_refresh_gallery_sessions.py`

### Web tests

Only after web integration:

- Capability detection recognizes the new options.
- Config/session restore preserves table-backed invocations.
- Existing non-table UI flows still emit valid legacy CLI args.

Run:

```bash
pytest tests/ -v -m "not slow"
ruff check gbdraw/ --select=E,F,W --ignore=E501,W503
```

For web changes, also check Node or Python Playwright availability as described
in `AGENTS.md`.

## Rollout Order

1. Add `gbdraw/io/cli_tables.py` and parser unit tests, including path
   dependency collection.
2. Add circular `--conservation_table` support and tests.
3. Add circular `--circular_track_table` support and tests.
4. Add `--records_table` support to linear and circular, including shared input
   validation and circular row-level selector/reverse/region support, without
   session/web changes yet.
5. Add session sidecar round-trip support for table files and every file
   referenced by table cells.
6. Update web capability detection and restore behavior.
7. Update documentation and recipes.
8. Run focused tests, then fast suite.

This order keeps each step reviewable and prevents the table parser, rendering
behavior, session IO, and web UI from changing all at once.

## Compatibility Policy

- Existing options remain supported.
- New table options are additive.
- New table options should have clear mutual exclusions where precedence would
  otherwise be ambiguous.
- Error messages should name the table path, row number, column name, and bad
  value whenever possible.
- Reference SVGs should not change unless row ordering or generated slot order is
  intentionally exercised by a new test.

## Open Decisions

- Whether `--records_table` should eventually support mixed GenBank and
  GFF3/FASTA rows. Recommendation: no for the first implementation.
- Whether the non-table circular CLI should eventually gain `--record_id`,
  `--region`, and `--reverse_complement`. Recommendation: keep that separate;
  `--records_table` should support these row-level fields from the first
  implementation.
- Whether `column` should preserve empty grid cells. Recommendation: defer; map
  `column` to within-row ordering first.
- Whether to add visible TSV import/export controls to the web UI. Recommendation:
  defer until the CLI format settles.
