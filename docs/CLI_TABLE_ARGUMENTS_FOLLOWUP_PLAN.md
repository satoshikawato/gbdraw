# CLI Table Arguments Follow-up Implementation Plan

Status as of 2026-06-07: planning document.

The first-pass `--input_table`, `--depth_track_table`, and `--track_table`
implementation is in place. This follow-up plan covers the remaining work that
was intentionally deferred:

- Add `--blast_table`.
- Make the web UI emit `--input_table` and `--track_table` from structured
  state, matching the existing web depth-track table staging.
- Keep docs, examples, and broader verification aligned after those changes.

This is not a new manifest design. The same table rules from
`docs/CLI_TABLE_ARGUMENTS_PLAN.md` apply: headered TSV files, table-relative
paths, explicit row-numbered validation, selector resolution through
`DisplayRecordContext`, and narrow CLI-boundary adapters.

## Goals

1. Replace positional BLAST comparison files with named table rows.
2. Keep web-generated CLI arguments readable and debuggable by staging TSV
   files instead of serializing compact slot strings.
3. Preserve stable user-facing IDs (`input_id`, `track_id`, and future
   comparison row IDs) until the renderer truly needs an ordered index.
4. Keep backward compatibility for existing CLI shortcuts and web settings.
5. Add focused tests first, then broaden verification enough to catch web and
   comparison regressions.

## Non-goals

- Do not introduce a project manifest or session file.
- Do not redesign linear comparison rendering in the first `--blast_table`
  pass. Existing linear comparison groups are still indexed by the gap between
  adjacent displayed records.
- Do not add per-row BLAST filtering columns initially. Existing global
  `--evalue`, `--bitscore`, `--identity`, `--alignment_length`, and
  `--pairwise_match_style` remain global defaults.
- Do not make circular conservation use `--blast_table` in the first pass unless
  the linear path is already stable. Circular conservation has separate
  `--conservation_blast` semantics and reference-side rules.

## Phase 1: `--blast_table`

### User Model

Users describe each BLAST outfmt 6/7 file as a row attached to two displayed
records. This removes positional guessing from `-b file1 file2 ...` and lets
tables join naturally with `--input_table` IDs.

Example:

```tsv
query_id	subject_id	file	order	enabled
input:mje	input:mela	MjeNMV.MelaMJNV.tblastx.out	1	true
input:mela	input:pemo	MelaMJNV.PemoMJNVA.tblastx.out	2	true
```

First pass behavior should be linear-only and adjacent-only:

- `query_id` and `subject_id` resolve against displayed records through
  `DisplayRecordContext`.
- The two records must be adjacent after `--input_table` ordering or legacy
  input loading.
- The row is assigned to the comparison gap between those two records.
- The BLAST file itself is not coordinate-remapped. Its query and subject names
  must still match the selected records in the same way existing `-b` files do.

### Columns

Required:

| Column | Meaning |
| --- | --- |
| `query_id` | Displayed-record selector for the BLAST query-side record. |
| `subject_id` | Displayed-record selector for the BLAST subject-side record. |
| `file` | BLAST outfmt 6/7 TSV file. |

Optional:

| Column | Meaning |
| --- | --- |
| `comparison_id` | Diagnostic/debug ID for the row. Does not affect rendering in the first pass. |
| `order` | Integer ordering key for rows that share a comparison gap. |
| `enabled` | Boolean. Disabled rows are parsed and validated for syntax but not rendered. |

Rejected for now:

- Per-row `evalue`, `bitscore`, `identity`, `alignment_length`, and style
  columns. Add these only after the global filtering behavior remains stable.
- Non-adjacent record pairs. Supporting arbitrary comparison bands requires
  renderer/layout changes, not only a table adapter.
- Circular conservation rows. Add a separate mode-specific design once the
  linear adapter has real use.

### CLI/API Changes

1. Add `--blast_table BLAST_TABLE` to linear mode.
2. Reject `--blast_table` with `-b/--blast`.
3. Treat `--blast_table` as precomputed comparison data and reject combinations
   that already conflict with `blast_files`, including `--protein_blastp_mode`
   values other than `none` and collinearity comparison inputs where the current
   API already requires a single comparison source.
4. Add a `BlastTableRow` adapter in `gbdraw/cli_utils/table_adapters.py`.
5. Add `load_blast_table(path, records)` that returns comparison data grouped
   by adjacent comparison gap.
6. Prefer returning already-loaded/filter-ready comparison DataFrames instead of
   forcing a fake positional file list. The internal linear assembler already
   accepts comparison DataFrames; expose or rename that bridge clearly if the
   public API needs a less protein-specific name.
7. Preserve existing `-b` behavior exactly when `--blast_table` is absent.

### Validation Details

- Resolve `query_id` and `subject_id` with the same selector semantics used by
  `--depth_track_table`.
- Reject selectors that resolve to the same displayed record.
- Reject non-adjacent pairs with a message naming both displayed indexes and
  recommending reordering inputs or using legacy `-b` until arbitrary bands are
  implemented.
- For adjacent rows, compute `gap_index = min(query_index, subject_index)`.
- Allow multiple enabled rows for the same gap only if they are explicitly
  ordered and their files load successfully; concatenate them in table order.
- Resolve file paths relative to the BLAST table.
- Keep global BLAST filters applied in the existing `load_comparisons` /
  `filter_comparison_dataframe` path.

### Tests

- Header parser keeps `#1` selectors in `query_id` and `subject_id`.
- `input:<id>` selectors work when `--input_table` is used.
- Adjacent rows produce the same rendered comparison as equivalent `-b`.
- Non-adjacent rows are rejected.
- Ambiguous bare selectors are rejected with a namespace hint.
- `--blast_table` is mutually exclusive with `-b`.
- Multiple rows for one adjacent gap concatenate in `order`.
- Minimal committed sample TSV runs through a real `gbdraw linear` invocation.

## Phase 2: Web UI `--input_table` Staging

### Current State

`gbdraw/web/js/app/run-analysis.js` stages uploaded files and then passes
legacy input arguments:

- Circular: `--gbk /input.gb` or `--gff /input.gff --fasta /input.fasta`.
- Linear: `--gbk /seq_0.gb ...` or paired `--gff` / `--fasta` lists.
- Linear transforms still rely on positional `--record_id`, `--region`,
  `--reverse_complement`, and `--record_label` style arguments where present.

### Target State

When the bundled wheel supports `--input_table`, the web UI should stage
`/web_input_table.tsv` and pass `--input_table /web_input_table.tsv`.

The web-generated table should use deterministic IDs:

- Circular single input: `input_id=record_1`.
- Linear inputs: `record_1`, `record_2`, ...
- Expanded multi-record input remains out of scope unless the UI already models
  individual source records structurally.

Suggested circular GenBank row:

```tsv
input_id	input_type	gbk	record_id	label	order
record_1	gbk	/input.gb	#1		1
```

Suggested linear GenBank rows:

```tsv
input_id	input_type	gbk	record_id	region	reverse_complement	label	order
record_1	gbk	/seq_0.gb	#1		false		1
record_2	gbk	/seq_1.gb	#1		true		2
```

### Implementation Steps

1. Add feature detection for `--input_table` in the web option-support helper.
2. Add `buildInputTable()` and `stageInputTable()` near the existing table
   builders in `run-analysis.js`, using the same `cleanTsvCell()` helper.
3. Stage input files first, then stage `/web_input_table.tsv` with paths to the
   staged files.
4. Replace legacy input arguments only when support is detected. Keep the
   legacy path as a fallback for older wheels if the current compatibility
   policy still requires it.
5. Update depth table staging to prefer `record_id=input:<input_id>` instead of
   `#1` / `#<index>` when an input table is staged.
6. Keep circular `region` and `reverse_complement` cells empty because circular
   CLI support currently rejects them.
7. Add static packaging tests that assert `/web_input_table.tsv` is staged and
   `--input_table` is pushed for both linear and circular paths.

## Phase 3: Web UI `--track_table` Staging

### Current State

The web UI has structured linear and circular track-slot editors, but
`run-analysis.js` serializes them back into compact CLI strings:

- `--linear_track_slot ...`
- `--linear_track_axis_index ...`
- `--circular_track_slot ...`
- `--circular_track_axis_index ...`

The CLI now accepts `--track_table` and `--track_table_axis_before`, so the web
path can avoid fragile string serialization.

### Target State

When Custom Track Slots are enabled and the bundled wheel supports
`--track_table`, stage `/web_track_table.tsv` and pass:

```text
--track_table /web_track_table.tsv
--track_table_axis_before <slot_id>
```

The axis boundary must be converted from the web editor's numeric axis index to
the enabled slot ID after order normalization and disabled-row filtering.

### Implementation Steps

1. Add feature detection for `--track_table` and
   `--track_table_axis_before`.
2. Add `buildTrackTable(slots, { circular })`.
3. Reuse normalized `linearTrackSlots` / `circularTrackSlots` from the existing
   editor logic. Do not rebuild layout order independently.
4. Emit one row per slot with these common fields:
   `slot_id`, `renderer`, `order`, `side`, `track_id`, `track_index`,
   `height`, `radius`, `width`, `spacing`, `z`, and `enabled`.
5. Emit renderer-specific fields already represented in web state:
   `nt`, `dinucleotide`, `source_index`, `lane_direction`, and
   `tick_label_layout`.
6. For depth slots, prefer `track_id` derived from the same depth logical index
   used by `/web_depth_track_table.tsv` (`depth`, `depth_1`, `depth_2`, ...).
   Keep `track_index` as a checked compatibility alias when it is present in
   web state.
7. If `--track_table_axis_before` is used, omit `side` cells or ensure they do
   not conflict with the derived side. The simpler first pass is to leave `side`
   empty for slots governed by the axis boundary.
8. Keep the existing compact slot-string fallback for older wheels only if
   compatibility support is still required.
9. Update `tests/test_web_packaging.py` to check staging, axis-boundary
   conversion by slot ID, disabled-row filtering, and depth `track_id`
   relationships.

## Phase 4: Optional Web UI `--blast_table` Staging

Do this only after Phase 1 is complete.

For linear mode, the web UI already stages uploaded BLAST files as positional
arguments. Replace that with `/web_blast_table.tsv` when support is available.

Suggested table rows:

```tsv
query_id	subject_id	file	order	enabled
input:record_1	input:record_2	/blast_0.tsv	1	true
input:record_2	input:record_3	/blast_1.tsv	2	true
```

Implementation notes:

- Use input IDs from `/web_input_table.tsv` when input-table staging is active.
- Preserve legacy `-b` fallback for older wheels if needed.
- Keep web-generated rows adjacent-only until CLI support expands beyond
  adjacent comparison gaps.

## Phase 5: Documentation and Examples

Update docs only after the CLI behavior and web staging tests are stable:

- `docs/CLI_Reference.md`: add `--blast_table` option and schema.
- `docs/RECIPES.md`: add one linear `--blast_table` example.
- `docs/TUTORIALS/2_Comparative_Genomics.md`: explain when to use
  `--blast_table` instead of `-b`.
- `docs/TUTORIALS/3_Advanced_Customization.md`: update web/table notes if the
  web UI begins staging input and track tables.
- Add committed sample BLAST table files under `examples/` or
  `tests/test_inputs/`, using existing BLAST examples where possible.

Avoid duplicating the full schemas in every document. Keep the reference in
`docs/CLI_Reference.md`, with short examples elsewhere.

## Verification Matrix

Run focused tests as each phase lands:

```bash
python -m pytest tests/test_cli_table_arguments.py tests/test_comparisons.py -q
python -m pytest tests/test_web_packaging.py -k "web_run_analysis or track_table or input_table or blast_table" -q
```

Before publishing the follow-up as complete, run:

```bash
python -m pytest tests/ -v -m "not slow"
python -m pytest tests/test_output_comparison.py -v
```

If web staging changes are broad, also run the offline web verification after
preparing the browser wheel:

```bash
python tools/prepare_browser_wheel.py
python tools/verify_gui_offline.py
```

## Suggested Implementation Order

1. Implement and test `--blast_table` for linear adjacent comparisons.
2. Add `--blast_table` examples and docs.
3. Stage `/web_input_table.tsv` in the web UI and update depth table selectors
   to use `input:<id>`.
4. Stage `/web_track_table.tsv` in the web UI and remove compact slot-string
   generation from the supported-wheel path.
5. Optionally stage `/web_blast_table.tsv` in the web UI.
6. Run the broader verification matrix and update this plan's status.

## Open Questions

- Should the public API grow a clearly named `comparison_dataframes` parameter
  so table-derived BLAST data does not travel through the existing
  `protein_comparisons` compatibility name?
- Should multiple BLAST rows for one adjacent pair concatenate silently, or
  require an explicit shared `comparison_id`?
- Should circular conservation get its own table argument later, or should
  `--blast_table` grow mode-specific columns for conservation reference-side
  selection?
- How long should the web UI keep fallback compact slot-string generation for
  older browser wheels?
