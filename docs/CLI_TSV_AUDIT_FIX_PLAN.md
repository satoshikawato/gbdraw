[Home](./DOCS.md) | [CLI Reference](./CLI_Reference.md)

# CLI TSV Audit Fix Plan

Date: 2026-07-12

This plan addresses four defects found in the table-backed CLI inputs:

- `--records_table`
- `--conservation_table`
- `--circular_track_table`

The fixes should preserve valid existing manifests, reject ambiguous manifests
at the table boundary, and keep session save/restore behavior aligned with direct
CLI execution.

## Goals

- Prevent `params` from silently replacing structural circular track columns.
- Keep every `records_table.region` value scoped to the record selected by its
  own row.
- Make `order` sorting correct for every accepted positive integer.
- Accept UTF-8 TSV files with or without a byte order mark (BOM).
- Add regression coverage at parser, CLI integration, and session round-trip
  boundaries.

## Non-goals

- Do not add new TSV columns or CLI options.
- Do not change the supported circular track renderers or their
  renderer-specific parameters.
- Do not add fixed empty-cell placement to the circular multi-record grid.
- Do not change the syntax or behavior of the non-table `--region` option.
- Do not change diagram geometry or reference SVG output.

## Confirmed Defects

### 1. Circular track `params` can replace structural columns

`read_circular_track_table()` converts each row to an inline slot string and
appends `params` after values from `id`, `renderer`, `side`, `r`, `w`, and the
other structural columns. The generic slot parser then interprets structural
keys found in `params` and lets the later values win.

For example, this row is accepted by the table reader:

```tsv
id	renderer	side	params
original	ticks	inside	id=replaced,renderer=spacer,side=outside
```

It becomes a slot named `replaced` with renderer `spacer` and side `outside`.
The table reader validates the slot without its computed Axis index, so an Axis
conflict can remain hidden until diagram assembly.

### 2. A row-scoped `region` can target another row

Selectorless region values are correctly qualified with the current loaded
record index. Already-qualified values are passed through unchanged. This
allows a cell in one row to crop a record loaded by another row:

```tsv
gbk	region
a.gbk	#2:1-10
b.gbk	
```

The documented row model says that `region` applies only to the record selected
by its own row. A cross-row selector therefore violates the TSV contract and
can silently crop the wrong record.

### 3. Large explicit `order` values cross the missing-value sentinel

Rows without `order` currently receive a synthetic sort value based on
`1_000_000_000`. Because positive `order` values have no upper bound, a valid
explicit value above that sentinel can sort after a row whose `order` is blank.

### 4. UTF-8 BOM-prefixed headers are rejected

CLI tables are opened with `encoding="utf-8"`. A BOM therefore remains attached
to the first header name, producing values such as `\ufeffgbk`, which are
rejected as unknown columns. The session restore path reads and rewrites TSV
headers separately and must use the same BOM policy.

## Required Behavior Decisions

Implement the following rules explicitly so parser and documentation behavior
remain deterministic.

### Circular track structural ownership

The dedicated table columns own slot structure. Reject generic slot keys in
`params` instead of applying precedence.

At minimum, reserve these keys and aliases:

- Identity and renderer: `id`, `renderer`, `type`
- Placement: `side`
- Geometry: `r`, `radius`, `w`, `width`, `spacing`, `inner_gap_px`,
  `outer_gap_px`
- Layering: `z`, `z_index`, `zindex`
- Generic state: `enabled`, `show`, `visible`, `strict`, `compress`, `reserve`

For `renderer=features`, also reject `lane_direction` and `lanes`. The table's
`side` column and Axis row already determine the feature lane. Continue to
allow renderer-specific parameters such as `nt`, `positive_color`,
`negative_color`, `legend_label`, and `tick_label_layout`.

Report a `ValidationError` containing the table path, row number, `params`
column, and conflicting key.

### Region scoping

Require `records_table.region` cells to use selectorless region syntax, such as
`1000-9000`, `1000..9000`, or `1000-9000:rc`. Reject record ID, record index,
and file selectors in a table cell. Do not silently strip or reinterpret an
explicit selector.

After validation, qualify every non-empty region with the loaded index of its
own sorted row. This preserves the one-row-one-displayed-record contract and
makes a copied cross-row selector fail loudly.

### Mixed explicit and missing `order` values

Sort rows with this tuple instead of using a numeric sentinel:

```python
(
    item.order is None,
    item.order if item.order is not None else item.row_index,
    item.row_index,
)
```

This defines the existing intended behavior without a hidden upper bound:

- Rows with explicit `order` come first and sort numerically.
- Equal explicit values retain source row order.
- Rows without `order` follow and retain source row order.
- If all values are blank, source row order is unchanged.

Do not introduce an arbitrary maximum for `order`.

### BOM handling

Read CLI manifest files with `encoding="utf-8-sig"`. This decoder accepts both
BOM-prefixed and ordinary UTF-8 files. Continue writing normalized restored
tables as UTF-8 without a BOM.

Apply the read policy to both:

- `gbdraw/io/cli_tables.py::_read_tsv_table()`
- `gbdraw/session_io.py::_rewrite_tsv_path_cells()`

Do not broaden this change to unrelated TSV formats unless their own parser
tests demonstrate the same requirement.

## Implementation Steps

### Phase 1: Harden the shared table parser

Update `gbdraw/io/cli_tables.py`:

1. Add a focused helper that parses and validates circular track `params` with
   `split_kv_list()` before constructing an inline slot specification.
2. Reject structural keys and feature lane keys according to the decisions
   above.
3. Keep the existing Axis-row checks, but reuse the parsed parameter list where
   practical to avoid parsing the same cell repeatedly.
4. Parse the completed slot specifications once, then validate them with
   `normalize_circular_track_slots_with_axis(..., axis_index)` before returning
   `CircularTrackTable`. Wrap failures in a table-level `ValidationError`.
5. Reject qualified `region` cells and qualify every accepted region with the
   current loaded row index.
6. Replace the `1_000_000_000` sorting sentinel with the tuple sort key defined
   above.
7. Change the shared table reader to `utf-8-sig`.

Keep these changes in the shared parser rather than adding defensive branches
to `circular.py` and `linear.py`. Both CLI modes should receive the same parsed
`RecordsTable` contract.

### Phase 2: Align session restoration

Update `gbdraw/session_io.py`:

1. Read restored TSV files with `utf-8-sig` before replacing dependency paths.
2. Preserve the current UTF-8, tab-delimited, LF-terminated output.
3. Verify that the first path column in a BOM-prefixed records or conservation
   table is rewritten, not left pointing to its original location.

No session schema migration should be necessary. The stored table and
dependency bindings remain unchanged.

### Phase 3: Add regression tests

Extend `tests/test_cli_tables.py` with focused parser tests:

- Accept a BOM-prefixed records table.
- Accept a BOM-prefixed conservation table and circular track table.
- Reject `id`, `renderer`, `type`, `side`, geometry, layering, and generic state
  keys in circular track `params`.
- Reject `lane_direction` and `lanes` for feature rows.
- Continue accepting representative renderer-specific parameters.
- Verify that Axis-aware normalization occurs inside the table reader.
- Reject record index, record ID, and file selectors in `region` cells.
- Confirm that selectorless regions are qualified after `order` sorting.
- Confirm correct ordering with an explicit value greater than
  `1_000_000_000` and a missing value.
- Confirm stable source order for equal explicit values and missing values.

Add CLI integration tests in the existing mode-specific files:

- Circular and linear `--records_table` runs apply selectorless regions to the
  intended row after sorting.
- A qualified region fails before rendering in both modes.
- A conflicting circular track `params` cell fails with table path, row, and
  column context.
- A valid circular track table reaches diagram assembly with the expected Axis
  index and slot order.

Extend `tests/test_session_io.py`:

- Capture and restore a BOM-prefixed records table whose first column contains
  a path dependency.
- Assert that the restored first-column path is rewritten to the materialized
  file.
- Parse the restored table again with `read_records_table()`.

Prefer small temporary TSV files and lightweight records. These regressions do
not require new reference SVGs.

### Phase 4: Update user-facing documentation

After the code and tests pass, update the input format descriptions in:

- `docs/CLI_Reference.md`
- `docs/TUTORIALS/5_Table_Driven_Inputs.md`

Document these points:

- UTF-8 files with or without a BOM are accepted.
- `region` cells are row-scoped and must not include a record or file selector.
- Structural circular track settings belong in dedicated columns and cannot be
  repeated in `params`.
- Explicit `order` values sort before blank values; blank values preserve table
  row order.

Check `docs/RECIPES.md` for examples affected by these clarifications. Do not
change examples that are already valid.

## Verification Matrix

Run the focused checks first:

```bash
pytest tests/test_cli_tables.py -v
pytest tests/test_session_io.py -v -k "table or tsv"
pytest tests/test_circular_track_slots.py -v
pytest tests/test_circular_multi_canvas.py -v
pytest tests/test_linear_selectors.py -v
```

Then run the fast suite and lint:

```bash
pytest tests/ -v -m "not slow"
ruff check gbdraw/ --select=E,F,W --ignore=E501,W503
```

Compare the documented options with current CLI help after the documentation
updates:

```bash
python -m gbdraw.cli circular --help
python -m gbdraw.cli linear --help
```

## Acceptance Criteria

- A valid CLI manifest produces the same row, ring, and track ordering as
  before these fixes.
- BOM-prefixed and ordinary UTF-8 manifests parse identically.
- Session restoration rewrites first-column dependencies in BOM-prefixed
  tables.
- Circular track `params` cannot alter slot identity, renderer, side, geometry,
  layering, enabled state, or feature lane.
- Circular track table errors are raised during table parsing with table, row,
  and column context rather than later during diagram assembly.
- Every non-empty `records_table.region` applies only to the record selected by
  its own row.
- Every accepted positive `order` value sorts consistently without a numeric
  cutoff.
- Focused tests, the fast suite, and Ruff pass.
- No reference SVG changes are required.

## Compatibility and Rollout Risks

- Manifests that currently place structural keys in `params` will become
  invalid. This is intentional because those keys can silently contradict the
  dedicated columns. The error should state which column to use instead.
- Manifests that use qualified `region` cells will become invalid. This aligns
  behavior with the documented row-scoped contract and prevents silent
  cross-row edits.
- Changing BOM handling is backward compatible for ordinary UTF-8 files.
- Replacing the numeric sort sentinel changes behavior only for mixed tables
  that use explicit `order` values at or above the old sentinel range. The new
  result follows the documented positive-integer model.

## Recommended Rollout Order

1. Add failing parser regression tests for all four confirmed defects.
2. Implement parser hardening and make those tests pass.
3. Add and pass CLI integration tests.
4. Align session restore BOM handling and pass session round-trip tests.
5. Run the focused test matrix.
6. Update CLI reference and table-driven tutorial wording.
7. Run the fast suite and Ruff.

Keep the parser, session, and documentation changes in separate commits when
possible so each compatibility decision remains easy to review.
