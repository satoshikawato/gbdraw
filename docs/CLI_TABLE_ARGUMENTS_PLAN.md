# CLI Table Arguments Implementation Plan

Implemented for the first-pass `--input_table`, `--depth_track_table`, and
`--track_table` scope. `--blast_table` remains intentionally deferred.

This document defines the first implementation plan for adding table-based CLI
inputs that replace fragile positional argument groups with headered TSV files.
The goal is to match the existing `--default_colors`, `--feature_table`,
`--label_table`, and `--qualifier_priority` style rather than introducing a new
`manifest` concept.

The agreed BLAST entry point is `--blast_table`. This plan focuses on the three
additional table arguments:

- `--depth_track_table`
- `--input_table`
- `--track_table`

## Current Status and Next Direction

Status as of 2026-06-07: the first-pass CLI-boundary implementation for
`--input_table`, `--depth_track_table`, and `--track_table` has landed and now
has committed example TSV files, minimal real CLI smoke tests, and user-facing
documentation in the CLI reference, recipes, and advanced customization
tutorial. This document remains the engineering tracking record; the detailed
sections below describe the intended model and should not be expanded into a
second user manual.

Implemented first-pass pieces:

- Shared headered TSV helpers were added in `gbdraw/cli_utils/tables.py`,
  including comment handling that preserves `#1` record selectors, required and
  optional column validation, row-numbered diagnostics, table-relative path
  resolution, primitive cell conversion, and displayed/source record selector
  contexts.
- Table adapters were added in `gbdraw/cli_utils/table_adapters.py` for input
  rows, depth-track rows, and custom track-slot rows.
- Linear and circular CLI commands now accept `--input_table`,
  `--depth_track_table`, `--track_table`, and `--track_table_axis_before`, with
  mutual-exclusion checks against legacy positional arguments that describe the
  same logical data.
- `--input_table` supports `gbk` and `gff` rows, row-local `record_id`
  selection, `input_id` aliases, exact one-row-to-one-record loading by default,
  explicit `expand_records`, labels and ordering, and linear row-local
  `region`/`reverse_complement` transforms. Circular mode rejects region and
  reverse-complement cells.
- `--depth_track_table` supports sparse per-record depth assignments keyed by
  logical `track_id`, wildcard rows, track-scoped metadata conflict detection,
  mode-specific `track_height`/`track_width` validation, and named depth-track
  resolution in linear/circular drawing and legends.
- Shared depth-axis calculation is now keyed by logical `track_id` for named
  depth inputs instead of compact record-local list position.
- `--track_table` converts rows to the existing linear and circular track-slot
  dataclasses, preserves named depth relationships through `track_id`, filters
  disabled rows before axis-boundary resolution, and supports
  `--track_table_axis_before`.
- Circular depth slots can inherit `track_width` from `--depth_track_table`
  metadata when no explicit track-table width is provided.
- The web UI depth-track staging path now writes `/web_depth_track_table.tsv`
  and invokes the CLI with `--depth_track_table`, removing the need for
  placeholder depth-file arguments and legacy per-depth-track metadata
  arguments in that path.
- Focused tests were added for parser behavior, table adapters, sparse named
  depth selection, CLI argument forwarding, circular width application, and web
  packaging expectations.
- Committed sample TSV files were added under `examples/` for linear and
  circular input tables, sparse depth-track tables, custom track tables, and
  small depth TSV inputs.
- Minimal real CLI smoke tests now run the committed sample TSV files through
  `gbdraw linear` and `gbdraw circular`.
- User-facing documentation now covers the new table arguments in
  `docs/CLI_Reference.md`, practical examples in `docs/RECIPES.md`, and a short
  conceptual explanation in `docs/TUTORIALS/3_Advanced_Customization.md`.

Remaining or intentionally deferred pieces:

- `--blast_table` is not implemented yet. It should be handled as a separate
  change after the selector and table conventions are stable.
- The web UI only stages a depth-track table. It does not yet emit
  `--input_table` or `--track_table` from structured UI state.
- If this feature needs additional hardening, move those items into smaller
  follow-up issues rather than expanding this plan. Keep this file as the
  implementation note and source of design rationale.

Recommended next direction:

1. Treat `--blast_table` as a separate feature after the selector conventions
   have had time to settle.
2. Consider a web UI cleanup that emits `--input_table` and `--track_table`
   from structured UI state, matching the existing depth-track table staging.
3. Keep this file as the engineering tracking record and maintain user-facing
   details in `docs/CLI_Reference.md`, `docs/RECIPES.md`, and the tutorials.

## Goals

1. Let users describe sparse per-record inputs without placeholder files or
   order-dependent empty slots.
2. Keep existing positional options as backward-compatible shortcuts.
3. Use readable, headered TSV files with stable column names.
4. Keep CLI and web UI mental models aligned: users add files and settings to a
   specific record or track, not to an implicit slot in a list.
5. Convert table inputs into the existing internal structures where possible,
   instead of rewriting diagram assembly in the first pass.
6. Preserve user-facing IDs such as `input_id` and `track_id` until the point
   where the renderer genuinely needs an index. Do not collapse named table
   relationships into positional lists earlier than necessary.

## Design Guardrails

Keep the table implementation aligned with SOLID, KISS, and DRY:

- Treat each table as a CLI-boundary adapter. A table parser should produce
  small row dataclasses, a table-specific normalizer should convert those rows
  to existing internal objects, and diagram assembly should keep owning layout
  behavior.
- Do not add a generic project manifest, global session object, or second layout
  resolver. The first pass should add narrow adapters around existing loaders,
  depth specs, and track slot dataclasses.
- Share only stable primitives: headered TSV parsing, path resolution, scalar
  conversion, boolean conversion, row-numbered diagnostics, and selector syntax
  parsing. Do not hide table-specific rules in a large generic parser.
- Keep selector parsing separate from selector resolution. A `#1` selector in a
  displayed-record table means the first displayed record, while `#1` in an
  `--input_table` row means the first record inside that row's source file. Do
  not force these two meanings through one resolver.
- Prefer explicit validation over warning-and-ignore behavior when a table cell
  looks meaningful but cannot be honored in the current mode.
- Avoid merging two sources that describe the same logical data. Table columns
  should own table-scoped metadata; legacy per-track or positional metadata
  should be rejected when it would compete with a table column.

## General Table Rules

All new `*_table` inputs should follow the same parser behavior.

- TSV files are headered.
- Blank lines are ignored.
- Comment lines are ignored, but `#1`-style record selectors are data, not
  comments. Treat a line as a comment only when the first non-whitespace text is
  exactly `#`, starts with `# `, starts with `#<TAB>`, or starts with `##`. Do
  not treat `#1<TAB>...` or `#2<TAB>...` as comments.
- Inline `#` comments are not supported. This avoids corrupting valid values
  such as `#RRGGBB` colors and `#1` record selectors.
- Required columns must be present, but optional columns may be absent.
- Unknown columns are rejected in the first implementation. This catches typos.
- Empty optional cells mean "use the CLI/global/config default" unless a
  table-specific normalizer explicitly defines inheritance from wildcard rows.
- File paths are resolved relative to the table file directory unless absolute.
- Boolean cells accept `1`, `0`, `true`, `false`, `yes`, `no`, `on`, and `off`.
- Numeric cells are validated with the same rules as the equivalent CLI option.
- Record selector syntax should be parsed through one shared CLI helper, but
  each table should resolve selectors against the correct context.
- For displayed-record tables, selectors should accept explicit namespace
  prefixes:
  - `#1` for a 1-based displayed-record index.
  - `input:<input_id>` for IDs attached by `--input_table`.
  - `record:<record.id>` for loaded record IDs.
  - `name:<record.name>` for loaded record names.
  - `file:<alias>` for source file aliases already attached by loaders.
  - `accession:<alias>` for accession-like aliases where practical.
- Bare displayed-record selectors remain supported as a convenience. They are
  resolved across the same namespaces, but ambiguity is defined by the final
  displayed record, not by the number of matching aliases. If `input:ref` and
  `record:ref` both point to the same displayed record, accept the selector. If
  the same bare value would resolve to two or more different displayed records,
  reject it and recommend an explicit namespace form, for example `input:ref` or
  `record:ref`.
- For `--input_table`, `record_id=#1` means the first record inside that row's
  input file before the row is converted to displayed records. It must be
  resolved by a row-local source-record context, not by
  `DisplayRecordContext`. The row-local context may support `#1`,
  `record:<record.id>`, `name:<record.name>`, and `accession:<alias>`, but it
  must not support `input:<input_id>` because displayed records do not exist yet.
- `*` is reserved as a wildcard selector where the table argument explicitly
  supports applying one row to every loaded record. When wildcard rows and
  specific selector rows describe the same logical target, the table-specific
  normalizer must define whether specific rows override wildcard defaults or
  whether the combination is rejected.

Implementation should add a small shared helper, for example
`gbdraw/cli_utils/tables.py`, with low-level parsing and conversion helpers:

- `read_headered_tsv_table(path, required, optional, table_name)`
- `resolve_table_path(table_path, value)`
- `parse_optional_float_cell(value, field_name, positive=False, nonnegative=False)`
- `parse_optional_int_cell(value, field_name, positive=False, nonnegative=False)`
- `parse_optional_bool_cell(value, field_name)`

Do not use ad hoc string splitting in each CLI module. Prefer Python's `csv`
module with `delimiter="\t"` so cells remain raw strings. If pandas is used for
an existing table path, disable automatic NA conversion and dtype inference so
empty cells, `NA`, accession-like values, and leading zeroes are not silently
rewritten. Keep this shared helper focused on syntax, path resolution, and
primitive cell conversion. Table-specific meaning should live in small
table-specific normalizers, for example `DepthTrackTableRow`, `InputTableRow`,
and `TrackTableRow`, rather than in one generic manifest parser.

## Shared Implementation Model

Keep the implementation small, but introduce a few narrow shared concepts so the
tables do not drift into separate selector and ID semantics:

- `TableParseError` or `ValidationError` messages should include the table name,
  row number, and column name where possible.
- `DisplayRecordContext` should wrap the loaded displayed records and expose a
  shared `resolve_display_record_selector(selector)` helper. This is where
  `input_id`, `record.id`, `record.name`, file aliases, accession aliases, and
  `#index` matching live. It should understand the namespace prefixes described
  above, accept multiple aliases that identify the same displayed record, reject
  bare selectors that identify multiple different displayed records, and produce
  diagnostics that name the matching namespaces. Table-specific code should not
  implement its own record matching.
- `SourceRecordContext` should wrap records parsed from one `--input_table` row
  before they become displayed records. It may share the same parsed selector
  dataclass as `DisplayRecordContext`, but it must not share displayed-record
  index semantics. This keeps `record_id=#1` row-local in `--input_table` while
  preserving `#1` as displayed-record-local in downstream tables.
- `InputTableRow` loading should attach the final displayed `input_id` to each
  displayed record's annotations, for example `gbdraw_input_id`, before later
  table normalizers run. For exact rows this is the row's `input_id`; for
  expanded rows it is the derived per-record ID. This lets
  `--depth_track_table` and the future `--blast_table` use stable table IDs
  instead of fragile GenBank IDs.
- Input loading should use a small row-loader protocol keyed by `input_type`,
  for example `load_source_records(row) -> list[SeqRecord]`. First-pass loaders
  should include `gbk` and `gff`, and each loader should own only source parsing
  plus source-file diagnostics. Selection, explicit expansion, labels, regions,
  and reverse-complement transforms remain in the table normalizer so row
  semantics stay consistent across loader types.
- Track lookup should have the same shape for every depth input source:
  table-derived sparse depth inputs and legacy dense depth inputs should both be
  converted at the CLI boundary into a stable ordered list of `track_id` values.
  Track-slot normalization keeps named depth relationships named until the
  selected diagram path genuinely needs a numeric index. Sparse table rows must
  not be compacted into record-local list positions before track identity is
  resolved.
- Use one narrow named-depth bridge for both table-derived and legacy dense
  depth inputs. Legacy `--depth_track` inputs can be adapted to the same
  `track_id` model with generated IDs such as `depth`, `depth_1`, and
  `depth_2`, while table inputs keep the user-provided `track_id` values.
- Add a shared depth-slot resolver, for example
  `resolve_depth_track(depth_tracks_by_id, slot_params)`, and use it in linear
  drawing, circular drawing, and legend generation. The resolver should prefer
  `track_id` when present and fall back to legacy `track_index` only for dense
  inputs by mapping the index through the ordered depth `track_id` list. Diagram
  code should not duplicate `track_index` parsing in several places for
  table-aware paths; this resolver is the compatibility boundary.
- Compute shared depth axes by logical `track_id`, not by record-local list
  position. This preserves the meaning of `--share_depth_axis` when some
  records are missing earlier logical tracks.

These helpers are not a project manifest. They are small adapters at the CLI
boundary. The public diagram assembly and rendering code should remain the
source of layout behavior.

## Compatibility Policy

Table arguments should be mutually exclusive with legacy positional arguments
that describe the same logical data.

Recommended first-pass rules:

- `--depth_track_table` cannot be combined with `--depth`, `--depth_track`, or
  legacy per-depth-track metadata arguments such as `--depth_track_label`,
  `--depth_track_color`, `--depth_track_height`,
  `--depth_track_large_tick_interval`, `--depth_track_small_tick_interval`, and
  `--depth_track_tick_font_size`.
- `--input_table` cannot be combined with `--gbk`, `--gff`, `--fasta`,
  `--record_id`, `--reverse_complement`, `--region`, or `--record_label`.
- `--track_table` cannot be combined with `--linear_track_order`,
  `--linear_track_slot`, `--circular_track_order`, or `--circular_track_slot`.
- `--track_table` cannot be combined with `--linear_track_axis_index` or
  `--circular_track_axis_index`. Use the table-specific
  `--track_table_axis_before` argument when a table needs the existing
  axis-boundary behavior.

Global styling options still apply as defaults when they do not compete with a
table column. For example, `--depth_min`, `--depth_max`,
`--share_depth_axis`, `--depth_log_scale`, `--depth_color`, `--depth_height`,
`--evalue`, and `--identity` remain global filtering or rendering defaults
unless the table schema later adds per-row overrides. If a table column is
non-empty, the table value wins over the global default.

For table-derived sparse depth inputs, `--share_depth_axis` keeps its global
meaning but must group data by logical `track_id`. It must not group by
record-local list position, because different records may be missing different
tracks.

Table arguments that describe different logical data may be combined. For
example, `--input_table` may be combined with `--depth_track_table` and
`--track_table`; after `--input_table` loads records, downstream record
selectors are resolved against the resulting displayed records and their
`input_id` aliases.

## `--depth_track_table`

### User Model

Users describe depth files as rows attached to records and logical depth tracks.
Missing rows mean the record has no data for that track. There are no placeholder
arguments.

Each row has two kinds of meaning:

- `record_id` + `file` assign one depth file to one displayed record and one
  logical `track_id`.
- `track_*` columns describe the logical track itself. They are intentionally
  track-scoped, not record-scoped. If users need different styling for two
  records, they should use different `track_id` values.

Linear example:

```tsv
record_id	track_id	file	track_label	track_color	track_height	track_large_tick_interval	track_small_tick_interval	track_tick_font_size
record:rec2	depth_A	rec2.A.depth.tsv	Sample A	#4A90E2	12			8
record:rec2	depth_B	rec2.B.depth.tsv	Sample B	#E45756	28	50	10	8
record:rec4	depth_A	rec4.A.depth.tsv	Sample A	#4A90E2	12			8
```

Circular example:

```tsv
record_id	track_id	file	track_label	track_color	track_width	track_large_tick_interval	track_small_tick_interval	track_tick_font_size
input:rec2	depth_A	rec2.A.depth.tsv	Sample A	#4A90E2	12			8
input:rec2	depth_B	rec2.B.depth.tsv	Sample B	#E45756	28	50	10	8
input:rec4	depth_A	rec4.A.depth.tsv	Sample A	#4A90E2	12			8
```

### Scope

Supported in both linear and circular mode.

For linear mode, `track_height` maps to the per-depth-track height currently
exposed by `--depth_track_height`.

For circular mode, use `track_width` for radial ring width. A non-empty
`track_height` cell is an error in circular mode, because accepting and ignoring
it would make a meaningful user setting disappear. A non-empty `track_width`
cell is an error in linear mode.

### Columns

Required:

| Column | Meaning |
| --- | --- |
| `record_id` | Displayed-record selector. Use `*` to apply the row to all loaded records. With `--input_table`, prefer `input:<input_id>`. |
| `track_id` | Logical depth track ID. Stable across records. |
| `file` | Depth TSV in samtools depth format. |

Optional:

| Column | Meaning |
| --- | --- |
| `track_label` | Legend label for the logical track. |
| `track_color` | Fill color for the logical track. SVG name or `#RRGGBB`. |
| `track_height` | Linear depth track height in pixels. |
| `track_width` | Circular depth ring width in pixels. |
| `track_large_tick_interval` | Large y-axis tick interval. |
| `track_small_tick_interval` | Small y-axis tick interval. |
| `track_tick_font_size` | Tick label font size. |
| `order` | Integer ordering key for logical tracks. Defaults to first appearance. |

Validation:

- `track_id` must not be empty.
- For each `track_id`, non-empty `track_label`, `track_color`, mode-valid size
  metadata (`track_height` for linear, `track_width` for circular), and tick
  metadata must be consistent across rows after wildcard fallback is resolved.
  Reject conflicting values instead of silently choosing one.
- Depth table metadata is track-scoped, not record-scoped. Wildcard rows may
  provide default files for records, but `track_*` metadata defines the logical
  track. A specific record row may repeat the same metadata for readability, but
  it must not silently override a different value. Conflict diagnostics should
  explain that per-record depth styling is not supported and suggest a distinct
  `track_id` when that appears to be the user's intent.
- Empty metadata cells do not erase a metadata value established by another row.
  If all rows for a `track_id` leave a metadata cell empty, use the CLI/global/
  config default for that field.
- Wildcard rows (`record_id=*`) are defaults. A specific record row for the same
  `track_id` overrides the wildcard file for that record. Duplicate wildcard
  rows for the same `track_id` are errors, and duplicate specific
  `(record_id, track_id)` pairs are errors.
- A `file` cell may not be empty.
- Missing records are an error unless `record_id` is `*`.
- Selectors that resolve to more than one displayed record are errors. Multiple
  namespace aliases that identify the same displayed record are accepted; aliases
  that identify different displayed records are rejected with diagnostics that
  recommend explicit selector prefixes.
- `track_id` values define the public logical depth-track names. If a custom
  `--track_table` contains a depth slot with `track_id`, it must resolve to one
  of these names.
- Sparse logical tracks must keep their `track_id` identity even when earlier
  tracks are missing for a displayed record. Do not represent a record that only
  has `depth_B` as a one-element list where `depth_B` becomes local index 0.
- Linear mode rejects non-empty `track_width` cells. Circular mode rejects
  non-empty `track_height` cells.

### Internal Mapping

After records are loaded:

1. Read and validate the table into `DepthTrackTableRow` values.
2. Resolve record selectors through `DisplayRecordContext`. Wildcard rows are
   expanded after the displayed-record set is known.
3. Split table meaning into two narrow concepts:
   track-scoped metadata keyed by `track_id`, and record-scoped file assignments
   keyed by `(displayed_record_index, track_id)`. This split should happen in the
   depth table normalizer, not in diagram assembly.
4. Resolve wildcard file defaults first, then apply specific record file rows as
   overrides for the same `track_id`. Reject duplicates within the wildcard
   layer or within the specific layer.
5. Determine logical track order by `order`, then first appearance.
6. Preserve the ordered `track_id` values as the internal depth track IDs. These
   IDs are user-facing and may be referenced by diagnostics, legend/debug output,
   default table-derived depth slots, and `--track_table` depth `track_id`
   fields.
7. Return a narrow bundle that keeps sparse data keyed by name. The bundle is
   the single CLI-boundary bridge for all CLI depth input; legacy dense
   `--depth` and `--depth_track` inputs should be adapted to the same shape with
   generated track IDs, rather than keeping a parallel depth-selection
   implementation.
   For example:

   ```python
   @dataclass(frozen=True)
   class DepthTrackMetadata:
       label: str | None
       fill_color: str | None
       height: float | None
       width: float | None
       large_tick_interval: float | None
       small_tick_interval: float | None
       tick_font_size: float | None


   @dataclass(frozen=True)
   class DepthTrackBundle:
       track_ids: tuple[str, ...]
       metadata_by_track_id: Mapping[str, DepthTrackMetadata]
       files_by_record_and_track_id: tuple[Mapping[str, str], ...]
   ```

   The bundle is a CLI adapter output, not a public manifest.
8. Materialize `DepthTrackSpec` values into record-indexed mappings, not compact
   record-local lists. A record that only has `depth_B` should be represented as
   `{"depth_B": DepthTrackSpec(...)}`, never as a one-element list where
   `depth_B` implicitly becomes index 0. Existing sequence-based public API paths
   may keep a compatibility adapter, but table-derived CLI paths must remain
   name-keyed until slot resolution.
9. Precompute depth data into the same record-indexed mapping shape, for
   example `tuple[Mapping[str, DepthTrackData], ...]`. Keep a small adapter for
   legacy public APIs that still expect sequences, but table-derived assembly
   should select tracks by `track_id`.
10. Add a shared resolver such as:

   ```python
   def resolve_depth_track(
       depth_tracks_by_id: Mapping[str, DepthTrackData],
       slot_params: Mapping[str, object],
       ordered_track_ids: Sequence[str],
   ) -> DepthTrackData | None:
       ...
   ```

   It should prefer `slot_params["track_id"]`, and use `track_index` only as a
   compatibility lookup into `ordered_track_ids`. Use this resolver in linear
   depth drawing, circular depth drawing, and depth legend generation.
11. When `--share_depth_axis` is enabled and no explicit `--depth_max` is set,
   calculate shared maxima per `track_id`. Do not group by local list index,
   because sparse table rows can place different logical tracks at the same
   compacted position.
12. For linear mode, store `track_height` on the `DepthTrackSpec` or equivalent
   table-derived depth data metadata. For circular mode, keep `track_width` as
   `DepthTrackMetadata.width` only until track slots are derived; then apply it
   to `CircularTrackSlot.width` for default depth slots or for a `--track_table`
   depth slot that leaves `width` empty. Do not add circular ring geometry to
   `DepthTrackSpec`, because depth specs should describe data and per-track
   rendering metadata, while circular slot objects own radial geometry.
13. Set `show_depth=True`.
14. Continue through the existing depth precomputation path.

This keeps the renderer changes small while preserving the user's logical track
names. The table helper should be the only place that turns sparse,
record-bound rows into named depth track objects. Track-slot normalization
should not discard `track_id`; the common resolver is the only place that maps a
legacy numeric `track_index` to a logical depth track.

### CLI Changes

Add to both `gbdraw/linear.py` and `gbdraw/circular.py`:

```text
--depth_track_table DEPTH_TRACK_TABLE
    Headered TSV assigning depth files and per-track metadata to records.
```

Update argument validation so `--show_depth` accepts `--depth_track_table` as an
input source. Also reject `--depth_track_table` with legacy `--depth_track_*`
metadata arguments listed in the compatibility policy.

### Web UI Changes

The web UI already behaves close to the desired model by letting users attach
depth files to sequence cards. Change the generated CLI arguments from
placeholder-based `--depth_track` groups to a staged TSV:

1. Stage each uploaded depth file at a stable virtual path.
2. Generate `/web_depth_track_table.tsv`.
3. Push `--depth_track_table /web_depth_track_table.tsv`.

The visible UI should continue to say "add depth to this record" rather than
"fill track slots".

### Tests

- Table parser accepts valid depth table with comments and relative paths.
- Table parser treats `#1` in the `record_id` column as a record selector, not a
  comment line.
- Linear CLI converts sparse rows into the expected named
  record-by-`track_id` bundle.
- Sparse rows keep logical `track_id` identity when earlier tracks are missing
  for a record; a `depth_B` slot must never render `depth_A` or vice versa due to
  compacted list positions.
- Table-derived depth data are precomputed and passed through assembly as
  record-indexed `track_id` mappings, not compacted record-local lists.
- `--share_depth_axis` shares maxima between rows with the same `track_id` and
  does not share axes between different logical tracks that happen to occupy the
  same compacted position.
- Wildcard rows provide defaults, and specific record rows override those
  defaults for the same `track_id`.
- Empty metadata cells inherit track metadata or global defaults according to
  the documented metadata rules; they do not erase non-empty track metadata.
- Circular CLI accepts the same logical table shape, preserves track labels and
  colors, accepts `track_width`, applies table-derived width to circular depth
  slots, and rejects non-empty `track_height`.
- Table `track_id` values are preserved in normalized depth specs and can be
  matched by custom depth slots through a `track_id` field.
- Linear drawing, circular drawing, and legend generation all resolve depth
  slots through the shared `track_id`/legacy `track_index` resolver.
- `record_id` values can match `input_id` aliases when `--input_table` is used.
- Ambiguous record selectors are rejected with a diagnostic that names the table
  row.
- Reject duplicate specific `(record_id, track_id)` rows and duplicate wildcard
  rows for the same `track_id`.
- Reject conflicting metadata for the same `track_id`.
- Reject use with `--depth`, `--depth_track`, or legacy `--depth_track_*`
  metadata arguments.
- Web packaging test checks `--depth_track_table` staging.

## `--input_table`

### User Model

Users describe the input records and per-record transforms in one table instead
of maintaining parallel `--gbk`, `--gff`, `--fasta`, `--record_id`,
`--reverse_complement`, `--region`, and `--record_label` lists.
The `input_id` is the stable CLI-table name for the displayed record and may be
used by later record-bound tables such as `--depth_track_table` and
`--blast_table`.

The default model is one table row to one displayed record, because that keeps
`input_id`, labels, regions, and downstream joins unambiguous. The first pass
also supports an explicit `expand_records=true` escape hatch for users who mean
"load every record from this source file", matching an existing CLI habit
without making expansion implicit.

Example for linear GenBank inputs:

```tsv
input_id	input_type	gbk	record_id	region	reverse_complement	label
ref	gbk	ref.gb	#1		false	Reference
sample_a	gbk	sample_a.gb	#1	1000-50000	true	Sample A
sample_b	gbk	sample_b.gb	#2		false	Sample B
```

Example for linear GFF3 + FASTA inputs:

```tsv
input_id	input_type	gff	fasta	record_id	region	reverse_complement	label
ref	gff	ref.gff	ref.fna	#1		false	Reference
sample_a	gff	sample_a.gff	sample_a.fna	#1	1000-50000	true	Sample A
```

Example mixing GenBank and GFF3 + FASTA rows:

```tsv
input_id	input_type	gbk	gff	fasta	record_id	label
ref	gbk	ref.gb			#1	Reference
sample_a	gff		sample_a.gff	sample_a.fna	#1	Sample A
```

Example expanding every record from one GenBank file:

```tsv
input_id	input_type	gbk	expand_records
plasmids	gbk	plasmids.gb	true
```

If the file contains records `pA` and `pB`, the displayed `input_id` aliases are
`plasmids:pA` and `plasmids:pB`. If a source record ID is empty or duplicated
within the expansion, fall back to row-local index aliases such as
`plasmids:#1`.

Example for circular GenBank inputs:

```tsv
input_id	input_type	gbk	record_id	label	order
ref	gbk	ref.gb	#1	Reference	1
sample_a	gbk	sample_a.gb	#1	Sample A	2
```

### Scope

First implementation:

- Support linear mode with `record_id`, `region`, `reverse_complement`, and
  `label`.
- Support circular mode with a minimal stable-ID subset: `input_id`,
  `input_type`, source file columns, `record_id`, `label`, `order`, and
  `expand_records`.
- Support mixed `gbk` and `gff` rows. Each row chooses a small source loader by
  `input_type`, and downstream row semantics are shared.
- Treat one table row as one displayed record by default. If a row without
  `record_id` would load more than one record from its source file or GFF3/FASTA
  pair, reject it and ask the user to add `record_id` or set
  `expand_records=true`.
- When `expand_records=true`, one row may produce multiple displayed records.
  The row's `input_id` becomes a stable prefix used to derive per-record
  `input_id` aliases.
- Allow the same source file to appear in multiple rows. Exact rows can select
  different records or regions from one multi-record file without relying on
  positional side effects, while explicit expansion rows can intentionally load
  every record from that file.

Follow-up implementation:

- Add circular `region`, reverse-complement, and multi-record placement support
  after the basic row-local loader and selector semantics are stable.
- Add label templating for expanded rows if users need custom labels such as
  `{input_id} {record_id}`. In the first pass, non-empty `label` is accepted only
  for one-row-to-one-record rows.
- Consider circular multi-record placement columns such as `row` after the basic
  loader path is stable. The existing `order` column controls display order in
  the first pass; it does not replace circular row placement.

The reason for keeping circular transforms minimal is that stable `input_id`
aliases are needed by downstream tables in both modes, while circular region and
reverse-complement semantics require loader/API changes and should be a
separate, explicit expansion.

### Columns

Required:

| Column | Meaning |
| --- | --- |
| `input_id` | User-facing displayed-record ID. Used for diagnostics and later table joins. |
| `input_type` | `gbk` or `gff`. |

Conditionally required:

| Column | Required when | Meaning |
| --- | --- | --- |
| `gbk` | `input_type=gbk` | GenBank/DDBJ flatfile. |
| `gff` | `input_type=gff` | GFF3 file. |
| `fasta` | `input_type=gff` | Matching FASTA file. |

Optional:

| Column | Meaning |
| --- | --- |
| `record_id` | Existing record selector syntax, including `#1`. |
| `region` | Linear-only region syntax without the record prefix, for example `1000-50000` or `1000-50000:rc`. |
| `reverse_complement` | Linear-only per-input reverse-complement flag. |
| `label` | Definition label, replacing positional `--record_label` where applicable. |
| `order` | Integer display order. Defaults to table row order. |
| `expand_records` | Boolean. When true, load every source record from this row and derive one displayed `input_id` per source record. Defaults to false. |

Validation:

- Non-expanded `input_id` values must be unique. Expanded rows use `input_id` as
  a prefix, and every derived displayed `input_id` must be unique across the
  whole table.
- `input_type` must be one of `gbk` or `gff`.
- Mixed `gbk` and `gff` rows are allowed. Each row is parsed by its
  `input_type` loader, and all rows then pass through the same selector,
  transform, label, and alias normalizer.
- Each row must resolve to exactly one displayed record unless
  `expand_records=true`. This keeps `input_id`, `label`, `region`, and
  `reverse_complement` attached to the intended record by default.
- `expand_records=true` cannot be combined with `record_id`, `region`, or
  `label` in the first pass. Use explicit rows when per-record selection,
  cropping, or labels are needed. A row-level `reverse_complement=true` may
  apply to every expanded source record in linear mode, matching the existing
  file-level CLI behavior.
- Expanded input IDs are derived as `<input_id>:<record.id>` when source record
  IDs are non-empty and unique within the expansion. Otherwise derive
  `<input_id>:#<row-local-index>`. Derived IDs are what downstream tables see
  through `input:<derived_id>` selectors.
- `record_id` is resolved only against the records parsed from that row's source
  file or GFF3/FASTA pair. It does not use displayed-record indexes and it does
  not see records produced by other rows.
- In circular mode, non-empty `region` and `reverse_complement` cells are errors
  in the first pass. This keeps stable ID loading available without implying
  transform semantics that have not been designed.
- In linear mode, `region` requires both start and end coordinates and must be
  valid according to existing coordinate rules.
- The `region` column is row-local and should not contain a record selector or
  file selector prefix. Accept `1000-50000` and `1000-50000:rc`; reject
  `rec1:1000-50000` in the first pass because the row already selects exactly
  one record.
- If `region` contains `:rc`, it should set region reverse-complement behavior.
  It should not implicitly set the row-level `reverse_complement` flag.
- In the first pass, reject rows that combine `reverse_complement=true` with a
  `region` ending in `:rc`. Supporting both is possible, but the transform order
  must be designed explicitly to avoid accidental double reverse-complementing.
- After loading, attach the final displayed `input_id` to the displayed record.
  If a later selector would match both an `input_id` and the same record's native
  `record.id`/`record.name`, accept it. If it would match different displayed
  records, reject it as ambiguous rather than preferring one namespace.

### Internal Mapping

For both modes:

1. Read and validate the input table before `validate_input_args()`, or update
   `validate_input_args()` to accept `args.input_table`. The linear command
   should branch to table loading before the legacy positional input loader, so
   table input is not rejected for missing `--gbk`/`--gff` arguments. The
   circular command should do the same when `--input_table` is present.
2. Resolve relative file paths.
3. Sort rows by `order`, then original row order.
4. Load source records row-by-row with a dedicated CLI helper such as
   `load_input_table_records(mode=...)`. Use a small row-loader registry keyed
   by `input_type` (`gbk`, `gff`) so mixed input types do not require a generic
   manifest or a second file-major loader. Reuse the existing GenBank/GFF3
   parsing, feature filtering, and mode-appropriate loader helpers where they
   fit, but keep row context until displayed records and `input_id` aliases are
   finalized.
5. Resolve `record_id` through `SourceRecordContext` when it is present. If
   `record_id` is empty and `expand_records=false`, the row must resolve to
   exactly one source record; otherwise raise a validation error asking for
   `record_id` or `expand_records=true`.
6. If `expand_records=true`, derive one displayed record per source record,
   attach derived `input_id` aliases, and reject unsupported per-record metadata
   cells (`record_id`, `region`, and `label`) in the first pass.
7. In linear mode, apply exact-row transforms in one explicit order:
   select the record from that row's source records, crop the row-local `region`
   on that selected record, then apply row-level `reverse_complement` if
   requested. For expanded rows, row-level `reverse_complement=true` applies to
   every expanded source record. Reject `reverse_complement=true` combined with
   `region:rc` in the first pass.
8. In linear mode, implement row-local region cropping through a small helper
   that operates on
   the already-selected record. Do not synthesize global `--region` strings or
   rely on displayed-record indexes, because table row order and duplicate
   source files make those indexes easy to misapply.
9. In circular mode, reject transform cells not supported by the minimal subset
   before loading proceeds. Do not silently ignore them.
10. Apply `label` to exact one-row-to-one-record rows using the existing
   `gbdraw_record_label` annotation.
11. Attach the final displayed `input_id` to each displayed record, for example
   `gbdraw_input_id`.
12. Return the loaded records in sorted table order, preserving expansion order
    within each expanded row.

The first pass should not change public API objects. Keep the table logic in CLI
code or a CLI helper. Avoid a generic loader manifest; the helper should be a
small adapter for the CLI input table only, with mode-specific behavior kept
explicit.

### CLI Changes

Add to linear and circular mode:

```text
--input_table INPUT_TABLE
    Headered TSV describing input files, record selectors, regions, reverse
    complement flags, and record labels.
```

In circular mode the first implementation accepts the same table argument but
rejects non-empty `region` and `reverse_complement` cells. Help text should make
that mode-specific limitation clear.

### Web UI Changes

No first-pass UI change is required. The web UI already holds per-sequence
state structurally. A later cleanup can generate `--input_table` from that state
instead of pushing positional input arguments.

### Tests

- Linear GenBank input table loads two files in table order.
- Linear GFF3 + FASTA input table loads matching pairs.
- Mixed GenBank and GFF3 + FASTA rows load in table order through row-local
  loaders.
- Circular GenBank input table loads two files in table order, attaches
  `input_id`, and applies `label`.
- Circular input table rejects non-empty `region` and `reverse_complement`
  cells in the first pass.
- `--input_table` is accepted without `--gbk`, `--gff`, or `--fasta` in both
  modes, while still rejecting combinations with those legacy positional inputs.
- A row with `record_id=#1` selects the first record within that row's input
  file.
- A row with `record_id=#1` does not select the first displayed record produced
  by the table.
- The same GenBank file may appear in two rows selecting different records or
  regions, and each row gets its own `input_id`.
- A row without `record_id` is rejected when its input file contains multiple
  records and `expand_records` is false.
- A row with `expand_records=true` loads every source record and derives stable,
  unique displayed `input_id` aliases.
- `expand_records=true` rejects non-empty `record_id`, `region`, and `label`
  cells in the first pass.
- `reverse_complement`, `region`, and `label` stay attached to the row's one
  displayed record.
- `region` rows do not get converted to ambiguous bare loaded-record indexes.
- Row-local region cropping is applied after row-local record selection and
  before row-level reverse complementation.
- A row-local `region` cell with a record/file selector prefix is rejected in
  the first pass.
- A later depth table can target the displayed record by `input_id`.
- Selectors involving `input_id` and native record IDs are accepted when all
  aliases point to the same displayed record and rejected when they point to
  different displayed records.
- Reject use with `--gbk`, `--gff`, `--fasta`, `--record_id`,
  `--reverse_complement`, `--region`, or `--record_label`.

## `--track_table`

### User Model

Users describe custom track layout as rows instead of compact
`<slot_id>:<renderer>@key=value` strings. This is mainly for advanced users and
for the web UI, where table generation is easier to debug than CLI string
serialization.

Linear example:

```tsv
slot_id	renderer	order	side	track_id	height	spacing	z	enabled
features	features	1	overlay					true
depth_A	depth	2	below	depth_A	12px	4px	0	true
gc	gc_content	3	below		20px	4px	0	true
skew	gc_skew	4	below		20px	4px	0	true
```

Circular example:

```tsv
slot_id	renderer	order	side	track_id	radius	width	spacing	z	enabled
features	features	1	inside				4px	0	true
ticks	ticks	2	inside				2px	0	true
depth_A	depth	3	inside	depth_A	78%	12px	2px	0	true
gc	gc_content	4	inside		68%	20px	2px	0	true
```

### Scope

Supported in both linear and circular mode.

The table should be converted into the existing `LinearTrackSlot` or
`CircularTrackSlot` objects and then passed through the existing normalizer. Do
not add a second layout resolver.

Axis-boundary layout remains table-level, not row-level. The table should
usually express placement directly with the `side` column. When a table needs to
reproduce the existing axis-boundary behavior, use
`--track_table_axis_before SLOT_ID`: enabled slots before that named slot render
above/outside the axis, and the named slot plus following enabled slots render
below/inside the axis. This avoids making users count raw TSV rows after `order`
sorting and disabled-row filtering.

### Columns

Required:

| Column | Meaning |
| --- | --- |
| `slot_id` | Unique slot ID. |
| `renderer` | Renderer name or existing alias. |

Optional common columns:

| Column | Meaning |
| --- | --- |
| `order` | Integer draw/layout order. Defaults to row order. |
| `side` | Linear: `above`, `below`, `overlay`. Circular: `inside`, `outside`, `overlay`. |
| `track_id` | Logical depth track ID for `depth` slots. Preferred over `track_index` for both table-derived sparse depth inputs and legacy depth inputs adapted to the named depth bundle. |
| `track_index` | Numeric depth track index for legacy dense depth inputs; conservation track index for circular conservation if explicitly needed. Use `source_index` to choose a conservation source. |
| `height` | Linear slot height, using the same pixel scalar syntax as `height=`. Prefer explicit `px`. |
| `radius` | Circular slot radius, using the same scalar syntax as `r=` (`px`, `%`, or unitless factor). |
| `width` | Circular slot width, using the same scalar syntax as `w=`. Prefer explicit `px` for fixed ring widths. |
| `spacing` | Slot spacing, using the same scalar syntax as existing slot strings. Prefer explicit `px`. |
| `z` | z-index. |
| `enabled` | Boolean. Disabled rows are accepted but not rendered. |
| `nt` | Dinucleotide for GC content/skew renderers. |
| `dinucleotide` | Alias for `nt`, matching existing slot-string behavior. |
| `source_index` | Circular conservation source selector, matching existing slot params. |
| `lane_direction` | Circular features lane direction (`inside`, `outside`, or `split`), matching existing feature slot params. |
| `tick_label_layout` | Circular tick label layout, matching existing tick slot params. |

Mode-specific validation:

- Linear table rejects non-empty `radius` and `width` cells.
- Circular table rejects non-empty `height` cells.
- `renderer` must normalize to one of the supported renderers for the current
  mode.
- `slot_id` must be unique after sorting.
- For `depth` rows, accept either `track_id` or `track_index`. If both are
  provided, validate that the numeric index points at the same logical
  `track_id`, then keep `track_id` for table-derived sparse depth paths.
- If a depth `track_id` is provided, some depth input must be present and the ID
  must match one of the current `DepthTrackBundle.track_ids`. The bundle may
  come from `--depth`, `--depth_track`, or `--depth_track_table`.
- Prefer `track_id` over `track_index` whenever both are present. Sparse depth
  rows can make compact per-record numeric indexes misleading, and legacy dense
  depth inputs have generated stable IDs such as `depth`, `depth_1`, and
  `depth_2`.
- If the existing slot normalizer inserts a default `track_index=0` for a depth
  slot that also has `track_id`, the downstream shared depth resolver must still
  use `track_id` first. Do not let the compatibility default override the table
  relationship.
- Generic layout fields must be stored on the slot object, not in `params`.
  Reuse existing validation for this.
- Do not add a per-row `axis_index` column. Axis-boundary support is table-level
  and comes from `--track_table_axis_before`, not from repeating a value on every
  row.
- Do not invent table-only scalar syntax. If `0.78x` is desired as an alias for
  `78%` or `0.78`, add it to the shared `ScalarSpec.parse()` behavior and test
  the existing slot-string parser too.
- Disabled rows are filtered before axis derivation and normalization and do not
  reserve layout space. Therefore `--track_table_axis_before` must name an
  enabled row after `order` sorting. If it names a disabled or missing slot,
  reject the table instead of silently moving the boundary.

### Internal Mapping

1. Read the table.
2. Sort rows by `order`, then original row order.
3. Convert each row directly to `LinearTrackSlot` or `CircularTrackSlot`; do not
   serialize table rows back into slot strings just to parse them again.
4. Store non-empty renderer-specific columns in `params`. For circular
   conservation, use `source_index`; if a future `source` alias is accepted,
   convert it to `source_index` before normalization.
5. Resolve depth `track_id` values against the ordered depth-track IDs from the
   current `DepthTrackBundle`, regardless of whether the bundle came from
   `--depth`, `--depth_track`, or `--depth_track_table`. Keep
   `params["track_id"]` so assembly can select the matching depth data by name
   for each displayed record through the shared depth resolver. Set
   `params["track_index"]` only as a checked compatibility alias that does not
   replace `track_id`.
6. In circular mode, if a depth slot selects a table-derived `track_id`, the
   slot leaves `width` empty, and the depth table metadata defines
   `track_width` for that track, apply that metadata as `CircularTrackSlot.width`.
   An explicit `width` cell in `--track_table` wins over depth-table metadata.
7. Validate every row's syntax and column compatibility, but filter disabled
   rows before applying `--track_table_axis_before` and before calling the
   existing slot normalizer. The axis boundary is resolved by `slot_id` against
   the exact active slot sequence that layout will see.
8. Run the existing `normalize_*_track_slots()` validation path on the active
   slot list. If `--track_table_axis_before` is provided, resolve the named
   boundary to the corresponding active-slot index, then use the existing
   `normalize_*_track_slots_with_axis()` path after the table has been converted
   and disabled rows have been removed. If `track_id` support is added to slot
   strings too, also cover `parse_*_track_slots()` with equivalent tests.
9. Pass the resulting slots to the existing diagram API arguments:
   `linear_track_slots` or `circular_track_slots`.

This table is intentionally close to the current slot string syntax. It should
not add layout behavior that cannot also be represented by the existing slot
dataclasses. The one acceptable table-only bridge in the first pass is
depth-track name selection through the named `DepthTrackBundle`; keep that
contained in the CLI adapter and diagram assembly path unless slot-string
`track_id` support is added deliberately.

### CLI Changes

Add to linear mode:

```text
--track_table TRACK_TABLE
    Headered TSV describing linear custom track slots.

--track_table_axis_before SLOT_ID
    Axis boundary for --track_table. Enabled slots before SLOT_ID render above
    the linear axis; SLOT_ID and following enabled slots render below. The
    boundary is resolved after order sorting and disabled-row filtering.
    Explicit side cells must not conflict with the derived side.
```

Add to circular mode:

```text
--track_table TRACK_TABLE
    Headered TSV describing circular custom track slots.

--track_table_axis_before SLOT_ID
    Axis boundary for --track_table. Enabled slots before SLOT_ID render outside
    the circular axis; SLOT_ID and following enabled slots render inside. The
    boundary is resolved after order sorting and disabled-row filtering.
    Explicit side or lane cells must not conflict with the derived side.
```

Use the same CLI name in both modes; mode determines the allowed renderer set
and column validation. `--track_table_axis_before` requires `--track_table` and
cannot be combined with legacy slot/order inputs.

### Web UI Changes

The web UI can keep its structured slot editor. In a later cleanup, replace
generated `--linear_track_slot` and `--circular_track_slot` arguments with a
staged `/web_track_table.tsv`. That should reduce quoting bugs and make debug
logs easier to inspect.

### Tests

- Linear table converts to the same `LinearTrackSlot` list as equivalent
  `--linear_track_slot` strings.
- Circular table converts to the same `CircularTrackSlot` list as equivalent
  `--circular_track_slot` strings.
- Circular `width=12px` is interpreted as pixels, while `radius=78%` is
  interpreted as a factor, matching existing `ScalarSpec` rules.
- Circular conservation rows use `source_index` to choose sources; a plain
  `source` column is rejected unless an explicit alias is implemented.
- Depth rows can select a logical depth track by `track_id` for both legacy and
  table-derived depth inputs; conflicting `track_id` and `track_index` values
  are rejected.
- Sparse depth rows are selected by `track_id`, not by compact per-record
  numeric list position.
- A depth slot with both a table-derived `track_id` and a defaulted
  `track_index=0` still renders the named track selected by `track_id`.
- Circular depth slots inherit `track_width` from `--depth_track_table` metadata
  only when the `--track_table` row leaves `width` empty; explicit track-table
  width takes precedence.
- Existing renderer-specific slot params such as `lane_direction` and
  `tick_label_layout` round-trip through the table representation.
- Reject duplicate `slot_id`.
- Reject invalid mode-specific columns.
- Reject use with legacy slot/order CLI arguments.
- Reject use with legacy linear/circular track axis-index arguments.
- `--track_table_axis_before` resolves a named enabled `slot_id` after order
  sorting and disabled-row filtering, then reproduces the same side derivation
  and conflict checks as the existing mode-specific axis-index options.
- `--track_table_axis_before` rejects missing and disabled slot IDs.
- Disabled rows are validated for syntax and unknown columns, but they are not
  counted by `--track_table_axis_before` and do not reserve layout space.

## Related `--blast_table` Decision

`--blast_table` should be the user-facing table argument for BLAST inputs in
both modes:

- In linear mode, it describes pairwise comparison files between displayed
  records.
- In circular mode, it describes conservation-ring BLAST files around displayed
  references.
- Record references in the BLAST table should use the same
  `DisplayRecordContext` selector resolver as `--depth_track_table`, including
  namespace-prefixed selectors such as `input:ref`, `record:NC_000001`, and
  `#1`.

The same option name is acceptable because the command mode determines the
schema. Internally, the implementation may still use separate helpers such as
`read_linear_blast_table()` and `read_circular_blast_table()`.

Do not add both `--comparison_table` and `--conservation_table` unless a future
use case requires both names. Keep the public CLI simpler.

## Implementation Order

1. Add shared headered TSV parsing helpers, primitive cell conversion, selector
   syntax parsing, and row-numbered diagnostics.
2. Add `DisplayRecordContext` for displayed-record selector resolution, including
   explicit namespace prefixes and bare-selector diagnostics. Ambiguity should
   be based on distinct displayed records, not distinct matching aliases.
   Existing positional inputs should create the same context as table-loaded
   inputs, even if no `input_id` aliases exist yet.
3. Add `SourceRecordContext` and the small `--input_table` row-loader registry.
   Reuse the parsed selector dataclass, but do not reuse displayed-record index
   semantics.
4. Add the named depth bridge: `DepthTrackBundle`, adapters for `--depth`,
   `--depth_track`, and future `--depth_track_table`, and the shared
   `resolve_depth_track(...)` helper used by linear drawing, circular drawing,
   and legend generation.
5. Update depth precomputation, depth slot layout, legend synchronization, and
   `--share_depth_axis` handling so named paths use record-indexed `track_id`
   mappings and shared axes are computed by logical `track_id`, not by compact
   record-local list index. Keep sequence-based adapters only for existing public
   API compatibility.
6. Implement `--input_table` with row-local loaders, mixed `gbk`/`gff` support,
   exact one-row-to-one-record behavior by default, explicit
   `expand_records=true`, full transform support for exact linear rows, and the
   minimal stable-ID subset for circular mode. This unlocks
   `input:<input_id>` selectors for all downstream table arguments in both
   modes.
7. Add the depth table CLI adapter that returns `DepthTrackBundle` with ordered
   logical `track_ids`, track-scoped metadata, and sparse
   record-by-`track_id` file assignments.
8. Implement `--depth_track_table` for linear and circular modes. Include the
   small depth assembly changes needed to select table-derived depth data by
   `track_id` rather than compact record-local list position. In circular mode,
   thread `track_width` into default or table-derived `CircularTrackSlot.width`,
   not into depth data objects.
9. Add web UI staging for depth table.
10. Implement `--track_table` by converting rows directly to existing track slot
   objects, preserving depth `track_id` for both legacy and table-derived depth
   data, applying circular depth metadata widths where appropriate, and
   supporting `--track_table_axis_before` after disabled rows are filtered.
11. Add `--blast_table` in a separate change, using the same parser and selector
   conventions.

The shared selector context, named depth bridge, shared depth resolver, and
`track_id`-based shared-axis behavior should land before the first depth table
implementation. That prevents the first table feature from baking in
position-only assumptions. `--input_table` should provide at least the circular
stable-ID subset before circular `--depth_track_table` documents `input:<id>`
selectors.

## Documentation Updates

After implementation:

- Update `docs/CLI_Reference.md`.
- Add examples to `docs/RECIPES.md`.
- Add a short table-schema section to
  `docs/TUTORIALS/3_Advanced_Customization.md`.
- Update web UI help text if generated CLI arguments switch to staged tables.

## Non-Goals

- Do not remove legacy positional arguments.
- Do not make `*_table` files a full project/session format.
- Do not support JSON/YAML for these inputs in the first pass.
- Do not allow arbitrary columns silently.
- Do not make table inputs mutate global config files.
