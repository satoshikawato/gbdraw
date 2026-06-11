# CLI Table Arguments Remaining Work Plan

Status as of 2026-06-11: planning document for hardening and extending the
already implemented CLI table arguments.

The table-argument feature is no longer a first-pass proposal. The current
branch already includes:

- `--input_table`
- `--depth_track_table`
- `--track_table`
- `--track_table_axis_before`
- linear `--blast_table`
- web staging for `/web_input_table.tsv`, `/web_depth_track_table.tsv`,
  `/web_track_table.tsv`, and `/web_blast_table.tsv`

This document defines the next implementation steps for the remaining gaps. The
goal is to keep the design aligned with KISS, SOLID, and DRY while avoiding a
second manifest system or a parallel layout engine.

## Design Principles

### KISS

- Keep table files as CLI-boundary adapters, not as a project manifest.
- Prefer narrow mode-specific extensions over a generic all-purpose table
  schema.
- Keep rejected states explicit until their semantics are clear.
- Add one user-visible capability at a time, with focused tests around that
  capability.

### SOLID

- Keep TSV parsing, selector resolution, row normalization, and diagram assembly
  as separate responsibilities.
- Add new behavior by extending small adapters instead of changing renderers
  unless layout semantics truly require renderer changes.
- Preserve existing CLI shortcuts and public API inputs as stable interfaces.
- Normalize compatibility aliases at the boundary, then pass one internal shape
  deeper into the code.

### DRY

- Reuse `read_headered_tsv_table`, `DisplayRecordContext`, and
  `SourceRecordContext` for all table-bound selectors.
- Reuse existing comparison filtering, depth normalization, and track-slot
  parsing instead of adding table-specific duplicates.
- In the web UI, derive table staging and legacy fallback arguments from the
  same normalized row models.
- Keep documentation status in one current plan and avoid expanding outdated
  implementation notes.

## Scope

### In Scope

1. Align stale documentation with the current implementation.
2. Harden linear `--blast_table` within its current adjacent-gap model.
3. Clarify and document the `comparison_dataframes` API path.
4. Reduce drift between web table staging and legacy argument fallback.
5. Prepare circular `--input_table` transform support only after its semantics
   are explicit and testable.

### Out of Scope

- A new project manifest format.
- Non-adjacent linear comparison bands in the same pass as BLAST table
  hardening.
- Circular conservation through `--blast_table` before the circular
  conservation data model is designed.
- Per-row BLAST styling before global filtering and current rendering remain
  stable.
- Rewriting existing track layout or comparison rendering solely to support
  tables.

## Phase 1: Documentation Alignment

### Problem

`docs/CLI_TABLE_ARGUMENTS_PLAN.md` still describes `--blast_table` and web
`--input_table` / `--track_table` staging as deferred, while the current branch
implements them.

### Implementation

1. Update the current-status section in `docs/CLI_TABLE_ARGUMENTS_PLAN.md`.
2. Keep the older first-pass rationale, but mark deferred items that are now
   complete as superseded by the follow-up implementation.
3. Keep user-facing schema details in `docs/CLI_Reference.md`,
   `docs/RECIPES.md`, and tutorials rather than duplicating full manuals in the
   plan file.
4. Add a short link from the old plan to this remaining-work plan.

### Acceptance Criteria

- No plan document says `--blast_table` is not implemented without clearly
  qualifying that statement as historical.
- The current limitation list is consistent across the plan, CLI reference, and
  tutorials.
- Documentation-only changes do not require SVG reference updates.

## Phase 2: BLAST Table Hardening

### Current Model

`--blast_table` is a linear-only adapter. Rows resolve `query_id` and
`subject_id` against displayed records, require adjacent displayed records, and
produce gap-aligned comparison DataFrames.

### Implementation

1. Keep adjacency as a hard validation rule.
2. Add targeted tests for ambiguous bare selectors, disabled malformed rows,
   empty BLAST files, and duplicate `comparison_id` diagnostics if those cases
   are not already covered.
3. Confirm that global `--evalue`, `--bitscore`, `--identity`, and
   `--alignment_length` are applied exactly once downstream.
4. Document that `comparison_id` is diagnostic only in the current model.
5. Keep reversed adjacent rows implemented as query/subject column swaps only;
   do not remap coordinates beyond that existing behavior.

### Acceptance Criteria

- Sparse comparison gaps remain gap-aligned.
- Multiple enabled rows in one gap require `order`.
- Non-adjacent rows fail with a row-numbered diagnostic.
- Legacy `-b/--blast` behavior is unchanged when `--blast_table` is absent.

## Phase 3: Public API Cleanup

### Problem

The API now has `comparison_dataframes` and the older
`protein_comparisons` compatibility path. Both names should not become parallel
concepts.

### Implementation

1. Treat `comparison_dataframes` as the canonical public API name for
   precomputed linear comparison DataFrames.
2. Keep `protein_comparisons` as a backward-compatible alias.
3. Normalize the alias once near the API boundary and pass one internal
   variable deeper into assembly.
4. Add or update API documentation and tests showing that passing both names is
   rejected.

### Acceptance Criteria

- Code paths below the API boundary reason about one precomputed-comparison
  variable.
- Existing callers using `protein_comparisons` continue to work.
- New docs and examples prefer `comparison_dataframes`.

## Phase 4: Web Argument Model Consolidation

### Problem

The web UI can stage table files when the bundled wheel supports the relevant
CLI options, but legacy fallback paths can still drift if they are maintained
as independent argument builders.

### Implementation

1. Keep one normalized input-row model for staged paths, `input_id`, row-local
   `record_id`, region, reverse-complement, label, and order.
2. Keep one normalized track-slot model for renderer, side, geometry, depth
   `track_id`, and axis placement.
3. Keep one normalized BLAST-row model for adjacent paths and input selectors.
4. Generate table TSVs from those models when support is detected.
5. Generate legacy CLI arguments from the same models when support is not
   detected.
6. Keep circular multi-record canvas on the legacy path until table semantics
   cover that mode.

### Acceptance Criteria

- Static packaging tests assert that table staging and legacy fallback both use
  the normalized model names.
- Depth table selectors use `input:<input_id>` when an input table is staged and
  displayed-record indexes otherwise.
- The web path never passes table and mutually exclusive legacy arguments for
  the same logical data.

## Phase 5: Circular Input Table Transforms

### Problem

Circular `--input_table` currently rejects `region` and
`reverse_complement`. Enabling those cells without a clear model could break
coordinate-sensitive features.

### Decision Gate

Do not implement circular transforms until these questions are answered:

1. Should circular region cropping behave like linear cropping or preserve full
   source coordinates in labels and metadata?
2. How should reverse-complement transforms interact with depth tracks and
   conservation rings?
3. Should multi-record circular canvas accept a table row per displayed record,
   or should it keep the existing positional options?

### Implementation After Gate

1. Reuse existing record transform helpers where possible.
2. Apply transforms before displayed-record selector contexts are created for
   downstream depth and track tables.
3. Add focused tests for labels, depth rows, and conservation inputs that depend
   on record identity and orientation.
4. Keep unsupported combinations rejected with clear row-numbered diagnostics.

### Acceptance Criteria

- Circular transforms are either fully documented and tested or remain
  explicitly rejected.
- No transform path duplicates linear parsing logic.
- Multi-record circular behavior is documented before it is changed.

## Phase 6: Future Extensions

Treat these as separate features, not as table-adapter cleanup:

1. Non-adjacent linear comparison bands.
2. Circular conservation table inputs.
3. Per-row BLAST thresholds.
4. Per-row comparison style or color mode.
5. A saved project/session manifest.

Each item needs its own layout/API design before implementation. The table
argument layer should only adapt user input into that design after it exists.

## Verification Plan

Run focused checks after each implementation phase:

```bash
python -m pytest tests/test_cli_table_arguments.py -q
python -m pytest tests/test_web_packaging.py -k "web_run_analysis or track_table or input_table or blast_table" -q
python -m pytest tests/test_api_library_usage.py -q
```

Run broader checks when changing circular transforms, comparison rendering, or
track layout:

```bash
python -m pytest tests/ -m "not slow" -q
```

Only update SVG references when a user-visible diagram output change is
intentional and reviewed.

## Implementation Order

1. Documentation alignment.
2. BLAST table hardening tests and diagnostics.
3. API documentation and alias cleanup.
4. Web normalized-row consolidation.
5. Circular transform decision and implementation, if approved.
6. Larger comparison/conservation extensions as separate plans.

This order keeps the low-risk correctness work first and avoids expanding the
table layer beyond its responsibility.
