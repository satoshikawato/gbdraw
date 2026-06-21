# Orthogroup Multi-HSP Coverage Fix Plan

## Status

Planned.

## Background

The `2026-06-21_majanivirus_Orthogroups.gbdraw-session.json` investigation
showed that `BDT62853.1` is a real CDS feature annotated as
`wsv343-like protein`, but it is not a member of `og_82:wsv343-like protein`.

The immediate reason is deterministic and reproducible:

- `BDT62853.1` is not part of any `anchor_core_v1` anchor edge.
- Its only candidate orthogroup is `og_82`.
- The highest-scoring support row to `og_82` is `BDT62853.1 -> BDT62565.1`.
- That row has `alignment_length=1172`, `query_length=4741`,
  `subject_length=4468`, so `min_coverage=0.2472`.
- `_row_domain_only()` treats `min_coverage <= 0.25` as domain-only.
- `_row_supports_membership()` then rejects the row because membership requires
  `min_coverage >= 0.30` and not domain-only.

The user's external BLASTP output shows four HSPs between the same protein pair:

```text
BDT62853.1 BDT62565.1 58.191 1172 ... q=1200..2318 s=1039..2075
BDT62853.1 BDT62565.1 55.136 1032 ... q=33..1034   s=48..964
BDT62853.1 BDT62565.1 64.078 824  ... q=2774..3577 s=2524..3242
BDT62853.1 BDT62565.1 48.235 340  ... q=2417..2709 s=2118..2436
```

Those HSPs cover much more than one local domain when merged by interval union.
The current implementation makes two independent assumptions that are too
coarse for this case:

1. Orthogroup membership coverage is calculated from one outfmt6 row, i.e. one
   HSP, not from all HSPs for a directional protein pair.
2. `_build_core_support_candidate()` first picks the highest normalized-score
   evidence between a candidate protein and a core, then tests whether that row
   supports membership. If the highest row is domain-only, a lower-scoring but
   membership-qualified row to the same core is not considered.

Separately, both CLI and web LOSATP blastp paths currently request one HSP per
subject:

- CLI: `run_losatp_blastp()` always passes `-max_hsps_per_subject 1`.
- Web: `buildLosatArgs()` always passes `--max-hsps-per-subject 1` for blastp.

That acquisition cap prevents the inference layer from seeing the evidence
needed to merge HSP coverage.

## Goals

- Correctly distinguish local domain-only hits from multi-HSP, broad protein
  similarity.
- Keep the conservative intent of `domain_only_max_min_coverage=0.25` and
  `min_membership_min_coverage=0.30`.
- Keep product names out of membership inference. `wsv343-like protein` remains
  display metadata, not clustering evidence.
- Preserve deterministic, explainable `anchor_core_v1` decisions.
- Share the fix between CLI, API, and web Pyodide paths.
- Avoid widening unrelated pairwise display behavior.

## Non-Goals

- Do not lower coverage thresholds to make this one protein pass.
- Do not assign by matching `product`, `gene`, or label strings.
- Do not introduce a new user-facing orthogroup mode.
- Do not add an external clustering dependency.
- Do not rewrite the full orthogroup inference model in this change.
- Do not use live NCBI BLAST queries in tests.

## Design Principles

### SOLID

- Keep HSP aggregation separate from parsing, threshold filtering, and
  orthogroup graph construction.
- Keep membership evidence selection separate from related-edge diagnostics.
- Avoid adding UI-specific logic to `gbdraw.analysis.protein_colinearity`.

### KISS

- Aggregate only what is necessary for this bug: interval-union coverage by
  `(query, subject)`.
- Use the best HSP as the representative row for scoring and display in the
  first implementation. This avoids recalibrating normalized-score thresholds.
- Store aggregate diagnostics as extra columns on the normalized evidence row
  instead of introducing a large new public object model.

### DRY

- Put the aggregation logic in one Python helper used by both orthogroup and
  collinear conversion paths.
- Let the web path continue to pass raw LOSATP output into Python conversion,
  so CLI/API/web share the same membership semantics.
- Reuse the existing comparison column schema for rendering, with added derived
  columns only inside inference.

## Proposed Implementation

### 1. Add HSP Interval Aggregation

Add internal helpers in `gbdraw/analysis/protein_colinearity.py`.

Proposed helpers:

```python
def _coverage_interval_from_hsp(start: object, end: object, length: int) -> tuple[int, int] | None:
    ...

def _merge_coverage_intervals(intervals: Sequence[tuple[int, int]]) -> tuple[tuple[int, int], ...]:
    ...

def _covered_length(intervals: Sequence[tuple[int, int]]) -> int:
    ...

def _aggregate_hsps_by_protein_pair(
    hits: DataFrame,
    protein_map: Mapping[str, CdsProtein],
) -> DataFrame:
    ...
```

Behavior:

- Group rows by directional `(query, subject)`.
- Drop groups whose query or subject is unknown.
- For each HSP, read `qstart/qend` and `sstart/send`.
- Treat BLAST coordinates as inclusive one-based positions.
- Normalize reversed coordinates with `min(start, end)` and `max(start, end)`.
- Clamp coordinates to protein length.
- Ignore an HSP interval if coordinates are missing, non-numeric, or empty after
  clamping.
- Merge intervals independently on the query and subject axes.
- Compute:
  - `hsp_count`
  - `query_covered_length`
  - `subject_covered_length`
  - `query_coverage`
  - `subject_coverage`
  - `min_coverage`
  - `representative_alignment_length`
  - `total_hsp_alignment_length`
  - `coverage_source="hsp_union"`
- Pick a representative row by a stable raw-HSP rank:
  1. Higher bitscore.
  2. Lower evalue.
  3. Higher identity.
  4. Longer alignment length.
  5. Stable query/subject coordinate order.
- Keep the representative row's outfmt6 columns for downstream display and
  edge metadata.

Important scoring decision:

- `normalized_score` should continue to use the representative row's bitscore
  for the initial fix.
- Aggregate coverage should replace single-HSP coverage for domain-only and
  membership coverage decisions.

Rationale:

- Coverage answers "how much of each protein is supported?"
- Score answers "how strong is the best aligned segment?"
- Combining HSP bitscores would require new calibration of local thresholds and
  could inflate repetitive low-complexity proteins. That can be considered later
  as a versioned inference change, not as part of this fix.

### 2. Integrate Aggregation Into Normalization

Update `_normalize_directional_hit_table()` so it aggregates HSP rows before
derived coverage fields are assigned.

Current flow:

```text
raw outfmt6 rows -> one normalized row per HSP
```

New flow:

```text
raw outfmt6 rows
-> threshold-filtered HSP rows
-> one aggregated evidence row per directional protein pair
-> normalized score and aggregate coverage columns
```

Preserve current threshold semantics:

- `filter_protein_hits_by_thresholds()` remains per-HSP.
- Aggregation uses only HSP rows that pass the user's `evalue`, `bitscore`,
  `identity`, and `alignment_length` filters.

This avoids surprising users who set a minimum HSP quality threshold.

### 3. Prefer Membership-Qualified Evidence When Assigning To A Core

Refactor `_best_evidence_between_protein_and_members()`.

Current behavior:

```text
pick highest normalized_score row
then decide whether it supports membership
```

New behavior:

```text
find best membership-supporting row
also keep best diagnostic row
```

Proposed internal shape:

```python
@dataclass(frozen=True)
class _BestCoreEvidence:
    support_score: float
    support_row: object | None
    support_query_id: str
    support_subject_id: str
    diagnostic_score: float
    diagnostic_row: object | None
    diagnostic_query_id: str
    diagnostic_subject_id: str
```

Rules:

- Membership scoring uses rows where `_row_supports_membership(row)` is true.
- Diagnostic related-edge reporting may still use the best row even when it is
  domain-only or otherwise below membership coverage.
- `_build_core_support_candidate()` should calculate `same_score` and
  `cross_score` from support rows, not from domain-only rows.
- If no support row exists, return a diagnostic candidate with
  `high_confidence_pass=False` and `low_confidence_pass=False`, so related
  evidence can still be emitted as `domain_only`, `related_homolog`, or
  `ambiguous_paralog`.

This fixes the second root cause: a domain-only high-score row can no longer
mask a lower-scoring valid membership row to the same core.

### 4. Remove The One-HSP Acquisition Cap For Orthogroup Inference

CLI/API:

- Add `max_hsps_per_subject: int | None = 1` to `run_losatp_blastp()`.
- Add the same parameter to `_run_losatp_search()`.
- Keep `max_hsps_per_subject=1` for pairwise display mode.
- Pass `max_hsps_per_subject=None` for:
  - `build_rbh_orthogroup_protein_blastp_comparisons()`
  - `build_orthogroup_collinearity_blocks()`

When the value is `None`, do not add `-max_hsps_per_subject` to the LOSATP
command. Let LOSATP/BLASTP return multiple HSPs under its default behavior.

Web:

- Update `buildLosatArgs()` in `gbdraw/web/js/app/run-analysis.js`.
- Keep `--max-hsps-per-subject 1` only for plain pairwise blastp display.
- Omit the HSP cap when `useOrthogroupBlastp` or `useCollinearBlastp` is true.

Cache behavior:

- The LOSAT raw cache key already includes the command args.
- Removing the HSP cap changes the cache key for orthogroup/collinear blastp,
  so existing one-HSP raw cache entries will not be reused for recomputation.
- No session schema bump is required for raw cache safety.
- Existing imported sessions will still display their saved SVG and
  `orthogroupState` until the user reruns analysis. That is acceptable and
  should be mentioned in release notes.

### 5. Keep Display Rendering Stable

The renderer expects one comparison row per displayed link. It does not need to
draw every HSP.

For orthogroup display edges:

- Use the representative HSP row's coordinates for the ribbon/link.
- Store aggregate coverage diagnostics in orthogroup metadata where useful.
- Do not duplicate ribbons for every HSP in the initial fix.

Potential follow-up:

- Add a popup detail section showing `hsp_count`, query union coverage, subject
  union coverage, and whether membership used aggregate coverage.

## Test Plan

### Unit Tests For HSP Aggregation

Add tests near existing protein/orthogroup tests in `tests/test_collinearity.py`
or `tests/test_protein_colinearity.py`.

1. `test_hsp_union_coverage_uses_merged_intervals`
   - Build two synthetic `CdsProtein` objects with lengths similar to the
     `BDT62853.1`/`BDT62565.1` case.
   - Provide four HSP rows for one `(query, subject)` pair.
   - Assert the aggregated row has:
     - `hsp_count == 4`
     - `query_coverage > 0.30`
     - `subject_coverage > 0.30`
     - `min_coverage > 0.30`
     - `domain_only is False`

2. `test_overlapping_hsps_do_not_double_count_coverage`
   - Provide overlapping HSPs for one pair.
   - Assert coverage equals interval union length, not the sum of
     `alignment_length`.

3. `test_invalid_hsp_coordinates_fall_back_safely`
   - Provide rows with invalid coordinates and a valid representative HSP.
   - Assert invalid intervals do not crash and do not inflate coverage.

### Anchor-Core Membership Tests

1. `test_anchor_core_assigns_member_with_multi_hsp_union_coverage`
   - Create an established two-record core.
   - Add an unassigned long protein whose best relationship to the core is split
     across multiple HSP rows.
   - Assert the protein becomes a member when aggregate `min_coverage >= 0.30`.

2. `test_anchor_core_prefers_membership_support_over_domain_only_top_hit`
   - Create an established core.
   - Add an unassigned protein with:
     - a higher normalized-score domain-only row to one core member;
     - a lower normalized-score but membership-qualified row to another core
       member.
   - Assert the protein is assigned to the core.
   - Assert the domain-only row may still appear as related diagnostic evidence,
     but it does not block membership.

3. `test_anchor_core_keeps_true_domain_only_hit_unassigned`
   - Provide only overlapping low-coverage HSPs.
   - Assert the candidate remains unassigned and related evidence is marked
     `domain_only`.

### Search Command Tests

1. CLI command construction:
   - Monkeypatch `subprocess.run`.
   - Assert pairwise blastp still passes `-max_hsps_per_subject 1`.
   - Assert orthogroup and collinear blastp do not pass that flag.

2. Web argument construction:
   - Add or extend a JS/module test if a suitable harness exists.
   - Assert `buildLosatArgs()` omits `--max-hsps-per-subject` for orthogroup and
     collinear modes.
   - Assert plain pairwise mode keeps the existing cap.

### Regression Fixture

Use synthetic data for CI. Do not require live BLAST.

Optional manual validation:

- Rerun the Majanivirus session after the fix with orthogroup blastp.
- Confirm `BDT62853.1` becomes a member of the `wsv343-like protein`
  orthogroup if LOSATP returns the additional HSPs.
- Confirm `og_82` remains separate from unrelated broad WSSV-like families.

## Implementation Sequence

1. Add HSP interval helpers and focused unit tests.
2. Integrate aggregation into `_normalize_directional_hit_table()`.
3. Update domain-only and membership-support tests to use aggregate coverage.
4. Refactor `_best_evidence_between_protein_and_members()` to separate support
   rows from diagnostic rows.
5. Add anchor-core regression tests for:
   - multi-HSP broad coverage;
   - domain-only top row masking valid support;
   - true domain-only non-membership.
6. Add `max_hsps_per_subject` plumbing for CLI/API LOSATP calls.
7. Update web blastp argument construction.
8. Run focused tests:

   ```bash
   python -m pytest tests/test_protein_colinearity.py tests/test_collinearity.py -q
   python -m pytest tests/test_web_packaging.py -k "web_run_analysis or web_config" -q
   ```

9. Run broader fast tests if the focused suite passes:

   ```bash
   python -m pytest tests/ -m "not slow" -q
   ```

## Risks And Mitigations

### Risk: More HSP output increases browser runtime or memory use

Mitigation:

- Remove the one-HSP cap only for orthogroup/collinear inference.
- Keep pairwise display mode unchanged.
- Measure LOSAT timing logs before and after on representative sessions.
- If needed, add an internal acquisition cap greater than one, but keep it
  separate from membership thresholds and document it as a performance limit.

### Risk: Aggregate coverage admits repeat-rich false positives

Mitigation:

- Use interval union, not raw summed HSP lengths.
- Keep score based on the best representative HSP for the initial fix.
- Continue requiring local threshold support and best/second-core separation.
- Keep domain-only diagnostic evidence visible without granting membership.

### Risk: Existing tests assume one row per protein pair

Mitigation:

- Aggregate only inside normalization/inference.
- Keep parsing and rendering row schemas stable.
- Add tests that assert display rows still contain standard comparison columns.

### Risk: Saved sessions contain old one-HSP cache entries

Mitigation:

- Cache keys include command args, so rerun will miss old capped entries and
  recompute.
- Saved SVG/orthogroup state remains unchanged until rerun.
- Add release-note text explaining that old sessions should be regenerated to
  receive improved multi-HSP orthogroup membership.

## Acceptance Criteria

- Multiple HSPs for one directional protein pair are represented as one
  aggregate inference evidence row with interval-union coverage.
- `_row_domain_only()` uses aggregate `min_coverage` when aggregate evidence is
  available.
- A high-score domain-only row cannot mask a lower-score membership-qualified
  row to the same core.
- Orthogroup/collinear blastp paths request multiple HSPs; plain pairwise
  display remains one-HSP by default.
- Existing anchor-core protections still pass:
  - same-record hits do not merge distinct cross-record cores;
  - ambiguous paralogs are not assigned;
  - true domain-only relationships remain non-members.
- The Majanivirus `BDT62853.1` case has a clear manual validation path and
  should be assigned to `og_82` when multi-HSP evidence is present and passes
  aggregate coverage.
