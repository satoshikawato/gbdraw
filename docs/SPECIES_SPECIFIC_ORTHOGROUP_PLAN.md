# Species-Specific Orthogroup Inference Plan

## Status

Planned.

This plan adds OrthoFinder-style species-specific orthogroups to the current
`anchor_core_v1` inference path without weakening the existing outparalog
protection.

The short version:

1. Keep cross-record orthogroup inference exactly as the authoritative core
   inference step.
2. Assign same-record inparalogs to established cross-record cores as today.
3. After that, look only at still-unassigned proteins.
4. Build record-local paralog groups from strong same-record reciprocal or
   near-reciprocal evidence.
5. Emit those groups as orthogroups with explicit `record_local` scope.

This follows the useful part of OrthoFinder's behavior: a cluster containing
multiple genes from one species can be an orthogroup, while a singleton remains
unassigned. It does not copy OrthoFinder's full MCL or HOG pipeline.

## Background

Current `anchor_core_v1` intentionally separates cross-record and same-record
evidence:

- Cross-record reciprocal or near-reciprocal evidence creates orthogroup cores.
- Same-record evidence can assign an unassigned protein to an existing core as
  an inparalog.
- Same-record evidence must not create cores or merge two existing cores.

That design protects DnaE/PolC-like old outparalogs from being collapsed by
same-record bridges. It also means a paralog family present in only one input
record is currently left without orthogroup membership, even when the internal
evidence is strong.

OrthoFinder treats this differently. Its legacy MCL output writes clusters with
two or more genes as orthogroups, including clusters whose genes all come from
one species. Single-gene clusters are written to `Orthogroups_UnassignedGenes`.
Its documentation calls these `Species-specific orthogroups`.

## Goals

- Add orthogroup membership for strong paralog clusters that occur in exactly
  one record.
- Keep singleton genes unassigned.
- Preserve all existing cross-record orthogroup behavior and outparalog guards.
- Make record-local groups deterministic, explainable, and clearly marked in
  Python and web metadata.
- Keep record-local inference independent from the existing `include_singletons`
  compatibility path: isolated proteins stay unassigned even when older callers
  request singleton-inclusive seed groups.
- Reuse existing hit normalization, ranking, coverage, edge, path, name, and
  serialization helpers where possible.
- Avoid a new user-facing orthogroup mode or setting in the first
  implementation.

## Non-Goals

- Do not implement OrthoFinder MCL.
- Do not implement OrthoFinder HOGs or tree reconciliation.
- Do not allow same-record evidence to merge two cross-record cores.
- Do not use product names or locus names as primary grouping evidence.
- Do not create orthogroups for singleton proteins.
- Do not add a broad UI control unless regression testing shows it is needed.

## Development Principles

### SOLID

- Keep the new behavior behind small helpers with one responsibility each:
  selecting record-local evidence, forming components, checking conflicts, and
  integrating accepted components into `OrthogroupResult`.
- Do not spread record-local special cases across unrelated rendering and UI
  code. Prefer group-level metadata such as `scope`.
- Put the collinearity exclusion at the one boundary where orthogroup display
  edges are converted into adjacent display edges, not in every downstream
  collinearity consumer.

### KISS

- Implement record-local groups as a post-processing phase after
  `anchor_core_v1` core membership is stable.
- Use connected components from accepted same-record evidence edges.
- Avoid a new clustering dependency or a new inference mode.
- Do not implement singleton-specific exceptions. Components must have at least
  one accepted same-record membership edge and at least two proteins.

### DRY

- Reuse existing normalized rows from `_dedupe_anchor_core_directional_rows()`.
- Reuse existing helpers such as `_normalized_score_from_row()`,
  `_row_supports_membership()`, `_row_min_coverage()`, `_row_domain_only()`,
  `_anchor_core_hit_rank()`, `_make_ortholog_edge()`,
  `_build_ortholog_paths()`, and `_orthogroup_result_from_member_ids()`.
- Reuse existing web orthogroup payload shape and extend it with small optional
  fields instead of creating a parallel payload type.
- Reuse existing snake_case web payload keys (`member_count`,
  `record_coverage_count`) and add new optional snake_case keys rather than
  introducing a second camelCase convention.

## Terminology

- `cross_record`: an orthogroup created from cross-record core evidence.
- `record_local`: an orthogroup whose accepted members all come from one record.
- `local_paralog`: member role for a protein assigned to a `record_local`
  orthogroup.
- `unassigned`: a protein with no orthogroup membership after all inference
  phases. Singletons remain unassigned.

## Current Code Touch Points

Primary file:

- `gbdraw/analysis/protein_colinearity.py`

Relevant existing functions:

- `_select_anchor_core_edges()`
- `_derive_anchor_core_thresholds()`
- `_build_core_support_candidate()`
- `_build_anchor_core_orthogroups()`
- `_orthogroup_result_from_member_ids()`
- `_copy_orthogroup_result_with_metadata()`
- `_build_ortholog_paths()`
- `_build_adjacent_display_edges_by_pair()`
- `select_rbh_orthogroup_edges_from_directional_hits()`

Collinearity conversion:

- `gbdraw/analysis/collinearity.py`
  - `orthogroup_edges_to_lossless_collinearity_anchors()`

Web serialization:

- `gbdraw/web/js/app/python-helpers.js`
  - `_serialize_orthogroups_payload()`

Web state/UI:

- `gbdraw/web/js/services/config.js`
  - `applyOrthogroupState()`
- `gbdraw/web/js/app/orthogroups.js`

Focused tests:

- `tests/test_collinearity.py`
- `tests/test_web_packaging.py`

## Data Model Changes

Add a group scope type:

```python
OrthogroupScope = Literal["cross_record", "record_local"]
```

Extend member roles:

```python
OrthogroupMemberRole = Literal[
    "anchor",
    "coortholog",
    "inparalog",
    "low_confidence",
    "local_paralog",
]
```

Extend `OrthogroupResult`:

```python
scope_by_orthogroup_id: dict[str, OrthogroupScope] = field(default_factory=dict)
source_record_index_by_orthogroup_id: dict[str, int] = field(default_factory=dict)
```

Rules:

- Existing groups default to `cross_record`.
- New species-specific groups use `record_local`.
- `source_record_index_by_orthogroup_id` is set only for `record_local` groups.
- Older saved sessions without `scope` are interpreted as `cross_record`.
- `_orthogroup_result_from_member_ids()` and
  `_copy_orthogroup_result_with_metadata()` must carry these maps forward so
  metadata is not lost when paths, related edges, or compatibility fields are
  rebuilt.

Why group-level scope instead of a separate result type:

- Orthogroup coloring, feature indexing, popup rendering, FASTA export, and
  saved-session handling already operate on orthogroup payloads.
- A small metadata extension keeps the data model simple while preserving the
  distinction needed for UI labels and future filtering.

## Proposed Internal Constants

Keep these internal at first:

```python
_RECORD_LOCAL_NEAR_RECIPROCAL_MIN_RATIO = 0.85
_RECORD_LOCAL_MIN_MEMBERS = 2
_RECORD_LOCAL_MIN_MEMBERSHIP_MIN_COVERAGE = _MIN_MEMBERSHIP_MIN_COVERAGE
_RECORD_LOCAL_CORE_COMPETITION_RATIO = 1.10
```

Reuse:

```python
_MIN_MEMBERSHIP_MIN_COVERAGE = 0.30
_DOMAIN_ONLY_MAX_MIN_COVERAGE = 0.25
```

Do not add a user-facing setting unless real datasets show that the default is
too aggressive.

## Algorithm

### Phase 1: Build Normalized Evidence

Use the existing `anchor_core_v1` path:

1. Normalize directional hits.
2. Deduplicate strongest directional rows per protein pair.
3. Keep same-record non-self rows in `best_by_direction`.
4. Keep cross-record rows as today.

No new parsing or normalization path is needed.

### Phase 2: Infer Cross-Record Cores

Keep the current behavior:

1. `_select_anchor_core_edges()` ignores same-record rows.
2. `_build_anchor_core_orthogroups()` unions only selected cross-record anchor
   edges.
3. Same-record evidence may assign unassigned proteins to an existing core as
   `inparalog`.

Do not alter existing tests that assert same-record hits do not merge distinct
cross-record cores.

### Phase 3: Identify Record-Local Candidates

Run this after existing inparalog assignment and before final related-edge
classification.

This phase must only inspect proteins that are still absent from
`group_by_protein`. It must not reuse the `include_singletons` branch in
`_build_anchor_core_orthogroups()` as a signal that isolated proteins should
become members. Record-local inference is edge-driven, not singleton-driven.

Candidate proteins:

- are present in `protein_map`;
- are not already in `group_by_protein`;
- belong to the same record as the same-record evidence being considered;
- have at least one same-record non-self hit in `best_by_direction`.

Excluded proteins:

- already assigned to a cross-record orthogroup;
- singletons with no accepted same-record edge;
- proteins whose best evidence is domain-only;
- proteins whose accepted evidence is below membership coverage;
- proteins with competing support for an existing cross-record core.

### Phase 4: Select Record-Local Evidence Edges

Add helper:

```python
def _select_record_local_paralog_edges(
    best_by_direction: Mapping[tuple[str, str], object],
    protein_map: Mapping[str, CdsProtein],
    thresholds: Mapping[str, _LocalThreshold],
    group_member_ids: Mapping[str, set[str]],
    group_by_protein: Mapping[str, str],
) -> tuple[_AnchorCoreEvidenceEdge, ...]:
    ...
```

The helper returns canonical same-record undirected evidence edges.

An edge `p - q` passes when:

1. `p` and `q` are in the same record.
2. `p != q`.
3. neither protein is already assigned to any orthogroup.
4. both directions exist in `best_by_direction`.
5. both directional rows pass `_row_supports_membership()`.
6. both rows have `min_coverage >= _MIN_MEMBERSHIP_MIN_COVERAGE`.
7. neither row is domain-only.
8. if either endpoint has a `_LocalThreshold`, the directional score from that
   endpoint must pass that endpoint threshold.
9. if an endpoint has no cross-record threshold, the directional score from
   that endpoint must be at least
   `_RECORD_LOCAL_NEAR_RECIPROCAL_MIN_RATIO` times that endpoint's best
   same-record score among unassigned candidates.
10. tie-breaking uses `_anchor_core_hit_rank()` and stable protein sort order.

This gives OrthoFinder-like behavior:

- when cross-record context exists, use the local ortholog-scale threshold;
- when no cross-record context exists, require reciprocal or near-reciprocal
  same-record evidence on the protein's own same-record score scale.

### Phase 5: Build Record-Local Components

Add helper:

```python
def _record_local_components_from_edges(
    edges: Sequence[_AnchorCoreEvidenceEdge],
    protein_map: Mapping[str, CdsProtein],
) -> tuple[tuple[str, ...], ...]:
    ...
```

Rules:

- Use connected components of accepted same-record evidence edges.
- Keep only components with at least `_RECORD_LOCAL_MIN_MEMBERS` members.
- Every component must contain proteins from exactly one record.
- Sort components deterministically by minimum `_protein_sort_key()`.
- Sort members deterministically by `_protein_sort_key()`.

Do not include isolated proteins. They remain unassigned.

### Phase 6: Score Record-Local Components

Add a small helper that derives local support from the already accepted
record-local evidence edges:

```python
def _record_local_support_by_member(
    edges: Sequence[_AnchorCoreEvidenceEdge],
) -> dict[str, float]:
    ...
```

Rules:

- For each member, store the strongest accepted directional score touching that
  member.
- For a component, `component_local_support` is the minimum member support in
  that component. This makes the weakest accepted member the limiting evidence.
- Do not reuse `_build_core_support_candidate().support` for local component
  strength. That helper mixes cross-record and same-record support for an
  existing-core assignment; record-local support must stay on the direct local
  evidence scale.

### Phase 7: Reject Components With Cross-Core Competition

Before accepting a record-local component, check whether any member is better
explained by an existing cross-record core.

Add helper:

```python
def _record_local_component_has_competing_core_support(
    member_ids: Sequence[str],
    group_member_ids: Mapping[str, set[str]],
    best_by_direction: Mapping[tuple[str, str], object],
    thresholds: Mapping[str, _LocalThreshold],
    protein_map: Mapping[str, CdsProtein],
    local_support_by_member: Mapping[str, float],
) -> bool:
    ...
```

Use existing `_build_core_support_candidate()` for each unassigned member
against existing cross-record groups, but use it as an evidence summarizer, not
as the final comparison score.

For each member/core candidate:

- `core_support = max(candidate.cross_record_score, candidate.same_record_score)`.
- `member_local_support = local_support_by_member[member_id]`.

Reject or leave unassigned when:

- a high-confidence cross-record assignment is available;
- low-confidence support to an existing core is within
  `_RECORD_LOCAL_CORE_COMPETITION_RATIO` of the member's local support, i.e.
  `core_support >= member_local_support / _RECORD_LOCAL_CORE_COMPETITION_RATIO`;
- credible support would connect different members of the candidate component
  to two or more existing cross-record cores;
- the strongest evidence is domain-only.

This guard keeps the current "same-record evidence cannot merge cores" rule.

### Phase 8: Integrate Accepted Components

Add helper:

```python
def _append_record_local_orthogroups(
    *,
    group_member_ids: dict[str, set[str]],
    group_by_protein: dict[str, str],
    representative_ids_by_group: dict[str, set[str]],
    member_roles: dict[str, OrthogroupMemberRole],
    member_confidence: dict[str, OrthogroupMemberConfidence],
    assignment_reasons: dict[str, str],
    supporting_edges: dict[str, tuple[str, ...]],
    best_core_support: dict[str, float],
    second_core_support: dict[str, float],
    ortholog_edges_by_group: dict[str, list[OrthologEdge]],
    scope_by_group: dict[str, OrthogroupScope],
    source_record_index_by_group: dict[str, int],
    ...
) -> None:
    ...
```

Integration rules:

- Continue the existing `og_N` sequence after cross-record groups.
- Set `scope_by_group[group_id] = "record_local"`.
- Set `source_record_index_by_group[group_id]` to the component record index.
- Set member role to `local_paralog`.
- Set confidence to `high` when every accepted support edge passed threshold
  and coverage rules.
- Assignment reason: `record-local reciprocal paralog cluster`.
- Choose one representative for the source record using the same rank strategy
  used for cross-record groups.
- Add accepted same-record evidence edges to
  `ortholog_edges_by_group[group_id]` with:
  - `edge_kind="record_local_paralog"`.
  - `render_role="display_edge"`.
- Set `rbh_orthogroups[group_id]` to the sorted member IDs for compatibility
  with existing representative/path code.

Add a dedicated edge kind:

```python
OrthologEdgeKind = Literal[
    "rbh",
    "coortholog",
    "same_record_inparalog",
    "record_local_paralog",
    "related_homolog",
    "ambiguous_paralog",
    "weak_bridge",
    "domain_only",
    "related_paralog",
]
```

Update `ORTHOLOG_DISPLAY_EDGE_KIND_RANK` with
`record_local_paralog` adjacent to `same_record_inparalog`. The exact rank only
affects display ordering because record-local edges are excluded from
collinearity anchor creation in a later phase.

### Phase 9: Final Related Edge Classification

Run the existing related-edge pass after record-local groups have been added.

Before adding related edges, build membership pair keys without `edge_kind`:

```python
membership_pairs = {
    tuple(sorted((edge.query_protein_id, edge.subject_protein_id)))
    for edges in ortholog_edges_by_group.values()
    for edge in edges
}
```

Skip any candidate related edge whose unordered protein pair is already in
`membership_pairs`. This prevents a `record_local_paralog` membership edge from
being re-added as `related_homolog` or another related edge kind.

Expected behavior:

- Same-record edges within a record-local group are membership edges, not
  related edges.
- Edges from a record-local group to a cross-record group are related evidence
  unless both endpoints already share membership, which they should not.
- Cross-group record-local to cross-record edges classify as `weak_bridge`,
  `related_homolog`, `ambiguous_paralog`, or `domain_only` using existing rules.
- Record-local groups do not create collinearity block anchors because block
  anchors require cross-record display pairs.

### Phase 10: Keep Record-Local Edges Out Of Collinearity

Record-local edges should remain available in `OrthogroupResult` for metadata,
inspection, popups, paths, and future UI work. They should not become
collinearity anchors.

Implement this at `_build_adjacent_display_edges_by_pair()`:

```python
if orthogroups.scope_by_orthogroup_id.get(orthogroup_id, "cross_record") == "record_local":
    continue
```

Apply the same exclusion to any direct-candidate path that uses
`_display_pair_orthogroup_id()` before adding rows to adjacent display edge
tables. This keeps the rule centralized at the orthogroup-display to
collinearity boundary.

## Web And Serialization

Update `_serialize_orthogroups_payload()` in
`gbdraw/web/js/app/python-helpers.js` to include:

```javascript
{
  id: "og_12",
  scope: "record_local",
  source_record_index: 0,
  member_count: 2,
  record_coverage_count: 1,
  ...
}
```

Member rows already serialize `role`, `confidence`, and `assignmentReason`.
They should show `role: "local_paralog"` for record-local groups.

Use the existing snake_case payload style. Do not introduce
`sourceRecordIndex` or `recordCoverageCount` unless the whole orthogroup state
pipeline is deliberately migrated to accept both styles.

Update `applyOrthogroupState()` only if needed. It currently stores group
objects and builds the feature index from members, so unknown group fields
should survive. UI helpers should interpret missing `scope` as `cross_record`
for older saved sessions.

UI polish:

- In `gbdraw/web/js/app/orthogroups.js`, display record-local groups as
  `Species-specific orthogroup` or `Record-local orthogroup` where group summary
  text is built.
- Do not add a new filter in the first implementation unless the panel becomes
  confusing.
- Saved sessions without `scope` should behave as existing cross-record groups.

## Test Plan

Add focused tests near the existing `anchor_core_v1` tests in
`tests/test_collinearity.py`.

### Python Unit Tests

1. `test_anchor_core_record_local_paralog_cluster_forms_species_specific_orthogroup`
   - One record with `a0` and `a1`.
   - Reciprocal same-record hits pass coverage and score thresholds.
   - No cross-record hits.
   - Expected: one orthogroup with `{a0, a1}`.
   - Expected group scope: `record_local`.
   - Expected member roles: `local_paralog`.
   - Expected record coverage: 1.

2. `test_anchor_core_record_local_singleton_remains_unassigned`
   - One record with one protein or with no accepted same-record edge.
   - Expected: no orthogroup membership for the singleton.
   - Expected: this remains true even when the caller's seed path previously
     requested singleton-inclusive groups.

3. `test_anchor_core_record_local_groups_do_not_merge_cross_record_cores`
   - Keep the current DnaE/PolC-style guard test behavior.
   - Add a record-local candidate edge that would be tempting to bridge.
   - Expected: cross-record cores stay separate.

4. `test_anchor_core_record_local_ignores_already_assigned_members`
   - `a0-b0` forms a cross-record core.
   - `a0-a1` same-record evidence exists.
   - Existing inparalog assignment should handle `a1` if threshold passes.
   - No separate record-local group should include `a0`.

5. `test_anchor_core_record_local_rejects_domain_only_component`
   - Same-record reciprocal hits have low coverage.
   - Expected: no record-local group.
   - Optional: related edge classified as `domain_only` if a related group
     context exists.

6. `test_anchor_core_record_local_ids_are_deterministic_after_cross_record_groups`
   - Cross-record `og_1` exists.
   - One record-local group is added.
   - Expected: local group ID is the next stable `og_N`.

7. `test_anchor_core_record_local_member_hit_limit_does_not_change_membership`
   - Run with different `orthogroup_member_max_hits` values.
   - Expected: record-local membership is unchanged.

8. `test_anchor_core_record_local_membership_edge_not_duplicated_as_related`
   - A record-local component has an accepted membership edge.
   - Expected: the same unordered protein pair is not also emitted as a related
     edge with a different `edge_kind`.

9. `test_anchor_core_record_local_edges_do_not_create_collinearity_anchors`
   - A record-local group has `record_local_paralog` display edges.
   - Expected: adjacent display edge tables and lossless collinearity anchors
     exclude those edges.

### Web Packaging Tests

Add or extend tests in `tests/test_web_packaging.py`:

- `_serialize_orthogroups_payload()` includes `scope`.
- `_serialize_orthogroups_payload()` includes `source_record_index`, not
  `sourceRecordIndex`.
- member rows can carry `role: "local_paralog"`.
- saved session restore preserves `scope`.
- no CSP or build step changes are needed.

### Suggested Commands

```bash
python -m pytest tests/test_collinearity.py -k "anchor_core or record_local or species_specific" -q
python -m pytest tests/test_web_packaging.py -k "orthogroup" -q
python -m pytest tests/test_collinearity.py tests/test_web_packaging.py -q
```

If SVG output changes only because orthogroup metadata changes, inspect the
affected reference outputs before updating them.

## Implementation Checklist

1. Extend type aliases and dataclasses in `protein_colinearity.py`.
2. Update result copy/build helpers to carry group scope metadata.
3. Add `record_local_paralog` to edge-kind typing and display ranking.
4. Add record-local helper functions near the anchor-core helpers:
   edge selection, component building, local support scoring, core competition
   checks, and component integration.
5. Call the record-local integration helper inside
   `_build_anchor_core_orthogroups()` after existing inparalog assignment and
   before final related-edge classification.
6. Change final related-edge duplicate suppression to use unordered protein
   membership pairs without `edge_kind`.
7. Exclude `scope="record_local"` groups when building adjacent display edges
   for collinearity.
8. Preserve existing behavior when no record-local component is accepted.
9. Add Python unit tests.
10. Add web serialization fields and tests.
11. Update UI display text only where group summary labels already exist.
12. Run focused tests.
13. Run broader orthogroup/collinearity tests.

## Acceptance Criteria

- Strong same-record paralog clusters with no cross-record membership become
  orthogroups with `scope="record_local"`.
- Singletons stay unassigned.
- Singleton-inclusive compatibility paths do not cause isolated proteins to
  become record-local members.
- Existing same-record inparalog assignment to cross-record cores still works.
- Same-record evidence still cannot create or merge cross-record cores.
- Record-local membership edges are not duplicated as related edges.
- Record-local groups do not create collinearity anchors or adjacent display
  rows used by collinearity.
- Orthogroup IDs are deterministic.
- Web payloads expose enough metadata to label record-local groups clearly.
- Web payloads use existing snake_case orthogroup keys, including
  `source_record_index`.
- Existing saved sessions without scope metadata continue to load.

## Risks And Mitigations

### Risk: Broad domain families become record-local orthogroups

Mitigation:

- require reciprocal or near-reciprocal same-record support;
- require membership coverage;
- reject domain-only rows;
- keep singletons unassigned;
- add regression tests for low-coverage domain-only hits.

### Risk: Record-local groups hide a valid cross-record assignment

Mitigation:

- run cross-record core inference and existing inparalog assignment first;
- consider only still-unassigned proteins;
- score record-local support directly from accepted local edges;
- reject components with competing support to existing cores using the direct
  local-vs-core support comparison.

### Risk: UI users confuse record-local groups with cross-record orthogroups

Mitigation:

- serialize `scope`;
- label record-local groups as `Species-specific orthogroup` or
  `Record-local orthogroup`;
- keep member role `local_paralog`.

### Risk: The implementation grows too large inside `protein_colinearity.py`

Mitigation:

- keep helpers small and pure where possible;
- if the record-local section becomes hard to read, move only the new helper
  functions to a dedicated internal module in the same PR, while keeping the
  public API unchanged.

## Notes For Future Refinement

- If real datasets show large chained local components, add an internal split
  step based on edge support gaps before accepting components.
- If users need to hide these groups, add a display-only filter in the
  orthogroup panel rather than changing inference.
- If a future tree-based path is added, `record_local` groups can be treated as
  preliminary groups and refined by gene-tree evidence.
