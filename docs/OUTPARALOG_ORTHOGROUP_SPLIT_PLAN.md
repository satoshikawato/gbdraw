# Outparalog Orthogroup Split Plan

## Status

Planned.

## Background

In LOSATP `blastp` orthogroup mode, broad homologous families can be merged into
one connected component. One concrete failure case is Hepatoplasmataceae
`AP027078.1`, where `BDU67452.1` (`DnaE`) and `BDU67697.1` (`PolC`) are assigned
to the same orthogroup. Both are DNA polymerases, but they represent different
genes that had already diverged before the common ancestor of the records in the
session. Under the OrthoFinder orthogroup definition, genes should be grouped by
descent from a single gene in the LCA of the analysed record set; these old
outparalogs should therefore be split into separate orthogroups.

References:

- OrthoFinder README: <https://github.com/davidemms/OrthoFinder#orthogroups-orthologs--paralogs>
- OrthoFinder 2015 score-normalisation paper: <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-015-0721-2>

## Goal

Add a tree-free, MCL-free orthogroup refinement step that prevents old
outparalog families from being merged only because weak homologous bridges exist
between otherwise dense subfamilies.

The first implementation adopts three ideas:

1. Normalize protein-hit similarity scores before clustering.
2. Derive local per-protein inclusion thresholds from reciprocal best normalized
   hits and each protein's own similarity distribution.
3. Recursively split broad connected components when their internal score
   distribution supports stable, dense subclusters separated by weak bridges.

## Non-Goals

- Do not call MCL or add an MCL binary dependency.
- Do not infer or require a rooted species tree.
- Do not claim phylogenetic proof of orthology from similarity alone. The result
  is a conservative graph-based estimate.
- Do not use annotation text such as `DnaE` or `PolC` as the primary split
  criterion. Annotation is useful for diagnostics and display names, but the
  clustering decision should come from hit-score structure.
- Do not redesign collinearity block calling. This plan only changes the
  orthogroup assignment feeding orthogroup and collinear protein comparison
  modes.

## Current Implementation Touchpoints

Primary Python path:

- `gbdraw/analysis/protein_colinearity.py`
  - `ORTHOGROUP_MEMBERSHIP_MODES`
  - `normalize_orthogroup_membership_mode()`
  - `select_rbh_orthogroup_edges_from_directional_hits()`
  - `build_orthogroups_from_protein_hits()`
  - `expand_orthogroup_membership_from_evidence()`

Browser/Pyodide path:

- `gbdraw/web/js/app/python-helpers.js`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/app/losat-settings.js`

Tests:

- `tests/test_protein_colinearity.py`
- `tests/test_collinearity.py`
- `tests/test_web_packaging.py`

## User-Facing Model

Add a new orthogroup membership mode:

```text
distribution_split
```

Aliases:

```text
local_split
density_split
outparalog_split
```

Initial rollout:

- Keep `rbh` unchanged.
- Keep `family_merge` available for backward compatibility.
- Add `distribution_split` as the new recommended mode.
- After regression testing on Hepatoplasmataceae and existing examples, make the
  web default `distribution_split` so new sessions avoid broad-family merges by
  default.

## Core Algorithm

### Phase 1: Normalize Hit Scores

Input:

- Directional LOSATP outfmt 6 tables keyed by `(query_record_index,
  subject_record_index)`.
- `protein_map`, providing protein lengths and record membership.

Tasks:

1. Add internal normalized-score columns without changing emitted comparison TSV
   columns:
   - `query_length`
   - `subject_length`
   - `length_product`
   - `query_coverage`
   - `subject_coverage`
   - `min_coverage`
   - `normalized_score`
2. For each directional record pair, fit a length-aware expected-best-score
   model:
   - sort hits by `query_length * subject_length`;
   - bin hits by rank, using fixed-size bins where possible;
   - select the top score fraction in each bin;
   - fit `log10(bitscore) = a * log10(length_product) + b`;
   - compute `normalized_score = bitscore / expected_bitscore`.
3. If a record pair has too few hits to fit a stable model, fall back to a
   deterministic ratio:
   - `normalized_score = bitscore / sqrt(query_length * subject_length)`;
   - mark the pair as fallback-normalized for diagnostics.
4. Clamp or drop pathological rows:
   - missing protein IDs;
   - nonpositive protein lengths;
   - nonpositive bitscores;
   - very low `min_coverage` after the existing visible-hit filters.

Rationale:

Raw bit score and e-value are length-biased. DnaE/PolC-like superfamily bridges
often survive broad filters because long proteins can produce strong-looking raw
scores even when they are not the same orthogroup. Normalization should make
within-DnaE, within-PolC, and DnaE-PolC scores comparable on the same scale.

### Phase 2: Build Local Inclusion Thresholds

Definitions:

- RBNH: reciprocal best normalized hit. This is the existing RBH idea, but the
  best-hit ranking uses `normalized_score`, then e-value, identity, alignment
  length, and stable protein IDs as deterministic tie-breakers.
- Local threshold: the minimum normalized score at which a protein accepts
  another protein as same-orthogroup evidence.

Tasks:

1. Select RBNH anchors from normalized directional tables.
2. For each protein `q`, collect:
   - normalized scores to all candidate hits;
   - normalized scores to RBNH anchors;
   - top-hit score gaps in descending score order.
3. Derive `threshold[q]`:
   - if `q` has RBNHs, start from the lower quantile or minimum RBNH score;
   - if a clear score valley exists above the RBNH floor, raise the threshold to
     the valley;
   - if `q` has no RBNH, use a conservative top-hit-relative fallback and mark
     it low confidence.
4. Keep an edge `q-h` as orthogroup evidence only when it passes local
   thresholds:
   - default policy: both endpoints pass their threshold when both thresholds
     are known;
   - fallback policy: one endpoint can accept only when the other endpoint has
     no usable threshold, and the edge remains non-anchor evidence.
5. Preserve RBNH anchors even when the local threshold would otherwise drop them.

Output:

- A weighted graph of local-threshold orthogroup evidence edges.
- Diagnostics keyed by protein:
  - threshold value;
  - threshold source (`rbnh`, `valley`, `fallback`);
  - RBNH count;
  - accepted edge count.

Rationale:

One global score cutoff cannot distinguish old paralogs from true orthogroup
members across records of different divergence and protein lengths. A DnaE node
should judge candidate members relative to its own DnaE-like best hits, not
relative to a session-wide threshold.

### Phase 3: Recursively Split Broad Components

Input:

- Existing broad components from RBNH/family evidence.
- The local-threshold weighted graph from Phase 2.

Tasks:

1. For each broad component, create an induced weighted graph using
   `normalized_score` as edge weight.
2. If the local-threshold graph already separates the component into multiple
   connected subcomponents, treat those subcomponents as candidate splits.
3. If it remains connected, run a deterministic threshold sweep:
   - sort unique edge weights from high to low;
   - union edges above each candidate threshold;
   - record partitions that persist across a score interval;
   - prefer partitions with high within-cluster density and low cross-cluster
     bridge strength.
4. Score each candidate split:
   - `within_median_min`: minimum median edge weight among child clusters;
   - `cross_p95`: 95th percentile of cross-cluster edge weights;
   - `separation_ratio = within_median_min / max(cross_p95, epsilon)`;
   - `conductance = cross_weight / (cross_weight + internal_weight)`;
   - `stability`: number or width of adjacent thresholds giving the same
     partition.
5. Accept a split only if all conservative criteria pass:
   - at least two child clusters have size >= 2, unless one child has strong
     multi-record RBNH support;
   - `separation_ratio` exceeds the configured minimum;
   - `conductance` is below the configured maximum;
   - the partition is stable across at least two adjacent threshold levels or a
     minimum score interval;
   - no child cluster is made only of weak fallback edges.
6. Recurse on accepted child clusters until no further split passes.
7. Renumber final orthogroups deterministically by genomic/protein sort order.

MCL-like behavior without MCL:

- Dense groups are reinforced by local-threshold mutual support.
- Weak bridge edges lose influence because they fail local thresholds or only
  appear at low threshold-sweep levels.
- Recursive splitting gives an inflation-like effect: broad components are
  repeatedly sharpened into denser subcomponents, but the implementation remains
  transparent and deterministic.

## Data Structures

Add private dataclasses in `gbdraw/analysis/protein_colinearity.py`:

```python
@dataclass(frozen=True)
class _NormalizedProteinHit:
    query_id: str
    subject_id: str
    query_record_index: int
    subject_record_index: int
    bitscore: float
    evalue: float
    identity: float
    alignment_length: int
    query_length: int
    subject_length: int
    query_coverage: float
    subject_coverage: float
    normalized_score: float


@dataclass(frozen=True)
class _LocalThreshold:
    protein_id: str
    score: float
    source: str
    rbnh_count: int
    accepted_edge_count: int


@dataclass(frozen=True)
class _SplitDecision:
    accepted: bool
    child_components: tuple[tuple[str, ...], ...]
    threshold: float | None
    separation_ratio: float
    conductance: float
    stability: int
    reason: str
```

Keep these private until the behavior stabilizes. Public API surface should stay
limited to the new membership mode and configuration parameters.

## Configuration Parameters

Mode-specific defaults should be conservative. These values apply when
`distribution_split` is selected; changing the CLI/web default to this mode is a
separate rollout decision after regression testing.

```text
orthogroup_split_min_separation_ratio = 1.35
orthogroup_split_max_conductance = 0.20
orthogroup_split_min_stability = 2
orthogroup_split_min_child_size = 2
orthogroup_split_top_fraction = 0.05
orthogroup_split_min_coverage = 0.30
```

Expose only the mode in the first web UI pass. Keep numeric tuning parameters
internal or CLI/API-only until there is enough real-data feedback to justify UI
controls.

## Implementation Phases

### Phase 1: Normalization Helpers

Files:

- `gbdraw/analysis/protein_colinearity.py`
- `tests/test_protein_colinearity.py`

Tasks:

1. Add helpers to attach protein lengths and coverage to hit rows.
2. Add pairwise length-aware normalized-score fitting.
3. Add deterministic fallback normalization for sparse pairs.
4. Add normalized best-hit ranking.
5. Add focused unit tests for:
   - long-protein bridge downweighting;
   - sparse-pair fallback;
   - deterministic ties.

Acceptance:

- Existing `rbh` and `family_merge` behavior is unchanged unless the new helpers
  are explicitly used.
- Normalized-score helpers do not alter output comparison columns.

### Phase 2: Local Threshold Edge Selection

Files:

- `gbdraw/analysis/protein_colinearity.py`
- `tests/test_protein_colinearity.py`

Tasks:

1. Add RBNH selection using normalized score.
2. Add per-protein threshold derivation.
3. Add local-threshold edge filtering.
4. Add diagnostics in debug logs, not user-facing output.
5. Add synthetic tests where:
   - A1-A2-A3 DnaE-like hits pass;
   - B1-B2-B3 PolC-like hits pass;
   - DnaE-PolC bridge hits fail local thresholds.

Acceptance:

- RBNH anchors are preserved.
- Weak cross-family bridges are dropped even if they pass the broad visible-hit
  filters.
- The algorithm is deterministic across Python versions supported by the
  project.

### Phase 3: Recursive Component Splitting

Files:

- `gbdraw/analysis/protein_colinearity.py`
- `tests/test_protein_colinearity.py`

Tasks:

1. Add recursive split evaluator for broad components.
2. Add split scoring and conservative acceptance criteria.
3. Rebuild `OrthogroupResult` from final split components while preserving:
   - `OrthogroupMember` metadata;
   - representative selection;
   - name/description metadata;
   - ortholog and related display edges where they still connect members of the
     same final orthogroup.
4. Ensure cross-split display edges keep no `orthogroup_id` or are downgraded to
   non-orthogroup pairwise evidence, depending on the existing rendering path.

Acceptance:

- The synthetic DnaE/PolC graph produces two orthogroups.
- No final orthogroup contains both members of an accepted split.
- Orthogroup IDs remain stable for unchanged inputs.

### Phase 4: Mode Plumbing

Files:

- `gbdraw/analysis/protein_colinearity.py`
- `gbdraw/web/js/app/python-helpers.js`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/app/losat-settings.js`
- CLI/API option plumbing files if they currently validate membership modes.
- `tests/test_web_packaging.py`

Tasks:

1. Add `distribution_split` to `ORTHOGROUP_MEMBERSHIP_MODES`.
2. Add aliases in Python and web config normalization.
3. Route `distribution_split` through
   `select_rbh_orthogroup_edges_from_directional_hits()`.
4. Keep saved sessions using `family_merge` loadable.
5. Add web packaging tests confirming the mode is accepted and serialized.

Acceptance:

- Existing sessions using `rbh` or `family_merge` still load.
- New web sessions can request `distribution_split`.
- Pyodide helper and Python package paths produce matching orthogroups for the
  same payload.

### Phase 5: Regression Cases

Files:

- `tests/test_protein_colinearity.py`
- optional curated fixture under `tests/test_inputs/`

Tasks:

1. Add a compact synthetic regression test for a duplicated ancestral family:
   - multiple records;
   - each record has one DnaE-like and one PolC-like protein;
   - within-family scores are high;
   - cross-family scores are nonzero and strong enough to overmerge under
     connected-component family merge.
2. Add a real-data regression when fixture size is acceptable:
   - session including Hepatoplasmataceae `AP027078.1`;
   - assert `BDU67452.1` and `BDU67697.1` receive different final
     `orthogroup_id` values in `distribution_split` mode.
3. Add a guard test that close inparalog/co-ortholog examples already covered by
   existing tests are not over-split.

Acceptance:

- `BDU67452.1` and `BDU67697.1` are not in the same orthogroup.
- Existing orthogroup popup and collinearity tests still pass after expected
  fixture updates.

## Diagnostics

Add debug-level summaries for `distribution_split`:

```text
orthogroup split: broad_component=og_12 members=18 accepted=true children=2 separation=1.92 conductance=0.08
orthogroup threshold: protein=BDU67452.1 source=rbnh score=0.71 accepted_edges=7
```

Do not emit these in normal CLI/web output unless a verbose/debug mode already
exists for the relevant path.

## Risks And Mitigations

- Risk: fragmented proteins can form artificially weak subclusters.
  Mitigation: require minimum coverage and avoid splitting child clusters made
  only from fallback evidence.
- Risk: recent inparalogs can be split if their within-copy scores form two
  dense groups.
  Mitigation: conservative separation/conductance thresholds and regression
  tests for existing co-ortholog examples.
- Risk: sparse sessions lack enough hits for robust distribution fitting.
  Mitigation: deterministic fallback normalization and no split unless stability
  criteria pass.
- Risk: browser helper code diverges from Python package behavior.
  Mitigation: keep `python-helpers.js` embedded Python updated with the same
  helper functions and add web packaging tests that execute the helper.

## Open Decisions

- Whether `distribution_split` should replace `family_merge` as the default in
  CLI as well as the web UI.
- Whether cross-split homologous hits should still be drawn as plain pairwise
  matches or omitted from orthogroup mode output.
- Whether to persist split diagnostics in session JSON for inspectability, or
  keep them debug-only.
