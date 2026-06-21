# Orthogroup Inference Redesign Plan

## Status

Planned.

This plan replaces the current `rbh`, `family_merge`, and
`distribution_split` orthogroup membership behaviors with one inference model.
The replacement is tree-free and browser-compatible, but it borrows one useful
idea from OrthoFinder's graph-based path:

> Use cross-record ortholog evidence to define the local scale of similarity,
> then allow same-record paralogs only when they are at least as convincing as
> the locally accepted ortholog evidence for an existing core.

This is not an implementation of OrthoFinder HOGs. OrthoFinder's current
recommended HOG output depends on rooted gene trees and a species tree, which is
outside the browser path and outside the scope of this feature.

## Background

The current LOSATP orthogroup modes solve different failure cases by changing
which evidence edges are allowed into the orthogroup graph:

- `rbh` is conservative, but misses many legitimate paralogous members.
- `family_merge` captures broad homologous families, but can merge old
  outparalogs such as DnaE and PolC into one orthogroup.
- `distribution_split` avoids many broad-family merges by relying on
  cross-record evidence and by excluding same-record paralog edges from
  orthogroup formation. This protects DnaE/PolC-like cases, but it also drops
  record-internal inparalogs such as `BDV02433.1` and `BDV02434.1` when
  `BDV02435.1` is the cross-record orthogroup anchor.

The conflict is not solved by exposing more user modes. The model needs to
separate these concepts:

- evidence that creates or merges cross-record orthogroup cores;
- evidence that assigns a protein to an already established core;
- related homologous evidence that is useful for display but should not affect
  membership.

## Goals

- Preserve DnaE/PolC-like old outparalogs as separate orthogroups.
- Assign record-internal inparalogs to the appropriate orthogroup when evidence
  passes local thresholds and is unambiguous.
- Keep weak homologous or domain-only relationships visible as related evidence
  without letting them collapse orthogroups.
- Remove the user-facing choice between `RBH`, `Family merge`, and
  `Distribution split`.
- Make every membership decision deterministic and explainable.

## Non-Goals

- Do not infer full gene trees or require a species tree in the browser path.
- Do not add an external MCL binary dependency.
- Do not use product names such as `DnaE` or `PolC` as primary clustering
  evidence.
- Do not use a fixed top-N hit cap as an inference rule. Hit caps may remain for
  display and diagnostics only.
- Do not preserve old saved-session semantics. If necessary, bump the session
  schema and recompute orthogroups under the new model.

## Design Decisions

### 1. Use Local Thresholds, Not Undefined Strength Labels

The algorithm must not depend on unimplemented words such as "strong",
"near-reciprocal", or "clearly better". Those terms are converted into explicit
internal parameters and deterministic comparisons.

Initial internal defaults:

```text
orthogroup_inference_version = anchor_core_v1
near_reciprocal_min_ratio = 0.85
core_bridge_min_ratio = 0.85
inparalog_min_same_record_ratio_to_local_anchor = 0.75
min_best_second_core_ratio = 1.25
min_membership_min_coverage = 0.30
domain_only_max_min_coverage = 0.25
low_confidence_min_best_second_core_ratio = 1.10
```

These values are intentionally internal at first. They can be adjusted during
regression testing, but the algorithm and tests must refer to named parameters
rather than prose-only thresholds.

Tie-breaking order for ranked hits:

1. Higher `normalized_score`.
2. Lower `evalue`.
3. Higher `min_coverage`.
4. Higher `identity`.
5. Longer `alignment_length`.
6. Stable `(query_record_index, query_protein_id, subject_record_index,
   subject_protein_id)` order.

### 2. Cross-Record Evidence Builds Cores

Cross-record reciprocal or near-reciprocal normalized hits create the initial
orthogroup cores. Same-record hits do not create new cores and do not merge two
existing cores.

This preserves the most important protection from `distribution_split`: a
same-record bridge between two old paralogs cannot collapse two independently
supported cross-record cores.

### 3. Same-Record Evidence Can Assign To A Core

Same-record non-self hits are retained and normalized. They may assign a protein
to an existing core as an inparalog when the candidate has support above local
thresholds for one core and clear separation from competing cores.

Same-record evidence has one permitted membership action:

```text
unassigned protein -> existing core
```

It never has these actions:

```text
new core creation
core A -> core B merge
```

### 4. Membership Roles And Pairwise Relations Are Separate

Orthogroup membership answers "which group is this protein in?" Pairwise
relation answers "what is the relationship between these two proteins?" The two
must be stored separately.

Membership roles:

- `anchor`: cross-record reciprocal or near-reciprocal ortholog evidence.
- `coortholog`: cross-record one-to-many or many-to-many supported core member.
- `inparalog`: same-record paralog assigned to an existing core.
- `low_confidence`: weaker but still unambiguous member.

Pairwise relation kinds:

- `one_to_one_ortholog`
- `one_to_many_ortholog`
- `many_to_many_ortholog`
- `same_record_inparalog`
- `related_homolog`
- `ambiguous_paralog`
- `weak_bridge`
- `domain_only`

Feature coloring and orthogroup popups use membership. Ribbons and diagnostic
views may also use related pairwise evidence.

## User-Facing Model

Remove the orthogroup membership mode selector from the web UI:

```text
RBH
Family merge
Distribution split
```

Replace it with one behavior:

```text
Orthogroups
```

The old `Member hits` setting must not limit inference. It can remain as a
display setting, for example:

```text
Shown related hits per feature
```

The orthogroup popup should expose:

- orthogroup ID;
- member count;
- record coverage;
- member role;
- confidence;
- assignment reason;
- key supporting evidence.

Related but non-member evidence should be inspectable separately so users can see
why a homologous hit did not become membership.

## Core Algorithm

### Phase 1: Collect Protein Hits

Input:

- Directional LOSATP/BLASTP outfmt 6 tables.
- Protein metadata: stable protein ID, source accession, record index, feature
  span, strand, sequence length, qualifiers.

Rules:

1. Exclude self-self hits where query and subject are the same protein.
2. Keep same-record non-self hits.
3. Keep cross-record hits.
4. Apply only basic validity filters:
   - known query and subject proteins;
   - positive protein lengths;
   - positive bitscore;
   - valid alignment coordinates when present.
5. For membership inference, keep the strongest directional hit per
   `(query_protein_id, subject_protein_id)` after applying the deterministic
   ranking order. Additional HSPs may be retained for diagnostics only.

### Phase 2: Normalize Hit Evidence

Each retained directional hit gets derived fields:

- `query_length`
- `subject_length`
- `query_coverage`
- `subject_coverage`
- `min_coverage`
- `length_product`
- `normalized_score`
- `score_rank_for_query`
- `record_pair_normalization_source`
- `domain_only`

Normalization is record-pair aware. The pair `(record_i, record_i)` is included
for same-record non-self hits.

Rules:

1. For each directional record pair, fit expected bitscore as a function of
   `query_length * subject_length` when enough hits are available:

   ```text
   log10(bitscore) = a * log10(length_product) + b
   normalized_score = bitscore / expected_bitscore
   ```

2. If the pair has too few hits for stable fitting, use:

   ```text
   normalized_score = bitscore / sqrt(length_product)
   ```

3. Mark the normalization source as `fit`, `sqrt_length_fallback`, or
   `empty_pair`.
4. Mark domain-only evidence when `min_coverage <= domain_only_max_min_coverage`.
5. Do not drop low-coverage hits at this phase. They can still explain related
   evidence even when they do not support membership.

### Phase 3: Compute Best-Hit Profiles And Local Thresholds

For each protein `q`, compute:

- best normalized hit per target record;
- reciprocal best normalized hits;
- near-reciprocal hits;
- cross-record anchor scores;
- same-record hit scores;
- local membership threshold.

Near-reciprocal definition:

A directional hit `q -> s` is near-reciprocal when:

1. `s` is in a different record from `q`;
2. `score(q, s) >= near_reciprocal_min_ratio * best_score(q, record(s))`;
3. `score(s, q) >= near_reciprocal_min_ratio * best_score(s, record(q))`;
4. both directions pass membership coverage and are not domain-only.

Local threshold for `q`:

1. If `q` has cross-record reciprocal or near-reciprocal anchors, use the
   weakest accepted anchor score as the local anchor floor.
2. If no anchor exists but cross-record hits exist, use the best cross-record hit
   as a low-confidence fallback floor.
3. If no cross-record hits exist, no membership core can be created for `q`;
   same-record hits may be emitted only as related evidence unless they assign
   `q` to another existing core through that core's members.

The local threshold is:

```text
threshold(q) = local_anchor_floor(q) * inparalog_min_same_record_ratio_to_local_anchor
```

When the threshold source is a fallback, any assignment from it must be marked
`low_confidence`.

### Phase 4: Build Cross-Record Cores

Use only cross-record anchor edges to create initial cores.

Core edge types:

- `rbh_anchor`: reciprocal best normalized hit.
- `near_rbh_anchor`: near-reciprocal normalized hit.
- `coortholog_anchor`: one-to-many or many-to-many cross-record evidence where
  multiple proteins independently pass local anchor criteria against the same
  cross-record core.

Rules:

1. Union proteins connected by `rbh_anchor` and `near_rbh_anchor` edges.
2. Add `coortholog_anchor` members when the support is cross-record and passes
   local thresholds.
3. Same-record hits are not considered in union operations.
4. A cross-core bridge can merge two cores only when it is cross-record,
   reciprocal or near-reciprocal, passes local thresholds on both endpoints, and
   each direction is at least `core_bridge_min_ratio` of the endpoint's best
   cross-record hit to the other core's record.
5. A single weak bridge is classified as `weak_bridge` or `related_homolog`, not
   as membership.

This phase is the only phase that can create or merge cores.

### Phase 5: Assign Record-Internal Inparalogs

After core clusters are stable, evaluate unassigned proteins against each core.

For candidate protein `p` and core `C`, compute:

- best same-record hit from `p` to a member of `C`;
- best cross-record hit from `p` to each record represented in `C`;
- profile similarity between `p`'s per-record hit profile and the core's member
  profile;
- coverage and domain-only penalties;
- score gap between the best-supported core and the second-best-supported core.

The core support score should be deterministic and inspectable. The first
implementation can use a simple weighted score:

```text
core_support(p, C) =
  max_cross_record_support(p, C)
  + 0.5 * best_same_record_support(p, C)
  + 0.5 * profile_similarity(p, C)
```

All terms must be normalized to the same score scale before combining. If this
weighted formula performs poorly in regression tests, replace it behind the same
interface rather than spreading formula-specific logic through the codebase.

Assign `p` to the best core as `inparalog` when:

1. `p` is not already an anchor member of a different core;
2. `p` has either:
   - same-record support to a member of the core above `threshold(p)` or the
     core member's threshold; or
   - cross-record profile support above the local fallback threshold;
3. `best_core_support / second_best_core_support >= min_best_second_core_ratio`;
4. the supporting evidence passes `min_membership_min_coverage`;
5. the relationship is not explained only by domain-only hits;
6. the assignment would not require merging two existing cores.

Assign `p` as `low_confidence` only when:

1. the best core is still separated from alternatives by
   `low_confidence_min_best_second_core_ratio`;
2. evidence is not domain-only;
3. no competing established core has comparable support;
4. the UI can display the low-confidence reason.

Leave `p` unassigned or mark it as related evidence when:

- two or more cores receive similar support;
- support comes only from low-coverage domain hits;
- the candidate would bridge DnaE/PolC-like distinct cores;
- there is no cross-record context and no existing core member to assign to.

This is the phase expected to place `BDV02433.1` and `BDV02434.1` with the
`BDV02435.1` core when their same-record and profile evidence are unambiguous.

### Phase 6: Classify Related Edges

After membership is fixed, classify all retained non-member evidence.

Related edge classification:

- `same_record_inparalog`: same-record relationship supporting an assigned
  inparalog.
- `related_homolog`: valid homologous hit that does not affect membership.
- `ambiguous_paralog`: candidate has comparable support for multiple cores.
- `weak_bridge`: cross-core bridge below merge criteria.
- `domain_only`: low-coverage/domain-only relationship.

Related edge classification must never mutate orthogroup membership.

### Phase 7: Emit Deterministic Results

Output should separate:

1. Orthogroup membership:
   - orthogroup ID;
   - member protein IDs;
   - member roles;
   - representative/core members;
   - confidence;
   - assignment reason;
   - supporting evidence references.
2. Pairwise relations:
   - query protein;
   - subject protein;
   - relation kind;
   - normalized score and coverage;
   - reason it did or did not affect membership.
3. Diagnostics:
   - normalization source;
   - local threshold source;
   - best and second-best core support;
   - score-gap ratio;
   - domain-only flag.

Orthogroup IDs must be deterministic. Use sorted member coordinates and stable
protein IDs, not construction order from hash maps.

## Proposed Data Structures

Add explicit internal models, even if implemented as dictionaries at first:

```text
ProteinHitEvidence
  query_protein_id
  subject_protein_id
  query_record_index
  subject_record_index
  bitscore
  evalue
  identity
  alignment_length
  query_length
  subject_length
  query_coverage
  subject_coverage
  min_coverage
  normalized_score
  score_rank_for_query
  record_pair_normalization_source
  domain_only

LocalProteinThreshold
  protein_id
  threshold
  local_anchor_floor
  source
  supporting_anchor_edges

OrthogroupCore
  core_id
  anchor_edges
  member_ids
  represented_records
  confidence

OrthogroupMember
  protein_id
  orthogroup_id
  role
  confidence
  assignment_reason
  supporting_edges
  best_core_support
  second_best_core_support

PairwiseProteinRelation
  query_protein_id
  subject_protein_id
  relation_kind
  normalized_score
  min_coverage
  membership_effect
  reason
```

Avoid using the same edge object for membership construction and display.

## Implementation Touchpoints

Primary Python path:

- `gbdraw/analysis/protein_colinearity.py`
  - remove or bypass `ORTHOGROUP_MEMBERSHIP_MODES`;
  - replace `select_rbh_orthogroup_edges_from_directional_hits()`;
  - replace `build_orthogroups_from_protein_hits()`;
  - replace mode-specific membership expansion.

Recommended extraction:

- `gbdraw/analysis/orthogroups.py`
  - hit normalization;
  - best-hit ranking;
  - local threshold derivation;
  - anchor selection;
  - core construction;
  - inparalog assignment;
  - related-edge classification.

Web/Pyodide path:

- `gbdraw/web/js/state.js`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/app/losat-settings.js`
- `gbdraw/web/js/app/python-helpers.js`
- `gbdraw/web/js/app/pairwise-match-popup.js`

UI tasks:

- remove the mode selector;
- remove or rename inference-limiting `Member hits`;
- expose member role/confidence in the orthogroup popup;
- expose related-only evidence separately from membership;
- keep saved session schema explicit so stale sessions do not silently imply old
  mode semantics.

## Migration Strategy

The old behaviors are intentionally replaced.

1. Remove old mode names from the web UI.
2. Keep parsing of old config/session names long enough to load existing data
   without crashing.
3. Normalize old names to the new inference behavior with an explicit warning in
   CLI/API contexts.
4. Bump any saved-session schema version that stores LOSATP settings.
5. Invalidate cached LOSATP orthogroup results when the new inference version
   changes.
6. Include `orthogroup_inference_version` in cache keys and result metadata.

The web UI should not silently resurrect old mode semantics from saved sessions.

## Development Phases

### Phase 1: Characterize Current Failures

- Add failing tests for the BDV02433/BDV02434/BDV02435 inparalog case.
- Add failing tests for the DnaE/PolC outparalog split case.
- Add a many-to-one or one-to-many coortholog fixture.
- Capture diagnostic expectations before changing the algorithm.

### Phase 2: Extract Orthogroup Inference

- Move mode-specific logic out of the main colinearity flow.
- Introduce a single internal inference entry point.
- Preserve current output shape temporarily so rendering does not change during
  extraction.

### Phase 3: Implement Evidence Normalization

- Collect same-record non-self hits and cross-record hits through one evidence
  path.
- Compute normalized scores and coverage once.
- Normalize `(record, record)` same-record pairs explicitly.
- Store normalization diagnostics for every retained hit.
- Remove inference dependence on `orthogroup_member_max_hits`.

### Phase 4: Implement Local Thresholds

- Rank normalized hits deterministically.
- Select reciprocal and near-reciprocal cross-record anchors.
- Derive per-protein local thresholds.
- Add diagnostics for threshold source and accepted edge count.

### Phase 5: Implement Anchor-Constrained Core Building

- Build initial cores from cross-record anchors only.
- Add coortholog core members from cross-record evidence.
- Prevent same-record edge-driven core creation and core merging.
- Classify cross-core bridges as related evidence.

### Phase 6: Implement Inparalog Assignment

- Score unassigned candidates against stable cores.
- Assign unambiguous same-record inparalogs.
- Mark ambiguous or domain-only candidates as related-only.
- Emit role, confidence, and reason for each member.

### Phase 7: Emit Pairwise Relations And UI Metadata

- Separate orthogroup membership from pairwise relation output.
- Update standalone SVG metadata and web popup payloads.
- Show member roles and confidence in orthogroup popups.
- Keep related-only evidence inspectable without coloring it as membership.

### Phase 8: Regression And Performance Pass

- Run protein colinearity tests and web packaging tests.
- Benchmark the Hepatoplasmataceae session.
- Confirm that browser runtime remains acceptable.
- Add diagnostics for any slow path before optimizing.

## Testing Plan

### Unit Tests

Add focused tests for:

- self-self hits are excluded;
- same-record non-self hits are retained;
- same-record pairs are normalized through the `(record, record)` path;
- sparse record-pair fallback normalization is deterministic;
- reciprocal and near-reciprocal anchors are selected deterministically;
- local thresholds are derived from anchor scores;
- same-record hits cannot create a core;
- same-record hits cannot merge two established cross-record cores;
- unambiguous same-record inparalogs are assigned to the correct core;
- ambiguous candidates remain unassigned or related-only;
- low-coverage domain-only hits do not create membership;
- one-to-many or many-to-many pairwise relation multiplicity is represented
  separately from member role;
- hit display caps do not change inferred orthogroups.

### Regression Fixtures

Add or update fixtures for:

- DnaE/PolC outparalogs remain split;
- `BDV02433.1`, `BDV02434.1`, and `BDV02435.1` resolve to the same OG when
  evidence is unambiguous;
- a broad multi-domain family with weak bridges remains separated;
- a many-to-one inparalog case produces one OG with multiple members from one
  record;
- a same-record-only duplicated pair is shown as related evidence rather than a
  newly inferred cross-record orthogroup.

### Web Tests

Verify:

- the old mode selector is gone;
- old session settings do not resurrect old mode behavior;
- orthogroup popups show member role/confidence;
- related-only evidence can be inspected separately from membership;
- SVG export remains deterministic for a fixed input.

## Acceptance Criteria

The redesign is successful when:

- users no longer need to choose `RBH`, `Family merge`, or `Distribution split`;
- DnaE/PolC-like outparalogs stay split without annotation-specific rules;
- BDV02433/BDV02434/BDV02435-like record-internal inparalogs are assigned to one
  OG when evidence is unambiguous;
- same-record evidence can assign proteins to stable cores but cannot create or
  merge cores;
- ambiguous or domain-only relationships are visible as related evidence but do
  not collapse OGs;
- changing display hit limits does not change OG membership;
- every orthogroup member has an explainable role and assignment reason;
- membership roles and pairwise relation kinds are emitted as separate data.
