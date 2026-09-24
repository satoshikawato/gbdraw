# Issue #561 — S01 domain model and resolver instruction prompt

Use the complete code block below as the instruction for a fresh implementation
session. The reader is not expected to know any earlier conversation.

```text
Implement Session S01 of gbdraw Issue #561: create the typed Similarity Group
alignment domain model and the shared pure anchor resolver. Complete focused
tests and update the execution ledger, but do not change rendering, Web UI,
Session/request schemas, CLI behavior, or public documentation in this session.

Repository and branch:
- Work in the gbdraw repository.
- Use the fixed branch issue-561-similarity-alignment.
- Do not create another runtime branch and do not commit to main or dev.
- Fetch origin and verify that the branch is based on a dev commit containing
  accepted PD-OI-026 through PD-OI-031 (or their final IDs for the same six
  concern keys). Candidate authority on this branch is insufficient. If the
  authority is absent from origin/dev, stop runtime work and report that exact
  boundary.

Read before editing:
1. AGENTS.md
2. CLAUDE.md
3. docs/internal/ISSUE_561_SIMILARITY_ALIGNMENT_MASTER_PLAN_2026-09-22.md
4. docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md, including the six
   similarity-alignment decisions
5. docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md
6. docs/internal/PRODUCT_IMPACT_RATCHET.md
7. gbdraw/layout/record_coordinates.py
8. gbdraw/diagrams/linear/orthogroup_alignment.py
9. gbdraw/analysis/protein_colinearity.py and the Orthogroup result/edge types
10. gbdraw/web_support/orthogroup_metadata.py
11. existing orthogroup, stable-identity, and alignment tests

Problem:
The current alignment accepts a group or feature string and lets renderer code
choose representatives or score-ranked members. The accepted behavior requires
the exact reference feature, one decision per displayed record, strict explicit
selection/single usable candidate/unique direct RBH priority, and no geometric,
representative, score, rendering, or multi-hop fallback.

Implement:
1. Add one small immutable domain boundary, preferably
   gbdraw/layout/similarity_alignment.py unless current ownership shows a more
   appropriate existing layout module. Keep it independent of Vue, SVG,
   filesystem I/O, Worker messages, and LOSAT execution.
2. Define validated types for:
   - canonical anchor identity keyed by record_key and biological_feature_id,
     with source_feature_index only for disambiguation and stable feature ID as
     consistency evidence;
   - candidate/member facts including group, displayed strand, mappability, and
     direct evidence;
   - explicit per-record choice;
   - resolution rationale enum;
   - resolved, skipped, and ambiguous record outcomes; and
   - immutable SimilarityAlignmentPlan with schema 1, mode, group, reference,
     and record decisions.
3. Validate uniqueness, identity consistency, exactly one reference, finite and
   supported values, allowed rationale/status combinations, and record-key
   coverage without relying on input order.
4. Implement one pure resolver with this exact priority:
   explicit valid selection; only usable candidate; one distinct direct RBH;
   ambiguous candidates; no usable candidate.
5. Normalize RBH direction. More than one distinct direct-RBH candidate stays
   ambiguous. Non-RBH, representative status, scores, confidence, edge count,
   coordinates, viewport center, ribbon state, and multi-hop paths do not rank.
6. Define usability as group membership plus uniquely resolvable canonical
   identity plus a center mappable into the current crop/display domain. Hidden
   is not unusable. Partial overlap with center outside is unusable.
7. Preserve deterministic output ordering by canonical record order and stable
   feature identity, independent of candidate/edge input order.
8. Keep current implicit selection code temporarily for the later legacy-reader
   adapter, but prevent new S01 types from depending on it. Do not delete or
   redirect current runtime in this session.

Testing:
- exact non-representative reference;
- zero, one, and several candidates;
- explicit choice precedence;
- one direct RBH, reversed query/subject RBH, multiple RBHs;
- non-RBH and multi-hop evidence ignored;
- hidden/mappable and partial-crop/unmappable cases;
- duplicate record key, conflicting aliases, duplicate decisions, stale
  explicit choice, and missing reference rejection;
- input permutation produces identical plan/resolution;
- no import or call to LOSAT execution, UI, or renderer code.

Design constraints:
- SOLID: resolver has one reason to change and returns data instead of causing
  UI or render effects.
- KISS: one priority list and a small closed enum set.
- DRY: the same resolver will be called by Python, CLI, and the Web Worker.
- YAGNI: no synteny, multi-hop, smart alignment, TSV, Collinear UI, general
  graph-ranking framework, or plugin system.

Architecture review:
- Record the before/after semantic owner and canonical path in the master plan
  ledger. This should be a non-increasing change: one shared resolver replaces
  future surface-specific resolvers.
- Do not change the accepted Product outcomes. If stable current identity cannot
  express a required case, report the exact mismatch before adding another
  identity system.

Finish:
- Run focused resolver/model tests and ruff on changed Python files.
- Review production and test diffs separately.
- Update section 16 of the master plan with base/head, files, commands, results,
  risks, and whether S02 may start.
- Provide an English proposed commit title and short summary. Do not push or
  create a PR unless the active task separately authorizes it.
```
