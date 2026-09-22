# Issue #563 baseline and preflight

- Status: `IMPLEMENTED_PENDING_EXACT_HEAD_ARCHITECTURE_REVIEW`
- Session: 07
- Recorded: 2026-09-22
- Branch: `issue-563-feature-popup-record-rotation-20260922`
- Session start HEAD: `3672a518e09339bffc2cc70cdee42add9103ccb2`
- Base: `origin/dev` @ `a9eaeadd105e0c26e46086626feaa695bdd33c94`
- Integration base: `origin/dev` @ `35b2716a33c726cac82842b839383f6ab294b3e6`
- Session 07 start HEAD: `10075e3d0f92ec7a69ab9d5680c2af35a12e5af2`
- PR target base observed after gates: `origin/dev` @ `1cb2a1b00ac97fe8f1213d12a63257ead0d56ded`
- Upstream: none
- Issue: `#563 Add feature-based rotation of circular records from the feature popup`

## Outcome

Sessions 01–07 are complete. The popup workflow, source-coordinate resolver,
per-record state, schema-7 target-only projection, shared candidate admission,
atomic History integration, Session 44/catalog 4 persistence, and browser
acceptance coverage are implemented. Session 07 documentation, generated
artifacts, separate reviews, and full local gates are complete.

The Product stop condition is resolved by `PD-OI-032` as recorded below. The
one remaining pre-merge authority boundary is the architecture exception for
the released Session 42/catalog 3 compatibility reader (`CB 0 -> 1`). Per the
architecture ratchet, that decision must be made manually against the exact
final head after CI; it cannot be supplied by the implementation agent.

## Authority resolution

The Product Decision Owner selected `A / POPUP-RECORD-ROTATION` in a complete
receipt dated `2026-09-22`. That receipt is serialized as `PD-OI-032`, merged
into `dev` by PR `#572` at
`322197a2e792ad9cabf28326bc21b7fafa537fa0`, and present in this branch through
merge commit `fbf195f2`. No conflicting authority appeared. The Product
authority stop condition is therefore cleared for Session 02 runtime work.

## Repository and change baseline

The required isolated worktree was already present at
`/tmp/gbdraw-issue-563-plan-20260922`. Start checks were:

```text
branch: issue-563-feature-popup-record-rotation-20260922
status: clean
HEAD: 3672a518e09339bffc2cc70cdee42add9103ccb2
last commit: 3672a518 docs: plan feature-based record rotation
upstream: none (expected)
```

The branch differs from `origin/dev` only by the eight saved implementation
plan documents. At Session 01 start:

- production diff: none;
- test diff: none;
- documentation diff: the saved Issue #563 implementation plan only;
- generated/ignored state: no tracked generated change. The browser wheel
  prepared for characterization remains ignored.

Git history contains no completed Session result commit and this directory had
no `BASELINE_AND_PREFLIGHT.md` or `FINAL_ACCEPTANCE.md`. Session 01 was therefore
the first incomplete Session.

## Actual implementation behavior

- `RecordPresentation.reverse_complement` and
  `RecordDisplayOptions.start_coordinate` are the typed absolute display
  values. Request schema 7 already carries both without a new Issue-specific
  request field.
- `RecordDisplayTransform` owns source/display base and boundary conversion for
  circular and reverse-complemented records.
- `record-display-options.js` owns the source-bound editable display draft,
  absolute orientation override, popup provenance, and existing sidebar
  shortcuts. Popup rotation uses the explicit open-popup
  `(recordKey, biologicalFeatureId)` and never falls back to global selection.
- `session-request.js` is the single canonical request owner. It materializes a
  `cardinality: all` record into exact-one records when per-record display
  values differ.
- `run-analysis.js` builds the canonical request and exposes the same candidate
  execution, diagram Worker, typed response, catalog validation, and Result
  admission path to normal Generate and popup Apply.
- `history.js::runUndoableArtifactReplacement()` remains the one artifact
  replacement transaction owner and now accepts bounded target-intent
  checkpoints so transform and Result share one Undo/Redo entry.
- The current writer is Session 44, the normalized feature catalog is schema
  4, and the canonical request remains schema 7. The only new compatibility
  reader promotes released Session 42/catalog 3 documents at the existing
  Session admission boundary; branch-only Session 43 is not supported.
- The characterization journey submitted four LOSATP jobs initially. Changing
  display start and then reverse complement left the cumulative count at four;
  a fresh Session Load submitted zero jobs. A target-only candidate is therefore
  needed for pending-edit isolation and atomic History, not to invent a LOSAT
  bypass.

## Product Impact developer preflight (historical Session 01 record)

### Trigger and user effects

The proposal adds a new popup entry point, new editable anchor/offset/end and
orientation choices, a target-only regeneration action, new disabled and stale
states, new Session provenance, and an atomic Result/intent History outcome.
These are material discoverability, state mutation, persistence, failure, and
next-action effects covered by the mandatory developer preflight.

Affected checkpoints are popup open, preview resolution, Apply submission,
fresh Result admission, failure/cancel/stale recovery, Undo/Redo, Session save,
fresh Load, and later regeneration.

### Authority search

| Source | Evidence | Result |
| --- | --- | --- |
| Product Impact map | Four hard concerns cover canonical request and current Result ownership/continuity. | They preserve boundaries but do not select a feature-popup rotation affordance. |
| Durable `BD-###` store | `tools/web-product-decisions.json` has no active decisions. | No authority for Issue #563. |
| Option Integrity Product Contract | Revision 10 is present on the base. OIPC-C06/C07 preserve explicit clearing and failure isolation. PD-OI-026–031 authorize Similarity Group alignment only. | Useful non-waivable continuations; they do not authorize extending that outcome to arbitrary popup rotation. |
| Static Web architecture contract | `gbdraw/web/CLAUDE.md` requires one request owner, live-edit invariants, sanitization, and failure preservation. | Implementation constraint, not selection of the new workflow. |
| Public compatibility/manual docs | Existing record rotation is sidebar-driven and applied by normal Generate. | Evidence of the current supported journey, not authority for the new popup action. |
| Scientific/domain rule | Source coordinates, strand traversal, and non-destructive transforms constrain correct calculation. | They do not choose whether to add the affordance, its controls, persistence, or atomic action. |
| Released compatibility | Main has Session 42 and catalog schema 3. | It constrains migration only. |
| Eligible current decision | None. The implementation request omits the required product rationale, preservation/retirement scope, residual risk, owner, and decision date. | Not an eligible receipt. |
| Issue #563 latest body | Updated 2026-09-21 by `satoshikawato`; contains the proposed semantics and ten scenarios. | Detailed proposal/evidence, but not base-branch durable authority by itself. |

### Classification

`PRODUCT_DECISION_REQUIRED`

Two materially different product-valid outcomes remain: add the complete
Issue #563 workflow, or retain the existing sidebar-only rotation workflow.
The first adds a new public affordance and persisted/History behavior; the
second preserves current behavior. No merged authority chooses between them.

- `IMPLEMENT_EXISTING_AUTHORITY` does not apply because existing authority
  selects only shared invariants and continuation behavior, not the complete
  new outcome.
- `EVIDENCE_REQUIRED` does not apply because the Issue already supplies the
  intended semantics and the baseline verifies technical feasibility; more
  deterministic measurement cannot choose whether the new affordance should
  exist.
- `NOT_ALLOWED` does not apply to the proposed outcome itself. It can be
  implemented within the existing typed request, transform, admission, and
  History owners after authority and architecture review.

The proposed outcome is a `PRODUCT_CHANGE` and requires durable authority on
the runtime implementation base. Candidate authority cannot authorize runtime
in the same candidate. The intended authority route is an authority-only
addition to the existing Option Integrity Product Contract (next available
decision ID, with concern
`diagram-generation.feature-popup-record-rotation`, scenario revision 1), or an
equivalent eligible durable decision selected by the Product Decision Owner.

## Issue and acceptance trace

The ten Issue scenarios map directly to the plan:

| Issue scenario | Contract | Actual implementation evidence |
| --- | --- | --- |
| 1 Popup-only | AC-01, AC-18 | Real rich/simple popup tests cover pointer, keyboard, search continuity, and 390 px layout with an explicit clicked-feature target. |
| 2 Both diagram modes | AC-02, AC-15 | Browser acceptance resolves the same source anchor in Circular and Linear; both use the existing typed transform. |
| 3 Independent chromosomes | AC-03, AC-09 | Same-file multi-record browser acceptance changes one exact record and preserves the other request and transform. |
| 4 Strand-aware offsets | AC-04 | Pure resolver tests cover plus/minus signed offsets and circular wrapping. |
| 5 Stable orientation | AC-05, AC-16 | Resolver tests cover absolute/idempotent orientation; record-state tests cover manual provenance clearing. |
| 6 Endpoint distinction | AC-06 | Pure resolver tests distinguish the 3′ covered base from the outgoing feature-end boundary. |
| 7 Circular/compound | AC-07 | Catalog 4 preserves capability facts; resolver tests cover multipart gaps, origin wrapping, and odd/even midpoint rules. |
| 8 Safe eligibility | AC-08, AC-12, AC-17 | Unit/browser tests cover operation-specific reasons plus cancel, stale, failure, and superseded preservation. |
| 9 Identity/isolation | AC-03, AC-09 | Exact source-bound identity and duplicate/split-fragment projections are covered in catalog, state, and browser tests. |
| 10 State/rendering | AC-10–AC-16, AC-19–AC-20 | Shared request/admission/History owners, Session round trip, zero added LOSAT jobs, and architecture contracts are exercised. |

| AC | Required outcome/evidence | Actual implementation result |
| --- | --- | --- |
| AC-01 | Explicit open-popup target, no global selection fallback | PASS: controller and real popup browser journey use the captured target only. |
| AC-02 | Same source coordinate in Circular and Linear | PASS: cross-mode browser scenario. |
| AC-03 | One record only, including same-file multi-record | PASS: projector unit and same-file browser scenarios. |
| AC-04 | Strand-relative signed offsets and wrapping | PASS: pure resolver vectors. |
| AC-05 | Absolute/idempotent orientation | PASS: resolver, state, and browser assertions. |
| AC-06 | 3′ base differs from feature-end boundary | PASS: pure resolver boundary vectors. |
| AC-07 | Covered traversal midpoint and compound/origin cases | PASS: metadata/catalog and resolver vectors. |
| AC-08 | Operation-specific disabled reasons | PASS: resolver and popup availability assertions. |
| AC-09 | Stable duplicate/split identity | PASS: catalog/state/projection coverage. |
| AC-10 | One History entry for transform and Result | PASS: transaction unit and browser Undo/Redo journeys. |
| AC-11 | Session round trip | PASS: Session 44 save and fresh-page load browser journey. |
| AC-12 | Failure/cancel/stale preserves prior state | PASS: candidate unit and popup browser scenarios. |
| AC-13 | Pending form isolation | PASS: committed-base unit/browser assertions preserve unrelated drafts. |
| AC-14 | Zero additional LOSAT executor jobs | PASS: instrumented browser count remains `4 -> 4`; fresh Load is `0`. |
| AC-15 | One transform for all geometry | PASS: browser geometry and existing transform contracts. |
| AC-16 | Manual edits clear provenance | PASS: record-state unit coverage. |
| AC-17 | Cancel is a complete no-op | PASS: popup browser state/History assertions. |
| AC-18 | Search/rebind/keyboard/mobile continuity | PASS: real popup browser journey at desktop and 390 px. |
| AC-19 | Request 7, no new Worker/renderer path | PASS: request and architecture contracts. |
| AC-20 | Reviewable Product/architecture evidence and full gates | PASS for the candidate: Product authority, exact changed-scope evidence, and full local gates are recorded. Merge remains blocked until the maintainer posts the required architecture decision against the exact PR head. |

## Semantic owner and privileged path inventory

| Capability | Current owner/path | Planned bounded change |
| --- | --- | --- |
| Source location facts | `gbdraw/web_support/feature_metadata.py`, compacted by `feature_catalog.py` | Add precision/operator/order/strand capability facts before Biopython semantics are lost. |
| Feature anchor calculation | `record-display-options.js::selectedFeatureDisplayStart()` | Move the one calculation owner to private `app/record-display/feature-anchor.js`; keep the sidebar as an adapter. |
| Record display effective state | `app/record-display-options.js` | Extend the one draft/effective-state owner with absolute orientation and provenance. |
| Canonical request | `services/session-request.js` | Extend the one schema-7 projector and all-to-exact materialization owner. |
| Candidate execution/admission | `app/run-analysis.js` -> existing diagram-generation service/Worker -> current Result admission | Extract an owner-internal helper shared by Generate and the popup candidate. |
| SVG admission | `services/svg-result-ingestion.js` through the current Result admission path | Unchanged; no popup-specific insertion or sanitizer. |
| Artifact History | `services/history.js::runUndoableArtifactReplacement()` | Add optional intent checkpoint hooks to the same transaction owner. |
| Record coordinate transform | `gbdraw/layout/record_coordinates.py::RecordDisplayTransform` | Unchanged. |
| Popup presentation | `index.html`, composed from `app/app-setup.js` | Add binding-only markup plus one focused private controller. |

Privileged operators remain the current canonical request builder, diagram
Worker client, current Result admission owner, SVG ingestion owner, and History
transaction owner. The popup/controller must not become a privileged request,
Worker, SVG, Result, or History owner.

## Architecture ratchet ledger

The planned runtime is not an ordinary no-architecture-change patch because it
adds one persisted compatibility path. It therefore requires an
`ARCHITECTURE_EXCEPTION` decision on the exact final runtime head. The current
Session 01 documentation change itself adds no runtime owner, path, or
compatibility reader.

Planned changed-scope owner rows:

- Anchor capability facts: owners `{}` ->
  `{gbdraw/web_support/feature_metadata.py}`; `O 0 -> 1`, `T 0 -> 1`,
  `OE 0 -> 0`, delta 0.
- Feature anchor resolution: owners
  `{app/record-display-options.js}` ->
  `{app/record-display/feature-anchor.js}`; `O 1 -> 1`, `T 1 -> 1`,
  `OE 0 -> 0`, delta 0. The superseded formula must be removed.
- Record display effective state, canonical request, candidate admission, SVG
  admission, History, and coordinate transform each remain one owner;
  `OE 0 -> 0` for every row.

Planned path rows:

- Normal Generate remains
  `run-analysis.js -> session-request.js -> diagram-generation.js -> Worker -> current Result admission`.
- Popup Apply joins that same candidate execution/admission path after a
  target-only projection. It is an adapter into the existing path, not a second
  meaningful render path. `P 1 -> 1`, `PE 0 -> 0`.
- SVG admission, History transaction, and record-coordinate transformation
  remain `P 1 -> 1`, `PE 0 -> 0`.

Compatibility row:

- Namespace: saved Session editor feature catalog.
- Stable ID before: none for catalog 3 -> 4.
- Stable ID after:
  `saved-session:v42-catalog3-to-v44-catalog4`.
- `CB 0 -> 1`, delta `+1` in this changed namespace.
- Positive fixtures:
  `gbdraw/web/gallery/sessions/lambda_basic_linear.gbdraw-session.json` and
  `gbdraw/web/gallery/sessions/hepatoplasmataceae_orthogroup.gbdraw-session.json.gz`
  are Session 42 documents with request schema 7, saved Result, and catalog
  schema 3 on current `main`.
- Reader location: the existing Session migration/admission owner only. Current
  runtime catalog admission remains schema 4 after normalization.
- Implemented compatibility behavior has one stable ID and one conservative
  rule, realized at the two existing language boundaries that can write or
  admit a Session. Browser admission uses
  `gbdraw/web/js/services/feature-catalog.js::migrateLegacyFeatureCatalog()`
  from `config.js::preflightSessionImport()`. Python replay/current-write and
  Gallery preparation use
  `gbdraw/web_support/feature_catalog.py::promote_legacy_feature_catalog()`
  from their existing owners. These are parity-constrained realizations of the
  same path, not independent migration chains. There is no Session 43 reader.
- Current writer and normalized runtime: Session 44, feature catalog 4, and
  canonical request 7. Catalog 4 stores `anchorProfile`; legacy single-part
  integer intervals are promoted conservatively, while compound legacy
  locations receive the explicit `precision: unavailable` regeneration state.
- Positive-fixture result: both named released Session 42/catalog 3 Gallery
  documents remain unchanged on disk and are admitted through the bounded
  reader without replacing their saved Results. A subsequent Save writes
  Session 44/catalog 4.
- Removal condition: remove after the repository's declared support window no
  longer includes any released Session capable of carrying catalog schema 3,
  with release/first-parent evidence and fixtures updated in the same change.

Exact changed-scope totals for the planned exception are:

```text
OE 0 -> 0; delta(OE) = 0
PE 0 -> 0; delta(PE) = 0
CB 0 -> 1; delta(CB) = +1
```

Superseded semantic owner: the anchor formula in
`record-display-options.js::selectedFeatureDisplayStart()`.
Superseded canonical paths: none. Superseded compatibility paths: none.

### Corrected compatibility premise

`main` @ `4556e04e` writes Session 42 and feature catalog 3. Session 43 exists
only on `origin/dev` at this preflight and no tag contains its writer commit.
Repository policy forbids creating a reader for a branch-only intermediate
format. The saved plan was minimally corrected so Session 03 will:

1. write Session 44/catalog 4;
2. migrate released Session 42/catalog 3 conservatively;
3. rewrite any branch-owned v43 artifacts to the current writer;
4. add no v43-only reader unless v43 is independently present in first-parent
   `main` or a release tag at the later implementation start.

## File-level implementation map

- `gbdraw/web_support/feature_metadata.py`: currently maps processed feature
  locations back to source intervals and detects fuzziness for sequence
  warnings. It must expose only capability facts, not resolved anchors.
- `gbdraw/web_support/feature_catalog.py`: schema constant, biological feature
  compaction, catalog validation, and emitted envelope. It is the Python-side
  schema-4 writer/validator change point.
- `gbdraw/web/js/services/feature-catalog.js`: current schema-3 admission,
  compact-part expansion, runtime projection, and stable identity indexes. It
  must admit only the normalized current catalog after Session migration.
- `gbdraw/web/js/app/record-display-options.js`: current draft validation,
  source freshness, effective topology, sidebar shortcut, and manual edit
  owner. It will own effective orientation/provenance state and delegate anchor
  math to the pure resolver.
- `gbdraw/web/js/services/session-request.js`: schema-7 request projection,
  display-difference exact-one materialization, current config/session
  projection, and future target-only committed-base projector.
- `gbdraw/web/js/app/run-analysis.js`: current canonical construction,
  comparison planning/cache, Worker execution, typed response, catalog
  admission, and Result commit orchestration. Generate and popup candidates
  must share one extracted helper here.
- `gbdraw/web/js/services/history.js`: existing artifact replacement before/
  after handles, rollback, Undo, and Redo. Optional record-intent hooks belong
  here.
- `gbdraw/web/index.html`: current rich/simple feature popup and bindings. It
  receives presentation-only Record actions markup.
- `gbdraw/web/js/app/app-setup.js`: dependency wiring and exported popup action
  controller.
- `gbdraw/web/js/services/config.js`: Session 44 writer and bounded released
  Session 42/catalog 3 migration owner.

## Fixture inventory

| Need | Reuse/minimal addition |
| --- | --- |
| Plus/minus, odd/even, multipart gaps, origin-spanning | Reuse the synthetic Biopython and JS constructions in `test_web_feature_metadata.py` and `record-display-options.test.mjs`; add focused in-memory cases, not a data file. |
| Unstranded, mixed, fuzzy, `order`/unknown | No current fixture retains the complete capability distinctions. Add the smallest in-memory `SeqFeature` cases to metadata/catalog tests. |
| Duplicate record IDs | Reuse the two `same` records in `record-display-options.test.mjs` and the existing duplicate-instance SVG/browser contracts. |
| Same-file multi-record and independent chromosomes | Reuse the `Vnig_TUMSAT-TG-2018` Gallery/session source and the synthetic complete-record browser upload. |
| Same-row multi-record | Reuse `linear-multi-record.playwright.spec.js` record-row fixtures. |
| Circular and Linear mode | Reuse `joint-display-placement.playwright.spec.js` setup and the current circular record-presentation fixtures. |
| LOSAT comparison/cache | Reuse the instrumented complete-record executor in `linear-multi-record.playwright.spec.js`; executor invocation, not cache size, is the metric. |
| Depth/statistics/ticks/labels | Reuse `visual-state-regressions.playwright.spec.js`, depth-track fixtures, and existing shared-transform SVG contracts; add only the rotation journey assertions. |
| Released Session/catalog compatibility | Reuse the two Session 42/catalog 3 Gallery sessions named in the compatibility ledger. |

No new persistent biological fixture is required unless a later focused test
cannot express an ambiguous Biopython location in memory.

## Session 01 verification

| Command | Result |
| --- | --- |
| `node --test tests/web/record-display-options.test.mjs` | PASS |
| `node --test tests/web/session-request.test.mjs` | PASS |
| `node --test tests/web/history.test.mjs tests/web/run-analysis-simple-path.test.mjs` | PASS |
| `node tests/web/architecture-contracts.test.mjs` | PASS, 137 tests |
| `npx playwright test ... --grep "LOSATP source jobs are reused"` | Environment failure before test discovery: `@playwright/test` is not installed. |
| `python /tmp/verify_issue563_session01_losat.py` | PASS equivalent Python Playwright journey after required sandbox escalation; counts `4, 4, 4, 0`, zero additional jobs after start/reverse and zero jobs after fresh Load. |
| `git diff --check` | PASS |

## Satisfied stop condition and resumption

Session 02 runtime work was allowed because all of the following were satisfied:

1. `satoshikawato` supplied a complete response from
   `PRODUCT_DECISION_PACKET.md`;
2. `PD-OI-032` serializes that choice without inferred fields in an
   authority-only change;
3. PR `#572` is merged and merge commit `fbf195f2` places that authority in the
   runtime branch ancestry through the explicitly authorized branch update;
4. no conflicting authority has appeared.

The compatibility exception still requires an exact-final-head architecture
decision before merge. The user subsequently authorized push, PR creation, and
merge for this branch; that authorization does not replace the mandatory
manual architecture decision.

## Session 07 verification

The complete command/result table and environment are recorded in
`FINAL_ACCEPTANCE.md`. The decisive full gates were:

- Python non-slow suite: `6132 passed`, `17 skipped`, `11 deselected` in
  `761.38s` after the required Chromium sandbox escalation;
- Web Node suite: `641 passed` in `74.636s`;
- architecture contracts: `137 passed` in `65.049s`;
- CI contracts: `59 passed` in `13.448s`;
- targeted Playwright: `1 passed` for linear multi-record and `3 passed` for
  the real popup journey;
- Ruff: PASS;
- LOSAT executor additions: `0`.

The exact-head change-budget report and manual architecture decision are
post-commit/PR gates because the reviewed SHA must identify the final commit.
