# Issue #563 baseline and preflight

- Status: `BLOCKED_PRODUCT_AUTHORITY`
- Session: 01
- Recorded: 2026-09-22
- Branch: `issue-563-feature-popup-record-rotation-20260922`
- Session start HEAD: `3672a518e09339bffc2cc70cdee42add9103ccb2`
- Base: `origin/dev` @ `a9eaeadd105e0c26e46086626feaa695bdd33c94`
- Upstream: none
- Issue: `#563 Add feature-based rotation of circular records from the feature popup`

## Outcome

Session 01 is complete. Sessions 02–07 are blocked. The latest Issue #563 body
specifies one detailed proposed workflow, but no authority already present on
the branch base selects that new public affordance and complete outcome. The
current implementation request also does not contain the complete Product
Decision Owner receipt required by the Product Impact Ratchet. Runtime files
must not change until the selected outcome has been recorded in an
authority-only change and that authority is present on the runtime branch's
base.

The ready-to-review choices and required response are in
`PRODUCT_DECISION_PACKET.md` beside this document.

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

## Confirmed current behavior

- `RecordPresentation.reverse_complement` and
  `RecordDisplayOptions.start_coordinate` are the typed absolute display
  values. Request schema 7 already carries both without a new Issue-specific
  request field.
- `RecordDisplayTransform` owns source/display base and boundary conversion for
  circular and reverse-complemented records.
- `record-display-options.js` owns the source-bound editable display draft and
  existing selected-feature 5′/midpoint shortcuts. Its current draft has only
  topology and start, and its shortcut depends on the global selection.
- `session-request.js` is the single canonical request owner. It materializes a
  `cardinality: all` record into exact-one records when per-record display
  values differ.
- `run-analysis.js` builds the canonical request and then uses the existing
  diagram Worker, typed response, catalog validation, and Result admission
  path.
- `history.js::runUndoableArtifactReplacement()` is the artifact replacement
  transaction owner.
- Session current writer is 43 on `origin/dev`; feature catalog current schema
  is 3; canonical request schema remains 7.
- The characterization journey submitted four LOSATP jobs initially. Changing
  display start and then reverse complement left the cumulative count at four;
  a fresh Session Load submitted zero jobs. A target-only candidate is therefore
  needed for pending-edit isolation and atomic History, not to invent a LOSAT
  bypass.

## Product Impact developer preflight

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

| Issue scenario | Planned contract | Current evidence |
| --- | --- | --- |
| 1 Popup-only | AC-01, AC-18 | Popup exists; no record action exists. |
| 2 Both diagram modes | AC-02, AC-15 | Shared typed display transform exists. |
| 3 Independent chromosomes | AC-03, AC-09 | Stable record keys and exact selectors exist. |
| 4 Strand-aware offsets | AC-04 | Current shortcut resolves anchors only; signed offset is new. |
| 5 Stable orientation | AC-05, AC-16 | Absolute typed reverse-complement value exists; popup operation is new. |
| 6 Endpoint distinction | AC-06 | No current 3′ or feature-end popup action. |
| 7 Circular/compound | AC-07 | Parts exist, but precision/operator/order facts are incomplete. |
| 8 Safe eligibility | AC-08, AC-12, AC-17 | Existing source replacement/crop checks are reusable; operation reasons are new. |
| 9 Identity/isolation | AC-03, AC-09 | Source-bound record and biological feature identities already exist. |
| 10 State/rendering | AC-10–AC-16, AC-19–AC-20 | Existing request, transform, Result, History, Session, and LOSAT cache boundaries are reusable. |

| AC | Required outcome/evidence | Session 01 result |
| --- | --- | --- |
| AC-01 | Explicit open-popup target, no global selection fallback | Traced; runtime blocked. |
| AC-02 | Same source coordinate in Circular and Linear | Existing shared transform confirmed; runtime test pending. |
| AC-03 | One record only, including same-file multi-record | Stable exact selectors confirmed; runtime test pending. |
| AC-04 | Strand-relative signed offsets and wrapping | Formula fixed by Issue; resolver pending. |
| AC-05 | Absolute/idempotent orientation | Existing absolute request field confirmed; resolver pending. |
| AC-06 | 3′ base differs from feature-end boundary | Product meaning specified; resolver pending. |
| AC-07 | Covered traversal midpoint and compound/origin cases | Existing parts insufficient for safe order/precision classification; catalog 4 pending. |
| AC-08 | Operation-specific disabled reasons | Existing crop/source checks confirmed; new reasons pending. |
| AC-09 | Stable duplicate/split identity | Existing `(recordKey, biologicalFeatureId)` confirmed. |
| AC-10 | One History entry for transform and Result | Existing artifact transaction owner confirmed; hook pending. |
| AC-11 | Session round trip | Existing draft/session projection confirmed; v44 work pending. |
| AC-12 | Failure/cancel/stale preserves prior state | OIPC-C07 and current Result admission authority confirmed. |
| AC-13 | Pending form isolation | Requires committed-base projector; pending. |
| AC-14 | Zero additional LOSAT executor jobs | Characterized: `4 -> 4 -> 4`, fresh Load `0`. |
| AC-15 | One transform for all geometry | `RecordDisplayTransform` confirmed as owner. |
| AC-16 | Manual edits clear provenance | New provenance behavior pending. |
| AC-17 | Cancel is a complete no-op | New controller behavior pending. |
| AC-18 | Search/rebind/keyboard/mobile continuity | Existing popup/search identity confirmed; new journey pending. |
| AC-19 | Request 7, no new Worker/renderer path | Feasible and required; architecture gate currently passes. |
| AC-20 | Reviewable Product/architecture evidence and full gates | Product authority unresolved; Session 01 gates pass as recorded below. |

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

## Stop condition and resumption

Do not begin Session 02 runtime work. Resume only after all of the following:

1. an identified Product Decision Owner supplies a complete response from
   `PRODUCT_DECISION_PACKET.md`;
2. that choice is serialized without inferred fields in an authority-only
   change;
3. the authority-only change is merged and is present on the runtime branch's
   base through an explicitly authorized branch update; and
4. no conflicting authority has appeared.

The later compatibility exception still requires an exact-final-head
architecture decision before merge. Push, PR creation, merge, rebase, and
branch update remain outside the current authorization.
