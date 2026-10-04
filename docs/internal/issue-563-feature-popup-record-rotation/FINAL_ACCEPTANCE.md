# Issue #563 final acceptance

## Candidate identity

- Branch: `issue-563-feature-popup-record-rotation-20260922`
- Session 07 parent HEAD: `10075e3d0f92ec7a69ab9d5680c2af35a12e5af2`
- Candidate HEAD: the Session 07 commit containing this report; its exact SHA
  is recorded in the PR and mandatory architecture decision because a commit
  cannot contain its own SHA.
- Integration base merged into the branch:
  `origin/dev` @ `35b2716a33c726cac82842b839383f6ab294b3e6`
- PR target base observed after local gates:
  `origin/dev` @ `1cb2a1b00ac97fe8f1213d12a63257ead0d56ded`
- Upstream before publication: none (expected)
- Product authority: `PD-OI-032`, merged by PR #572 and present through
  `fbf195f2`

## User-visible outcome

The feature popup now has a separate **Record actions** workflow. A user can
place the open popup's feature 5′ end, covered-traversal midpoint, or 3′ base at
the display start with a signed biological-direction offset, optionally orient
the feature forward, or place the outgoing feature-end boundary at the start.
The preview reports the source coordinate and absolute orientation before
submission.

Apply regenerates only the explicitly captured `(recordKey,
biologicalFeatureId)` record from the committed base. It leaves unrelated form
drafts and other records unchanged, and commits the record transform and Result
as one Undo/Redo entry. Cancel, failure, stale source, and superseded completion
keep the previous Result and transform. Save Session and fresh Load preserve the
successful transform and provenance. Circular records work in Circular and
Linear diagram modes; unsafe topology, crop, location, or freshness states are
disabled with an operation-specific reason.

## Architecture and ownership

- Source capability facts: `gbdraw/web_support/feature_metadata.py`, compacted
  and validated by `gbdraw/web_support/feature_catalog.py`.
- Anchor calculation: one pure owner,
  `gbdraw/web/js/app/record-display/feature-anchor.js`; the superseded formula
  was removed from `record-display-options.js`.
- Effective record state and provenance:
  `gbdraw/web/js/app/record-display-options.js`.
- Canonical request: the existing schema-7 owner,
  `gbdraw/web/js/services/session-request.js`.
- Candidate execution and Result admission: the existing
  `gbdraw/web/js/app/run-analysis.js` path through the existing Worker and SVG
  admission owners. Popup Apply adds no request, Worker, renderer, sanitizer,
  SVG admission, or Result path.
- Atomic artifact/intent transaction: the existing
  `history.js::runUndoableArtifactReplacement()` owner.
- Coordinate projection: the unchanged Python `RecordDisplayTransform` owner.
- Popup markup/controller: binding-only markup in `index.html` and the focused
  controller wired by `app-setup.js`.

The catalog-3 compatibility behavior has one stable ID and conservative rule:
`saved-session:v42-catalog3-to-v44-catalog4`. Because both browser admission and
Python replay/current-write can independently encounter released Session 42,
the rule has parity-constrained JavaScript and Python realizations at those
existing boundaries. This avoids a new migration chain while preserving the
same exact-single-only promotion and compound-unavailable outcome.

Exact changed-scope evidence:

```text
OE rows: all changed capabilities 0 -> 0; total delta(OE) = 0
PE rows: Generate/popup candidate, SVG admission, History, and transform 0 -> 0; total delta(PE) = 0
CB row saved-session-feature-catalog: {} -> {saved-session:v42-catalog3-to-v44-catalog4}; 0 -> 1; delta(CB) = +1
```

Superseded semantic owner: the feature-anchor formula formerly in
`record-display-options.js::selectedFeatureDisplayStart()`. Superseded
canonical paths: none. Superseded compatibility paths: none. The positive CB
delta requires the manual exact-head architecture exception decision before
merge.

## Issue scenarios

| # | Scenario | Result and evidence |
| --- | --- | --- |
| 1 | Popup-only operation | PASS — real rich/simple popup pointer and keyboard journeys use the captured target without sidebar/global selection. |
| 2 | Circular and Linear modes | PASS — browser acceptance resolves the same source coordinate in both modes. |
| 3 | Independent chromosomes/records | PASS — same-file multi-record acceptance changes one record and preserves the other. |
| 4 | Strand-aware offsets | PASS — resolver vectors cover plus/minus signed offsets and circular wrapping. |
| 5 | Stable orientation | PASS — orientation is absolute/idempotent; manual edits clear old provenance. |
| 6 | 3′ base versus feature end | PASS — resolver and popup tests distinguish the covered base from the outgoing boundary. |
| 7 | Circular compound locations | PASS — catalog capability facts and resolver vectors cover multipart gaps, origin wrapping, and odd/even midpoints. |
| 8 | Safe eligibility and recovery | PASS — operation-specific disabled reasons plus cancel/failure/stale/superseded preservation are covered. |
| 9 | Stable identity and isolation | PASS — duplicate IDs, split fragments, and exact record projection retain source-bound identity. |
| 10 | State, rendering, and persistence | PASS — shared request/admission/History owners, fresh Load, geometry assertions, and zero added LOSAT jobs are covered. |

## Acceptance criteria

| AC | Result | Evidence |
| --- | --- | --- |
| AC-01 | PASS | Open-popup target only; no global-selection fallback. |
| AC-02 | PASS | Same source coordinate in Circular and Linear. |
| AC-03 | PASS | Target-only same-file multi-record projection. |
| AC-04 | PASS | Plus/minus signed offsets and wrap vectors. |
| AC-05 | PASS | Absolute, idempotent orientation. |
| AC-06 | PASS | 3′ covered base and outgoing boundary differ. |
| AC-07 | PASS | Compound/origin and odd/even covered traversal. |
| AC-08 | PASS | Operation-specific unavailable reasons. |
| AC-09 | PASS | Duplicate/split stable source identity. |
| AC-10 | PASS | Transform and Result share one History entry. |
| AC-11 | PASS | Session 44 save and fresh-page Load round trip. |
| AC-12 | PASS | Failure/cancel/stale/superseded preserve prior state. |
| AC-13 | PASS | Unrelated pending edits remain pending. |
| AC-14 | PASS | Instrumented LOSATP count stays `4 -> 4`; fresh Load uses `0`; additional jobs = `0`. |
| AC-15 | PASS | Feature/label/tick/depth/statistics/comparison geometry shares the existing transform. |
| AC-16 | PASS | Manual start/orientation edits clear feature provenance. |
| AC-17 | PASS | Cancel changes no draft, Result, or History state. |
| AC-18 | PASS | Search/rebind/keyboard and 390 px acceptance. |
| AC-19 | PASS | Request 7 reused; no new Worker/renderer/admission path. |
| AC-20 | PASS for candidate | Product authority, full gates, and exact architecture evidence are reviewable; merge still requires the manual exact-head exception decision. |

## Verification

Environment: Python 3.13.3 in `/tmp/gbdraw-issue563-testenv` with the candidate
installed editable, Node v26.8.2, and Playwright 1.61.1. Browser gates were run
outside the agent sandbox as required by repository guidance.

| Command | Result |
| --- | --- |
| `python -m pytest tests/test_web_feature_metadata.py tests/test_web_feature_catalog.py tests/test_session_io.py -v` | PASS — 280 passed in 3.09s. |
| `node --test tests/web/feature-anchor.test.mjs tests/web/record-display-options.test.mjs` | PASS — 2 files in 86.58ms. |
| `node --test tests/web/session-request.test.mjs tests/web/history.test.mjs tests/web/run-analysis-simple-path.test.mjs` | PASS — 3 files in 299.03ms. |
| `npx playwright test tests/web/linear-multi-record.playwright.spec.js --grep "record rotation\|LOSATP source jobs" --workers=1 --retries=0` | PASS — 1 passed in 24.2s. |
| `npx playwright test tests/web/interactive-svg-v3.playwright.spec.js --grep "record rotation" --workers=1 --retries=0` | PASS — 3 passed in 35.6s. |
| `node --test tests/web/*.test.mjs` | PASS — 641 passed in 74.636s. |
| `python -m pytest tests/ -v -m "not slow"` | PASS — 6132 passed, 17 skipped, 11 deselected, 17 warnings in 761.38s. |
| `ruff check gbdraw/` | PASS. |
| `node tests/web/architecture-contracts.test.mjs` | PASS — 137 passed in 65.049s. |
| `node --test tests/ci/*.test.mjs` | PASS — 59 passed in 13.448s. |
| `node tools/check-web-change-budget.mjs --base origin/dev --head HEAD` | Post-commit exact-head gate; result appended to the PR handoff. |
| `git diff --check` | PASS before Session 07 commit; repeated after final documentation. |

One unprivileged full-suite attempt was stopped at 45% after expected Chromium
sandbox failures. It is not acceptance evidence; the complete escalated rerun
above is the decisive result.

## Schema and compatibility outcome

- Canonical request schema: unchanged at 7.
- Current Session writer: 44.
- Current normalized feature catalog: 4 with source anchor capability facts.
- Released reader: Session 42/catalog 3 to Session 44/catalog 4 only.
- Branch-only Session 43 reader: none.
- Legacy exact single intervals are promoted; compound/ambiguous legacy
  locations remain explicitly unavailable for feature-based rotation.
- Standalone interactive SVG payloads and Gallery assets now consume catalog 4.

## Artifact disposition and residual risk

Tracked documentation and Gallery/recipe artifacts were regenerated because
their embedded current catalog payload changed to schema 4. They were validated
by the recipe, Gallery, packaging, standalone-browser, and full regression
tests. No `tests/reference_outputs/`, `dist/`, or `gbdraw.egg-info/` change is
included. The prepared browser wheel is generated and gitignored and is not
committed. No new dependency, build step, CDN, or CSP change was introduced.

Accepted residual product risk is the bounded popup/catalog compatibility UI
and maintenance cost recorded in `PD-OI-032`. Open merge boundary: the
architecture maintainer must manually approve the `CB +1` exception for the
exact final PR head after CI. Other unresolved product, scientific, runtime, or
test issue: none.
