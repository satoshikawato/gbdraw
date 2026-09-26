# Session 02 instruction prompt — remove orientation from the alignment plan

Paste this entire file into a new agent session. It is self-contained.

## Mission

In gbdraw, a Linear Similarity Group alignment plan
(`gbdraw/layout/similarity_alignment.py`) currently stores, per target:

- a requested `orientation_policy` (`preserve` or `match_reference`);
- an `effective_reverse_complement` value.

Python applies that value late, in `materialize_similarity_alignment_display()`
(`gbdraw/api/record_planning.py`). By then the browser has already projected
protein-comparison rows into displayed coordinates from each record's Reverse
setting. Any plan-driven reversal therefore fails with `ValueError: Comparison
source feature index conflicts with its view feature ID.`

The selected design makes each record's Reverse setting the only orientation
owner. This session makes the plan and the render path orientation-free across
Python and the Web contract readers. The per-target orientation checkbox is
also removed, because its data no longer exists. Session 03 then adds the
single draft-level **Match reference direction** option, which works through
record state.

## Branch and gate

1. Use exactly `fix/similarity-alignment-orientation-owner-20260926`, with
   upstream `origin/fix/similarity-alignment-orientation-owner-20260926`.
2. Fetch `origin`. Verify that the approved authority for
   `PD-OI-027/028/029/031/034` (contract revision 16 in
   `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`) is merged into
   `origin/dev`.
3. Merge `origin/dev` into the fix branch. Never rebase or recreate the branch.
4. If the authority is not merged, stop and report the gate. The Decision Pack
   and any unmerged authority branch are not runtime authority.
5. Confirm the Session 01 error-reporting commit is present.

## Read

Read:

- `AGENTS.md`, `CLAUDE.md`, `gbdraw/web/CLAUDE.md`;
- the merged Product records;
- `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md`;
- `docs/internal/similarity-alignment-orientation-owner-20260926/01_MASTER_PLAN.md`
  (sections 1, 3, 4, 5, and 6).

Run the reproduction in master-plan section 1 before changing code.

## Implementation

Apply the contracts in master-plan section 4.3. In the same change, delete the
superseded code instead of leaving it unused.

### Python

1. **Domain.** In `gbdraw/layout/similarity_alignment.py`:
   - remove `AlignmentOrientationPolicy`, `AlignmentOrientationEffect`,
     `_orientation_result()`, `orientation_policy` from choices and decisions,
     `effective_reverse_complement` from decisions and candidates, and the
     `match_reference_*` review fields;
   - add `AlignmentStrandRelation` (`same`, `opposite`, `unknown`) and store one
     `strand_relation` per `AlignmentReviewCandidate`, computed once from the
     reference and candidate displayed strands;
   - update validation and `__all__`;
   - leave anchor resolution, priority, and recommendations unchanged.
2. **Render.** Replace `materialize_similarity_alignment_display()` with
   `project_similarity_alignment_centers(collection, plan)`. It returns only the
   projected anchor centers and never modifies records. Update
   `plan_linear_request()` in `gbdraw/api/request_render.py` and every other
   caller so that the resolved collection is used unchanged.
3. **Adapters and codecs.**
   - `gbdraw/web_support/similarity_alignment.py`: choices without
     `orientationPolicy`; candidates with `strandRelation` replacing
     `orientation`; decisions and plan without orientation fields.
   - `gbdraw/session_request_codec.py`: encode and decode the plan without
     orientation fields.
   - `gbdraw/api/session_compat.py`: the legacy plan builder emits none.
   - `gbdraw/api/record_planning.py`: the CLI candidate builder drops
     `effective_reverse_complement`.
   - `gbdraw/api/__init__.py`: remove the `AlignmentOrientationPolicy` export.
     Keep plan and helper `schema: 2`, and add no reader for the removed fields
     (master-plan section 5).

### Web contract readers

1. In `gbdraw/web/js/app/similarity-alignment.js`:
   - validate the new candidate, decision, and plan shapes;
   - remove `orientationPolicy` and `orientationEffect` from review rows and
     choices (`reviewRows()`, `editRow()`, `applyDraft()`, `planChoices()`);
   - remove the `orientation` edit branch, the `setOrientation` export, the
     plan-driven parts of `orientationsFromRequest()` and
     `requestWithOrientations()`, and the use of plan orientation in
     `successfulSummary()` and `inspectActivePlan()`;
   - keep the summary's `reversed` count, which reports `0` until Session 03.
2. Make the same shape change in `gbdraw/web/js/services/session-request.js`
   (plan validator) and `gbdraw/web/js/services/legacy-similarity-alignment.js`.
3. In `gbdraw/web/index.html`, remove the per-target Match reference direction
   checkbox and its `Effective orientation` line. Session 03 adds the
   replacement UI.
4. In `gbdraw/web/js/app/app-setup.js`, change `linearRecordOrientationValue()`
   to stop reading the plan. Session 03 removes it.

### Generated artifacts and documentation

1. Run the unfiltered Gallery owner command
   `python tools/refresh_gallery_sessions.py`, then
   `python tools/gallery_artifact_manifest.py`.
   - Confirm that only expected fields and bytes changed.
   - The Gallery session is generator-owned and may be overwritten.
2. Update documentation this change makes inaccurate:
   - `docs/REFERENCE/python-api.md`: rewrite the executable example
     `S07-PY-01` so it reads the target's `strand_relation` from the resolution
     review rows, reverses the target through
     `RecordPresentation(reverse_complement=True)`, and renders with the plan.
     Change its printed line and the expectation in
     `tests/test_documentation_reference_contracts.py` together.
   - `docs/REFERENCE/typed-requests.md`.
   - `docs/REFERENCE/session-and-request-compatibility.md`.
   - `docs/SESSION_COMPATIBILITY.md`.
   - Do not edit historical plan directories.

## Tests

Update or replace orientation-policy tests in:

- `tests/test_similarity_alignment.py`
- `tests/test_similarity_alignment_rendering.py`
- `tests/test_similarity_alignment_web_adapter.py`
- `tests/test_session_request_codec.py`
- `tests/test_session_compat.py`
- `tests/test_api_session.py`
- `tests/web/similarity-alignment-actions.test.mjs`
- `tests/web/session-request.test.mjs`
- `tests/web/run-analysis-simple-path.test.mjs`
- `tests/web/gallery-session-publication.test.mjs`
- `tests/web/similarity-alignment-ui.playwright.spec.js`

They must cover:

- `strand_relation` for same, opposite, and unknown strands on either side;
- plan and choice validation rejecting the removed fields;
- a rendering test showing that a plan never changes record orientation and
  that anchor centers align exactly for a target reversed through
  `RecordPresentation.reverse_complement`;
- Session and legacy round trips without orientation fields.

Remove browser assertions that expect a per-target checkbox; Session 03 adds
the new journey.

Run:

- the focused Python and Node tests;
- the alignment Chromium spec (check both Playwright installations; use Python
  Playwright if Node's runner is missing);
- `ruff check gbdraw/`;
- `node tests/web/architecture-contracts.test.mjs`;
- `node tools/check-web-change-budget.mjs --base origin/dev`;
- `git diff --check`.

Review production, tests, documentation, and generated artifacts as separate
diffs. Commit on the fix branch and push only to its same-named remote branch.

## Handoff

Report:

- the merged authority SHA and branch state;
- commits and changed files;
- test results and the Gallery refresh result;
- the line-count delta for production code;
- concise owner/path evidence (master-plan section 9).

At the end, print the entire contents of
`docs/internal/similarity-alignment-orientation-owner-20260926/SESSION_03_ORIENTATION_WORKFLOW.md`
in a copyable code block as the prompt for the next session.
