# Session 01 instruction prompt — show the underlying alignment failure

Paste this entire file into a new agent session. It is self-contained.

## Mission

When a Linear Similarity Group alignment fails during generation in the Web
app, the review palette and the global error banner show only
`Alignment generation failed. Review the draft and retry Apply.`, and the
palette then appends a second `Review the choices and retry Apply.` The real
cause is lost:

1. `runAnalysis` (`gbdraw/web/js/app/run-analysis.js`) stores the normalized
   error in `errorLog` and returns `{ status: 'error' }` without it.
2. `applyPlan()` in `gbdraw/web/js/app/similarity-alignment.js` then calls
   `publishError()` with a generic `Error`.
3. The `onError` wiring in `app-setup.js` overwrites `errorLog` with that
   generic error.

This session makes both surfaces show the same underlying failure once. It
changes no alignment, orientation, or plan behavior. Existing authority
`PD-OI-031` revision 3 already requires actionable errors, so the
classification is `IMPLEMENT_EXISTING_AUTHORITY`. This session does not wait
for Session 00.

## Branch

Use exactly `fix/similarity-alignment-orientation-owner-20260926`:

1. Fetch `origin`.
2. Check out the existing branch.
3. Verify that its upstream is `origin/fix/similarity-alignment-orientation-owner-20260926`.
4. Merge newer `origin/dev` into it if needed. Never rebase or recreate it.

Push only to that same-named remote branch.

## Read and map

Read:

- `AGENTS.md`, `CLAUDE.md`, `gbdraw/web/CLAUDE.md`;
- `docs/internal/similarity-alignment-orientation-owner-20260926/01_MASTER_PLAN.md`
  (sections 1, 2, 3 item 9, and 4.4).

Then map every `return { status: 'error' }` in `runAnalysisInternal()` and
`runAnalysis()`, and every caller that consumes a `runAnalysis` outcome.
Identify the single normalizer, `normalizeUserFacingError()` in
`gbdraw/web/js/services/error-normalization.js`, and every place that writes
`errorLog`.

## Implementation

1. Wherever `runAnalysis` publishes an error to `errorLog`, return that same
   normalized value in the outcome: `{ status: 'error', error }`. Do not create
   a second normalization.
2. In the alignment controller, when an outcome carries `error`:
   - show that error in the review;
   - do not overwrite the banner `runAnalysis` already set;
   - keep the generic text only as the fallback when no error is supplied.
   Keep resolver and validation errors, which `runAnalysis` never saw, on the
   existing `publishError()` path.
3. In `gbdraw/web/index.html`, remove the appended
   `Review the choices and retry Apply.` sentence. The message needs to appear
   only once.
4. Preserve everything else:
   - draft retention and retry;
   - stale, cancel, and superseded guards;
   - last-Result preservation;
   - History behavior.

Aim for a net reduction or a minimal increase in production lines.

## Verification

- Add or extend unit tests in `tests/web/similarity-alignment-actions.test.mjs`:
  - an automatic-path render failure and an Apply render failure both expose
    the supplied error message;
  - both keep the draft for retry;
  - neither overwrites the banner with the generic text.
- Add a focused `run-analysis` test asserting that the error outcome carries
  the same normalized error stored in `errorLog`.
- Confirm the root cause in a real browser, and do not commit a test that
  asserts the current failure:
  1. Run `gbdraw gui` and load the Gallery example `BGC0000708-BGC0000713`.
  2. Click `CAG38712.1` in the first record and choose
     **Review alignment options…**.
  3. Enable Match reference direction for *Streptomyces fradiae* ATCC 10745 and
     Apply.
  4. Record that both the palette and the banner now read
     `Comparison source feature index conflicts with its view feature ID.`
  - Check both Playwright installations as described in `CLAUDE.md`. Use Python
    Playwright if Node's `@playwright/test` is missing, and rerun a Chromium
    sandbox failure with escalation.
- Run the focused Node tests, the existing alignment Chromium spec,
  `node tests/web/architecture-contracts.test.mjs`,
  `node tools/check-web-change-budget.mjs --base origin/dev`, and
  `git diff --check`.
- Review production and test diffs separately.

Commit the verified change on the fix branch. Push it only to the same-named
remote branch.

## Handoff

Report the branch and upstream, commit, changed files, test results, and the
observed browser message.

Also state the authority gate. Sessions 02–04 start only after the Session 00
authority for `PD-OI-027/028/029/031/034` has merged into `origin/dev` and
`origin/dev` has been merged into this branch.

At the end, print the entire contents of
`docs/internal/similarity-alignment-orientation-owner-20260926/SESSION_02_PLAN_CONTRACT.md`
in a copyable code block as the prompt for the next session.
