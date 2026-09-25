# Session 03 instruction prompt — one Match reference direction option through record orientation

Paste this entire file into a new agent session. It is self-contained.

## Mission

In gbdraw's Web app, a Linear Similarity Group alignment plan now stores
anchors only. Each record's **Reverse complement** setting
(`linearSeqs[].region_reverse`) is the only source of orientation. The
previous session removed the per-target Match reference direction checkbox.
Each review candidate now carries a Python-computed `strandRelation`
(`same`, `opposite`, or `unknown`) relative to the reference anchor.

This session delivers the selected workflow:

- one draft-level **Match reference direction** option;
- Apply reverses eligible targets by changing their record orientation, in one
  generation and one History entry;
- a manual Reverse after alignment keeps the plan;
- Reset Align leaves orientation unchanged.

## Branch and gate

1. Use exactly `fix/similarity-alignment-orientation-owner-20260926`, with
   upstream `origin/fix/similarity-alignment-orientation-owner-20260926`.
2. Fetch `origin`.
3. Confirm the merged authority for `PD-OI-027/028/029/031/034` (contract
   revision 16) is on the branch, and that the Session 01 and Session 02
   commits are present.
4. Merge newer `origin/dev` if needed. Never rebase or recreate the branch.

## Read

Read:

- `AGENTS.md`, `CLAUDE.md`, `gbdraw/web/CLAUDE.md` (especially the reactive
  availability, live-edit, and History rules);
- the merged Product records;
- `docs/internal/similarity-alignment-orientation-owner-20260926/01_MASTER_PLAN.md`
  (sections 3, 4, and 8).

## Implementation

### Draft option

In `gbdraw/web/js/app/similarity-alignment.js`:

- Add `matchReferenceDirection: false` to every newly created draft and a
  controller action `setMatchReferenceDirection(enabled)`. It is a local edit
  and starts no Worker job.
- Add one pure function, `matchedOrientations(response, request, enabled)`. It
  returns the complete `{recordKey, reverseComplement}` list:
  - current orientation for every record;
  - flipped only for aligned targets whose selected anchor has
    `strandRelation === 'opposite'`, and only when `enabled` is true.
  This is the only place the "reverse only when known opposite" rule is
  applied.
- Derive the draft's view facts from the same function or the same facts:
  - the reversal list, the unknown-strand list, and each row's direction line;
  - the disabled state and its reason.
  Do not keep a second copy of the rule.

### Apply

- After the existing batch Python validation, compute orientations from the
  **validated** response with `matchedOrientations()`.
- Call `runAnalysis` once. Its `canonicalStateOverride` contains the validated
  plan, the translations, and those orientations.
- The automatic resolved path passes current orientations; nothing changes.
- The summary's `reversed` count reports how many records this Apply flipped.
- Retain the draft, including the checkbox state, on failure.

### One orientation input per generation

In `gbdraw/web/js/app/run-analysis.js`, apply
`canonicalStateOverride.linearRecordOrientations` once at the entry of the run,
as a run-local value.

- Every orientation reader reached by that run must use it:
  - the comparison display projection (`buildRegionSpec()` and
    `getViewTransform()`);
  - canonical serialization, including any injected serializer that reads
    state directly;
  - request construction;
  - any other `region_reverse` reader in the run.
- Delete the late patch that rewrites only the serialized `linearSeqs`.
- Commit the orientations through the existing generated-artifact owner set on
  success. Leave state untouched on failure, cancel, or stale completion.
- Add a unit test proving that the projection and the serialized request
  receive the same overridden orientation.

### Manual Reverse keeps the plan

- Remove `setManualOrientation()` from the controller.
- Remove `setLinearRecordOrientation()` and `linearRecordOrientationValue()`
  from `app-setup.js`.
- Bind the record checkbox in `index.html` with `v-model="seq.region_reverse"`,
  as on `origin/main`.
- No `Alignment cleared` notice appears. Plan validation before the next
  Generate re-aligns the same anchors in the new orientation.
- Keep the other clearing triggers: record drag, crop, selector, and source
  replacement.

### Reset Align

- `resetAlignment()` renders with no plan and the current base translations,
  and passes no orientation override.
- Notice: `Alignment reset: record positions restored; record directions unchanged.`
- Remove the orientation parts of `baseline()` and `installBaseState()`.
  Record drag still materializes translations only.

### Review UI (`gbdraw/web/index.html`)

- Below the exact-reference card, add the single **Match reference direction**
  checkbox with an accessible name and its status line:
  - checked: `Apply reverses N record(s): <record labels>.`
  - unchecked: `Record directions stay unchanged.`
  - any selected target with an unknown strand:
    `Unchanged because a strand is unknown: <record labels>.`
- When disabled, show `All selected anchors already face the reference direction.`
- In each target row, show one direction line for the selected candidate:
  `Direction: same as reference`, `Direction: opposite to reference` (plus
  ` — reversed on Apply` when the option is on), or
  `Direction: unknown strand — unchanged`.
- Keep keyboard access, focus handling, and 390 px reachability.

## Tests

Add unit tests in `tests/web/similarity-alignment-actions.test.mjs` for:

- `matchedOrientations()` with opposite, same, unknown, skipped, and reference
  records;
- local edits start no Worker job;
- the disabled reason;
- Apply passes the matched orientations and commits one History entry;
- failure, cancel, and stale completion leave orientations unchanged;
- Reset passes no orientation override;
- a manual Reverse keeps the plan.

In `tests/web/similarity-alignment-ui.playwright.spec.js`, add the
real-rendering regression journey on the Gallery example
`BGC0000708-BGC0000713`, with comparisons enabled:

1. Reference `CAG38712.1`, then **Review alignment options…**, then enable
   Match reference direction.
2. Apply succeeds and the *Streptomyces fradiae* ATCC 10745 record is reversed.
3. Its comparison ribbons are present, and the anchor centers match the
   reference within 0.5 px.
4. Reopening the review shows `Direction: same as reference` for that record.
5. Undo restores orientation and positions.
6. Reset Align restores positions and keeps the reversal.
7. A manual Reverse keeps the plan and re-aligns after Generate.
8. Save and load round-trip the result.

Check both Playwright installations; use Python Playwright if Node's runner is
missing, and rerun a Chromium sandbox failure with escalation.

Then run:

- the focused Node and Python tests;
- `node tests/web/architecture-contracts.test.mjs`;
- `node tools/check-web-change-budget.mjs --base origin/dev`;
- `git diff --check`.

Inspect desktop and 390 px screenshots of the review. Review production and
test diffs separately. Commit on the fix branch and push only to its
same-named remote branch.

## Handoff

Report:

- commits and changed files;
- test results and the measured anchor-center offsets;
- screenshot observations;
- the line-count delta for production code;
- owner/path evidence.

At the end, print the entire contents of
`docs/internal/similarity-alignment-orientation-owner-20260926/SESSION_04_DOCS_GALLERY_ACCEPTANCE.md`
in a copyable code block as the prompt for the next session.
