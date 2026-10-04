# Session 04 instruction prompt — documentation, Gallery, and integrated acceptance

Paste this entire file into a new agent session. It is self-contained.

## Mission

This session finishes and verifies record-owned orientation for gbdraw Linear
Similarity Group alignment. The intended result is:

- each record's Reverse complement setting is the only orientation owner;
- an alignment plan stores anchors only;
- the review offers one **Match reference direction** option that reverses
  eligible targets through record orientation;
- a manual Reverse keeps the active plan;
- Reset Align restores positions but not orientation;
- Undo restores everything;
- a failed generation shows its underlying error once.

The full behavior is section 3 of
`docs/internal/similarity-alignment-orientation-owner-20260926/01_MASTER_PLAN.md`.

## Branch and prerequisites

1. Use exactly `fix/similarity-alignment-orientation-owner-20260926`, with
   upstream `origin/fix/similarity-alignment-orientation-owner-20260926`.
2. Fetch `origin`. Confirm the merged authority for
   `PD-OI-027/028/029/031/034` (contract revision 16) and the Session 01–03
   commits.
3. Merge newer `origin/dev` if needed. Never rebase or recreate the branch.
4. Preserve unrelated working-tree content.

## Read

Read:

- `AGENTS.md`, `CLAUDE.md`, `gbdraw/web/CLAUDE.md`;
- `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md` and
  `docs/internal/PRODUCT_IMPACT_RATCHET.md`;
- the master plan;
- `.agents/skills/web-gallery-screenshot-maintenance/SKILL.md` before editing
  Gallery tutorials or screenshots.

Inspect every in-scope production, test, documentation, and generated diff
once, separately.

## Documentation

- `docs/REFERENCE/web-app.md`: describe the following:
  - the single Match reference direction option and its status line;
  - per-target direction lines;
  - that a manual Reverse keeps the alignment;
  - Reset Align without orientation restoration;
  - Undo.
  Remove the per-row wording, and remove manual orientation from the list of
  edits that clear alignment.
- `docs/RELEASE_NOTES_0.14.0.md`: describe the current behavior only. Do not
  mention unreleased intermediate plan fields.
- Search current public docs and Gallery tutorials under
  `gbdraw/web/gallery/tutorials/` for `Match reference direction`,
  `orientationPolicy`, `effectiveReverseComplement`, `match_reference`, and
  "orientation clears". Correct each current statement. Leave historical plan
  directories unchanged.
- Recapture a Gallery tutorial image only if its instruction shows the changed
  review. Use `tools/capture_gallery_tutorial_screenshots.py` and its strict
  check.

## Gallery

- Run `python tools/refresh_gallery_sessions.py` (unfiltered) if any render or
  Session byte changed since Session 02. Then run
  `python tools/gallery_artifact_manifest.py`.
- Confirm the BGC Gallery session contains no removed plan field.
- Render the BGC example and inspect it at readable scale.

## End-to-end acceptance

Using real rendering, demonstrate each item of master-plan section 3,
including the regression journey in section 8, and record measurements:

- **Anchor offsets:** after Match on the BGC example, the maximum
  anchor-center offset is at most 0.5 px.
- **Comparison ribbons:** the reversed *Streptomyces fradiae* ATCC 10745 record
  has comparison ribbons.
- **Idempotence:** a second review shows `same as reference`.
- **Undo, Reset, and manual Reverse:** each behaves as specified.
- **Session:** save and load, followed by regeneration, reproduces the figure.
- **Exports:** they contain no transient UI.
- **Default Align:** orientation is unchanged and the summary reports
  `0 reversed`.
- **CLI:** `--align_orthogroup_feature` output is unchanged.

Inspect desktop and 390 px review screenshots for legibility and reachability.

## Required gates

Run focused tests first and fix concrete failures. Then run:

- the focused Python and Node commands listed in master-plan section 8;
- `npx playwright test tests/web/similarity-alignment-ui.playwright.spec.js
  tests/web/linear-multi-record.playwright.spec.js --project=chromium
  --workers=1` when Node Playwright is installed; otherwise run equivalent
  Python Playwright checks;
- `ruff check gbdraw/`;
- `node tests/web/architecture-contracts.test.mjs`;
- `node tools/check-web-change-budget.mjs --base origin/dev`;
- `python tools/update_cli_reference_help.py --check`;
- `python -m pytest tests/ -v -m "not slow"`;
- `python -m pytest tests/test_output_comparison.py::TestOutputComparison -v`;
- `python -m build`;
- `git diff --check`.

Prepare the gitignored browser wheel with `python tools/prepare_browser_wheel.py`
when needed. Allow at least 30 minutes for the full pytest run and monitor it
incrementally. Do not update `tests/reference_outputs/` to silence a
comparison. Run an offline bundle audit only if dependencies, privacy, bundle
composition, or the Worker lifecycle changed.

## Architecture confirmation and record

Confirm that each of the following exists exactly once:

- **Orientation:** one owner, `region_reverse`, with one application path in
  the typed request and one run-local orientation input per generation.
- **Anchor resolution:** one Python resolver and validator.
- **Workflow:** one Web controller.
- **Persistence:** one request and Session writer.
- **Commit:** one Result and History admission path.

Confirm that none of the superseded code remains:

- plan orientation fields;
- late render reversal;
- the partial serialization patch;
- per-target controls;
- `setManualOrientation()`.

Provide concise before/after owner/path evidence, and full exception sets only
if a ratchet exception condition applies. Check every jointly required Product
effect, not only decision IDs.

Append the acceptance record to section 10 of the master plan:

- exact commits and the authority SHA;
- commands and results;
- measured geometry and visual observations;
- documentation, Gallery, and reference-output dispositions;
- remaining limitations;
- a proposed English pull-request title and summary.

Commit on the fix branch and push only to its same-named remote branch. Open
or merge a pull request only with explicit authorization for that step.

This is the final implementation session. If every criterion and gate is
complete, state that no further session is needed. Otherwise, write a
self-contained continuation prompt naming the exact remaining task, branch,
evidence, and next gate, and print it in a copyable code block.
