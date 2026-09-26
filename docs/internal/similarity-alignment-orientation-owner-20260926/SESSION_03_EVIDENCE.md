# Session 03 — record orientation workflow evidence

Branch: `fix/similarity-alignment-orientation-owner-20260926`.
Upstream and authorized push target: `origin/fix/similarity-alignment-orientation-owner-20260926`.

## Authority and scope

`origin/dev@4d1cf93514d0f75fa7a0ee32c1c4b2f176e03682` is an ancestor of this branch.
It includes authority commit `61eec6c5`, Product contract revision 16, and the
accepted `PD-OI-027/028/029/031/034` receipts. Session 01 is `5bcdd4f5`;
Session 02 is `92880e99`. No merge, rebase, or branch recreation was needed.

Preflight: `IMPLEMENT_EXISTING_AUTHORITY`. The merged receipts select the
single review option, record-owned orientation, manual Reverse retention,
position-only Reset, and atomic History and failure behavior. No unresolved
Product outcome, new compatibility reader, or authority change was introduced.

## Owners and paths

- Python remains the only owner of candidate eligibility, strand facts, and
  batch validation. `matchedOrientations()` in `similarity-alignment.js` alone
  applies the known-opposite reversal rule. Review availability, reversal
  labels, unknown labels, and direction lines derive from these same facts.
- `run-analysis.js` takes one run-local sequence orientation input at entry.
  Comparison projection, record display rows, the injected file serializer,
  and canonical request construction all read that input. The late serialized
  `linearSeqs` patch was removed. The real Gallery check exposed a live
  `recordDisplayRows` getter that still read old flags; its run-local projection
  is now covered by the unit regression.
- The existing generated-artifact owner set admits orientations, translations,
  plan, and Result through the existing History replacement transaction.
  Failed, canceled, and stale runs have no orientation commitment.
- The alignment controller no longer owns manual Reverse. Its
  `setManualOrientation()` and the two app-setup orientation adapters were
  removed; the checkbox uses `v-model="seq.region_reverse"`.
- Reset and record-drag baselines retain translations only. The orientation
  baseline and its installation, plus `requestWithOrientations()`, were removed.

The run's competing live/late-patch orientation paths converge on one input.
No additional semantic owner or canonical production path remains; OE and PE
are non-increasing (the run paths decrease), and CB is unchanged. This is the
ordinary concise ratchet case; no exception is required. The change-budget
checker reports Gate PASS and Review REQUIRED, with no blocking violation.
Review signals include the earlier Session 02 changes relative to `origin/dev`.

Session 03 production delta: **+22 lines** (124 added, 102 removed) across four
existing files. The new global option, derived direction/status view, and native
checkbox markup account for the growth after Session 02 had already removed
per-target controls. No module, dependency, persisted field, compatibility
branch, or Worker lifecycle was added.

## Changes reviewed separately

Production:

- `gbdraw/web/js/app/similarity-alignment.js`
- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/web/js/app/app-setup.js`
- `gbdraw/web/index.html`

Tests:

- `tests/web/similarity-alignment-actions.test.mjs`
- `tests/web/run-analysis-simple-path.test.mjs`
- `tests/web/similarity-alignment-ui.playwright.spec.js`

Current documentation corrections:

- `docs/REFERENCE/web-app.md`
- `docs/RELEASE_NOTES_0.14.0.md`

Generated artifacts: none changed. Existing Gallery input was used unchanged;
tracked reference outputs and the owner-maintained social image were untouched.
The initially untracked `docs/internal/WEB_GUI_DEV_BEHAVIOR_AUDIT_2026-09-26.md`
was excluded from this task and from staging.

## Geometry and visual inspection

Real Chromium rendering used the Gallery session `BGC0000708-BGC0000713`,
reference protein `CAG38712.1`, with comparisons enabled. Four selected targets
had opposite strands, so Apply reversed all four and reported `4 reversed`.
The exact *Streptomyces fradiae* ATCC 10745 record (`record-2`) reversed.

After Match, measured horizontal offsets in CSS pixels were:

| Target | Offset |
| --- | ---: |
| record-2 | 0 |
| record-3 | -0.0001220703125 |
| record-4 | -0.0001220703125 |
| record-5 | +0.0001220703125 |

There were **77 comparison ribbons**, including **40 attached to record-2**.
After manual Reverse of record-2 and Generate, the offsets were respectively
`0`, `-0.0001220703125`, `-0.0001220703125`, and `0`, with the same ribbon counts.
Both maxima are below the required 0.5 px.

The browser journey also proves reopened `same as reference` with the Match
option off and disabled, one History entry for Apply, exact Undo/Redo of Result
and orientations, Reset restoring the original five record-group transforms
while keeping the reversals, manual Reverse retaining the same plan and anchors,
and a reversed result's Save/Load round trip in a fresh app.

Final desktop and 390 x 740 review screenshots were visually inspected. The
single checkbox follows the exact-reference card; reversal labels wrap; target
direction text is readable; the footer and Apply/Cancel controls remain visible;
there is no horizontal overflow. Space toggles the checkbox and Escape dismisses
the review. Existing Gallery definitions contain literal `<i>` markup, which
was already displayed as text in review record labels; it remains visible in
the status list as well. This session does not alter label formatting.

Artifacts and command logs are retained under `/tmp/session03-*`; final review
screenshots and the saved Session are under `/tmp/session03-accepted-browser-results`.

## Verification

Both Playwright installations were checked: CLI 1.61.0, Python Playwright
available, and Node `@playwright/test` resolvable. Node Chromium ran the real
journeys. The workspace sandbox failed to initialize with the known host-mount
error, so local reads, edits, and checks used the reviewed sandbox escalation.

- `node --test tests/web/similarity-alignment-actions.test.mjs
  tests/web/session-request.test.mjs tests/web/run-analysis-simple-path.test.mjs
  tests/web/gallery-session-publication.test.mjs`: **41 passed**.
- `python -m pytest tests/test_similarity_alignment.py
  tests/test_similarity_alignment_rendering.py
  tests/test_similarity_alignment_web_adapter.py tests/test_session_request_codec.py
  tests/test_session_compat.py tests/test_api_session.py
  tests/test_documentation_reference_contracts.py
  tests/test_gallery_session_semantics.py -q`: **304 passed**, 10 existing
  Biopython qualifier-length warnings.
- `python -m pytest tests/test_documentation_reference_contracts.py -q` after
  current documentation corrections: **12 passed**.
- `npx playwright test tests/web/similarity-alignment-ui.playwright.spec.js
  --workers=1 --output=/tmp/session03-accepted-browser-results`: **8 passed**.
- `npx playwright test tests/web/linear-multi-record.playwright.spec.js
  --grep 'Linear region controls|No comparison completes a real render'
  --workers=1 --output=/tmp/session03-linear-results`: **2 passed**.
- `node tests/web/architecture-contracts.test.mjs`: **137 passed**.
- `node tools/check-web-change-budget.mjs --base origin/dev`: **Gate PASS**,
  **Review REQUIRED**, no blocking violations.
- `git diff --check`: passed.

The first real-render check found the display-row reader mismatch described
above; the corrected regression passes. A full alignment browser run also
caught an existing test reading the error box immediately after viewport
resize. It now awaits the owner's resize clamp and still asserts that the
entire box fits within 390 px; the complete final spec passes.

Session 04 remains responsible for the public/Gallery consistency sweep,
regeneration round trip, exports, and the integrated full-test/build gates in
its instruction prompt. Session 03 introduces no public Gallery replacement.

Commit title: `Match similarity alignment directions through record orientation`

Summary: Add one local Match reference direction option, pass validated target
orientations through every generation reader, preserve the plan on manual
Reverse, and restore positions only on Reset. Cover atomic History, rollback,
comparison ribbons, exact anchor alignment, and Session round trips.
