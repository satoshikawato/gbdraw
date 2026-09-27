# PR #618 — synchronization after concurrent dev advancement

Author / Product Decision Owner: `satoshikawato <kawato@kaiyodai.ac.jp>`.
Date: 2026-09-27 (Asia/Tokyo). Same dedicated clone and work branch as
[original review](REVIEW_20260927.md). This report does not contain its own
containing commit SHA; final candidate identity is reported in the handoff.

## Publication, initial green CI and strict-base boundary

The owner approved publication and merge of the prepared original candidate
`731dfdff1fe6d3ff6f484c0a24bcff0d03938747`. PR
https://github.com/satoshikawato/gbdraw/pull/618 was created to dev with the
exact language-checked title/body. Published head, initial base, title/body
and all 23 paths matched local preparation. No existing same-head PR existed.

The original PR's required `Web base policy (trusted base)` and `PR / gate`
both completed SUCCESS. Core PR, Recipes, Gallery, Web contracts, Web PR smoke,
Lint, Web change budget and CodeQL also completed SUCCESS; other jobs were
properly skipped by PR routing. Tests run:
https://github.com/satoshikawato/gbdraw/actions/runs/36287108981 . Trusted-base
run: https://github.com/satoshikawato/gbdraw/actions/runs/36287108658 . Its log
reports Gate PASS / Review CLEAR and four HARD architecture rules CONFORMING.
These hosted results bind to the original head only.

Dev protection requires both named statuses with `strict: true`, no forced
pushes and no approving-review count. During the CI wait, separate approved
PR #617 merged into dev:
`f5f86634459e0dcd46c1a452e9219fbba635d429`. PR #618 became BEHIND despite green
checks. No admin bypass, protection change or merge of the old candidate was
attempted. Remote status observations were more than five minutes apart.

## Inspected dev changes and authority

The c922fc38→f5f86634 diff has 89 paths (5170 additions/1309 deletions), all
from the independently admitted #600 work. Its completed writer evidence was
read from Git and after incorporation:
`issue-600-implementation-20260926/SESSION_RESULTS/S05.md` and
`S05_GUARD_INTEGRATION_20260927.md`. Historical failing diagnostics in S05
are explicitly superseded by the guard-integration follow-up (including a
supervised 6642-pass native run and separate browser/Node/guard success).
Those are #600 evidence, not fresh executions by this session.

Reviewed intersections: run-analysis candidate validation now carries safe
annotation warnings; Generate color preparation uses a shallow candidate
state and native normalized captions; existing History checkpoints carry
rules/warnings and suppress semantic file watchers while restoring; Session
save/load retains warnings; prepared native annotation resolution and caption
normalization remain with their established owners. No #600 runtime was added
independently by this agent: only the already-admitted dev was normally merged.

The complete Product Contract is byte-identical across both bases, SHA-256
`5451db7da8e4de970f05afd5d19b003b7de77d52eb84c05d71ff85a533d4b1fb`.
Therefore all nine fields of PD-OI-027 r5, 029 r3, 031 r5, 034 r5 and 039 r2,
and the entire independent PD-OI-035 r3 section retain the exact original
review match. EXCLUSIVE_DIRECTIONS_WITHOUT_MATCH is not reopened. The four
required authority/dev ancestors remain ancestors of f5f86634. Alignment
controller, direction projection/receipt admission, capture owner, source
pins and restored alignment browser test remain unchanged by this dev delta.

Preflight: IMPLEMENT_EXISTING_AUTHORITY. No new unresolved material Product
outcome, domain-policy selection, authority change, compatibility migration
or BUG-17 action. Existing accepted conditions and all independent canvas,
keyboard/focus/Editor contributions remain necessary together.

## Normal incorporation and fresh scoped verification

`git merge --no-commit --no-ff origin/dev` completed without conflicts, retaining
all S00–S04 and original review history. Planned parents are original review
HEAD `731dfdff` and actual dev `f5f86634`. Working runtime equals f5f86634
byte-for-byte, excluding Gallery tutorial prose. No production change is made
relative to that admitted base. JS runtime tree:
`68a21e70be1dda423bcdaee7458bebd80ad1a5c1`.
The binary c922fc38→f5f86634 diff SHA-256 is
`b5ea7ee240bc4efa0d9fae8c4a3c7034cff9bbb634131905105d9379f186fc26`.

The old S04 broad results remain historical and do not certify all changed
runtime. Fresh executions target the actual shared-boundary delta instead:

| Current check | Result |
| --- | --- |
| Focused Node alignment/receipt/record/candidate/run-analysis/History/Session/palette owners | 118 PASS, exit 0, 3.54s |
| Native alignment/projection/adapter/rendering | 83 PASS |
| Documentation capture/reference contracts | 41 PASS |
| Read-only tracked SVG comparison | 16 PASS |
| Combined native/doc/SVG command | 140 PASS, exit 0, 75.12s |
| Current browser alignment + related History/palette + region | 16 PASS, exit 0, 5.6m |
| Ruff | PASS |
| Revised exact title/body language check | PASS, exit 0 |
| Working-tree actual-base Web policy | Gate PASS / Review CLEAR, ordinary, no blocker/review reason, four HARD rules CONFORMING |
| Actual raw-input T-GUI-04 | six real PNGs, 155 features; SVG 249702 bytes; TSV 14159 bytes/232 rows; direction/both Reset/Undo semantic checks PASS |
| Final exact screenshot check | PASS: all six committed images match a fresh capture, exit 0 |

Browser geometry remains logical: reference center before/after
1816.4087447808777 / 1816.4087625383818; sparse actual correction
870.6455831649031. Inverse screen CTM and automatic composition translation
are removed; no screen-coordinate promise is introduced. Comparisons retain
77 ribbons (40 target ribbons), maximum anchor offset 0.0001221 logical px.
Fresh-load both Reset scopes start zero LOSAT jobs. Full browser journeys also
cover source resource/strand preservation, pending form, rollback/retry,
History/receipt consumption, keyboard, compact canvas and Editor lifecycle.
Region tests cover the changed common History owner; old Rotate/other broad
acceptances are not relabeled as fresh executions.

## Screenshot failure diagnosis and finished artifacts

The first `T-GUI-04 --check` completed actual operations and SVG/TSV checks but
reported stale 03/04/05. This failure is recorded, not hidden as PASS.
The existing capture owner was rerun normally. Against S04, 03 changed 25640
pixels within comparison bands with maximum channel delta 4; 04 changed 587
pixels with maximum channel delta 2. Both updated captures were visually
inspected at 1440×900 and subsequently reproduced exactly. They retain the
finished five-record biological/presentation context.

One regenerated 05 differed in rendered figure pixels. The next strict run
and a read-only instrumented replay both produced a 05 byte-identical to S04.
The replay used the existing `_capture_scenario`/real GUI journey, saving
viewport observations and actual downloads; it injected no completed state or
private generation call. Final viewport was zoom 0.4, pan (-782,0), wrapper
x 367.45001220703125 / y 157, dimensions 1044.60009765625×332.6312561035156.
Final admitted SVG SHA-256:
`44b3a289cc2d0238f2fbbe9bd21761259feb628f1f6ccea991aee2f6ee3793bc`.
The differing capture is a diagnostic artifact; its underlying cause is not
claimed conclusively resolved as a runtime bug. Restore only 05 from the real
repeat capture output (equal to its S04 bytes), retain reproducible 03/04,
then run the unchanged strict `--check` again: all six match. Final PNG total
1057022 bytes. No comparison tolerance, raster exception, harness code,
bitmap hand edit or output baseline was used to force success. The observed
single divergent capture remains a bounded reproducibility limitation.

Changed PNG hashes:
- 03-losatp-settings.png: `95be2dabcbb1c30c8f03142070e0a092387da6f3703d87d0430c8e6e66a4bf43`
- 04-align-og1.png: `76ab73145c1ed5bf4a2036ca50f669a7d4b512b0bf42232ea8f7c366a1fdf057`

All other PNGs retain the original-review hashes. Reference SVGs, social
preview, Gallery source/session/thumbnail/manifest and vendored assets remain
unchanged. Generated browser wheel is ignored; no cache-bust refresh.

## Environment, commands and documentation skill

The same isolated .venv / clone-local npm / pinned Chromium environment from
original review is retained. No other checkout, environment, wheel or server
is operated. Fresh clone-local wheel command:
`.venv/bin/python tools/prepare_browser_wheel.py --no-build-isolation`.
Wheel SHA-256:
`c87e3805e7889f2c548a0e12469b3d3e5c39dec18aaa147848413a11fde1ed57`.
Node browser server used dedicated kernel-selected port 41477. Each existing
documentation harness run selects its own loopback port.

```sh
node --test tests/web/alignment-direction-projection.test.mjs tests/web/alignment-reset-receipt.test.mjs tests/web/similarity-alignment-actions.test.mjs tests/web/record-display-options.test.mjs tests/web/run-analysis-simple-path.test.mjs tests/web/candidate-render.test.mjs tests/web/history*.test.mjs tests/web/session-active-config-contract.test.mjs tests/web/session-request.test.mjs tests/web/session-authority.test.mjs tests/web/session-feature-metadata.test.mjs tests/web/palette-history-preservation.test.mjs
.venv/bin/python -m pytest tests/test_alignment_direction_projection.py tests/test_similarity_alignment.py tests/test_similarity_alignment_web_adapter.py tests/test_similarity_alignment_rendering.py tests/test_documentation_capture_contracts.py tests/test_documentation_reference_contracts.py tests/test_output_comparison.py::TestOutputComparison -q
PATH="$PWD/.venv/bin:$PATH" PYTHONPATH=. GBDRAW_WEB_TEST_PORT=41477 node node_modules/@playwright/test/cli.js test --config=playwright.functional.config.js tests/web/similarity-alignment-ui.playwright.spec.js tests/web/history-comparison.playwright.spec.js tests/web/history-generated-authority.playwright.spec.js tests/web/palette-history-preservation.playwright.spec.js tests/web/history-region.playwright.spec.js --workers=1 --retries=0 --output=/tmp/gbdraw-issue598-review-Gb6Qjp-review-artifacts/f5-browser
PATH="$PWD/.venv/bin:$PATH" PYTHONPATH=. .venv/bin/python docs/capture/run_all.py --scenario T-GUI-04 --tier extended
PATH="$PWD/.venv/bin:$PATH" PYTHONPATH=. .venv/bin/python docs/capture/run_all.py --scenario T-GUI-04 --tier extended --check
.venv/bin/python -m ruff check gbdraw/ docs/capture/flows/bgc_losatp.py
node tools/check-web-change-budget.mjs --base origin/dev
# After commit:
node tools/check-web-change-budget.mjs --base origin/dev --head HEAD
```

Raw logs and read-only diagnostic frames/exports are isolated under
`/tmp/gbdraw-issue598-review-Gb6Qjp-review-artifacts/`, with f5-* names.
Original failed checks retain their own logs. No test was externally timed
out or its test-owned timeout changed.

`love-me-love-my-docs` applied for this failed capture boundary. Existing S04
page decisions are reused: keep the English end-user Web reference, five-BGC
Tutorial and capture README; keep source acquisition/provenance, five original
GenBank inputs and full comparison figure. No new page/harness is introduced.
The existing reader-route authority is docs/DOCS.md and the Tutorial manifest.

Docs Progress:
- [x] Step 1: Frame — existing English end-user Web/Markdown owner
- [x] Step 2: Flow census — keep existing owners, no new pages or harness
- [x] Step 3: Demo data — same five checksum-pinned public MIBiG inputs
- [x] Step 4: Smoke proof — complete actual GUI recipe and real diagram
- [x] Step 5: Execution harness — existing committed T-GUI-04 owner unchanged
- [x] Step 6: Evidence run — real GUI/Reset/Undo, PNG/SVG/TSV checked
- [x] Step 7: Manual written — existing prose remains accurate; only two PNG refreshes
- [x] Step 8: Verify + report — original strict check PASS; regeneration command already documented

The existing S04 selector findings remain unchanged (six changed .locator
calls, eight in the full BGC file, plus its label-ancestor scroll); no new
selector or private state injection is added to committed automation. No auth
or seed is used. Base host, viewport, locale/theme and font versions stay with
existing config.py. The failed --check exposed why this existing capture
should also be exercised when shared rendering changes.

## Review, limits and approval scope

Production, tests, docs and generated diffs were reviewed separately. Compared
with actual new dev, production runtime remains zero-delta and alignment test
diff remains the original review digest. No architecture exception/label is
needed; OE/PE/CB scope is unchanged. Explicitly stage only this report, updated
PR_BODY.md and the two refreshed PNGs; already-admitted dev paths are staged by
the normal merge. Do not use git add -A.

S04 broad results are historical; no blanket unchanged-runtime reuse is claimed
across #617. The #600 writer's native broad verification is prior base evidence.
Current focused checks cover changed shared paths and this PR's restored browser
assertions. Whole current Python/Node/browser suites, every supported version,
slow/performance, physical zoom, live assistive-technology speech and release
qualification are not newly executed. Original green remote CI does not certify
the synchronized head; it gets its own CI. Final commit Gate, push and new-head
hosted results are reported with the final SHA externally.

Initial authorization still permits required local integration/verification/
commit and normal same-name branch push. The original publication/merge approval
was tied to 731dfdff and its exact body. The supplied AGENTS rule says:
“A new candidate does not inherit approval tied to an earlier candidate.”
The new candidate and revised local PR body therefore need a new explicit
approval before body publication and merge; do not treat the old green checks
or CLEAR as that approval. No deployment/tag/release/direct dev/main push.

English commit title: `Merge latest dev and verify Issue 598 shared workflows`.
Summary: Incorporate admitted annotation/style changes, verify alignment and
History continuity, refresh two reproducible screenshots and prepare updated
PR review evidence.

## Retained differential and source identifiers

Alignment browser diff SHA-256: `af7aed11a44a93be878d5e3e8d9001cfba5b374ac90c106a10cbde3b9c8521a6`.
Native test/reference and source baseline checks remain with the focused commands above.

Revised local PR body SHA-256: `20735d2d0c5d45922318e28739da7b6bad410fd4d0f9eb31a0b89b324399c53d`.
Language-check log: `f5-pr-language.log`; no revised-body publication yet.
