# Issue #602 final-review Status corrections

Local follow-up on 2026-09-27 (Asia/Tokyo). F1 and F2 are corrected in the existing
projection/import owners. This record supplements S00–S07; their historical
receipts and validation records remain unchanged. No commit or remote publication
is part of this task. Verification is bound to the four exact source hashes in the local evidence receipt.

## Candidate, base, and remote boundary

- Dedicated worktree: `/tmp/gbdraw-issue602-s00-results-20260926`.
- Branch: `fix/issue-602-linear-live-edit-20260926`.
- Upstream: `origin/fix/issue-602-linear-live-edit-20260926`.
- HEAD and fetched same-named remote HEAD:
  `9a14db3af8fafd6f8374cfe2a8197a8927cf43f7`; ahead/behind **0/0**.
- Reviewed base and current merge-base:
  `c922fc38aac78da9be83342c09ac0164ecef6ff6`.
- Fetched current `origin/dev`: `f5f86634459e0dcd46c1a452e9219fbba635d429`.
  HEAD/dev ahead/behind **18/15**; neither tip is the other's ancestor.
- Fetch discovered Issue #600 integration (PR #617), 89 changed files from the
  reviewed base. In the two affected owners, it adds annotation-warning Session
  persistence and Circular slot gap validation/projection, outside this correction.
  Shared History/watchers and editor changes also require preservation during
  later integration. No synchronization merge was made.
- Read-only remote branch query matched these fetched tips; PR lookup returned
  `[]`. Local edits are not present on the remote branch.
- Existing branch/history and the separate dev synchronization merge
  `52fda2061111ff7d1c6f69b46792a4a68a0ad5fc` remain intact. Shared checkout and
  existing stashes are preserved. No reset, branch recreation, checkout switch,
  commit, push, PR creation/edit, merge, deployment, or tag was performed.

## Authority and preflight

Classification: **IMPLEMENT_EXISTING_AUTHORITY**. PD-OI-037 / OIC-024 already
require truthful generation/live/History/Session observations and uncertainty
without Status-only Worker, genome reads/hash, or SVG/checkpoint cloning.
The correction realizes that complete existing outcome. Evidence exists for the
reported defects; no competing Product outcome remains, so neither
EVIDENCE_REQUIRED nor PRODUCT_DECISION_REQUIRED applies. The implementation
introduces no prohibited parallel owner, schema, reader, or admission path.

Formal Contract revision **22**, SHA-256
`5451db7da8e4de970f05afd5d19b003b7de77d52eb84c05d71ff85a533d4b1fb`,
is byte-identical in HEAD, the worktree, the reviewed base, and latest fetched dev.
The audit rechecked the four retained Issue #602 receipts (36 fields), all five
original receipts in historical S03 authority (45 fields), and the four Issue
#598 source receipts (36 fields). PD-OI-027/029/031/034 revisions **5/3/5/5**,
PD-OI-035 revision **3**, PD-OI-038, and PD-OI-039 revision **2**,
`A / EXCLUSIVE_DIRECTIONS_WITHOUT_MATCH`, retain their independent AND
contributions. Historical Match receipts remain historical. No authority bytes
were edited. PD-OI-046/047 are retained; Issue #601 runtime is outside scope.
Read-only audit: `/tmp/issue602-status-fixes-audit.py` and `.log`.

## Causes, owners, and resulting behavior

| Finding | Cause | Correction and owner |
| --- | --- | --- |
| F1: restored live width falsely Pending | Live commits update the artifact's in-memory intent, but import rebuilt that intent from the original render request alone. | `config.js::importSession` passes normalized saved editor evidence to `session-request.js::projectAppliedGenerationIntent`. Only explicit global block width is recovered from `originalSvgStroke.width`; saved draft alone never becomes applied. Missing/invalid evidence produces width-specific uncertainty. Automatic width remains request-owned; zero is valid. |
| F2: undiscovered selector falsely Pending | A placeholder selector was compared as a proven difference, then generic uncertainty was appended. | `projectGenerationIntent` identifies unresolved selection domains; `compareGenerationIntent` and its existing recursive comparator exclude those domains from definite differences while retaining known differences. Single-record selection/key is uncertain; undiscovered multi-record membership is uncertain. Explicit resolved changes and invalid selectors remain distinguishable. |

A subsequent proven live commit removes only its own uncertainty through the
existing `projectAppliedGenerationFields`; it cannot promote other unknown
fields. Scale, source, and other independent proven changes still report Pending.
Circular fresh Load may correctly show **Unknown** while discovery is deferred;
this correction does not force Applied or launch discovery to obtain it.

Owners/paths before and after: one canonical request projection and comparator in
`session-request.js`; one applied artifact handle and import coordinator in
`config.js`; the existing live SVG action supplies committed fields. Session
reader, Result admission/sanitizer, History and Worker are unchanged. The former
request-only import assignment and unconditional unknown-domain comparison are
replaced in place, with no parallel fallback left behind. There are no new
modules, exports, reactive stores, watchers, resource signals, persisted keys,
compatibility namespaces or privileged importers. Ordinary non-increasing
architecture evidence applies; no OE/PE/CB exception condition is triggered.
Rollback of this follow-up consists of reverting its two production-file edits
and accompanying regressions; retain authority and independent dev merges.

## Regression evidence

New browser tests extend `history-generated-authority.playwright.spec.js` using
real Gallery Sessions, Generate, the public Block Stroke Width control, Save
Session downloads, and fresh contexts with public file upload. Before the runtime
fix, all **3 cases failed** at the reported differences (Linear width short/long;
Circular selector plus widths; independently unedited Circular selector).
After the initial correction, all **3 passed**, **40.9 s**. The final extended
suite also checks genuinely Pending scale save/restore and failed-import rollback.
Baseline: `/tmp/issue602-status-fixes-before-browser.log` and report JSON.
Initial passing run: `/tmp/issue602-status-fixes-target-browser.log` and report JSON.

The new Node regression matrix in `session-request.test.mjs` covers both modes,
missing/invalid width evidence, zero/automatic values, draft width distinct from
artifact width, independent scale Pending plus unknown, resolving width
uncertainty with a later live commit, explicit selector changes, invalid selectors,
unknown discovery, multi-record uncertainty, and source replacement. Existing
Status spies throw on Worker construction, Blob reads, base64 conversion, digest,
and SVG clone. No timeout, checker, hit oracle, authority or existing acceptance
assertion was weakened. The fixture explicitly selects the single-record journey;
its first draft incorrectly left the normal grid default enabled and was corrected
without changing production semantics.

Mandatory SVG sanitization on import is compared against the actual sanitized
saved Result. Request identity and restored width remain exact; fresh import
constructs no diagram Worker. Failed import restores the exact canonical request,
resource table and applied-intent handles, draft, History counts and block width.
Existing Undo/Redo, Generate failure/cancel/stale, retry and Result/draft tests are
rerun against this comparator, rather than inheriting their older successes.

## Commands and final results

Environment: Node **26.8.2**, Python **3.13.3**, Python Playwright **1.61.0**,
shared Node Playwright **1.61.1**. Node root lookup was unavailable in the dedicated
worktree; shared absolute CLI and NODE_PATH resolved it. Sandbox commands could
not start because of the WSL bubblewrap mount configuration; the same authorized
commands ran with required escalation. Existing services were preserved; test-owned
ports 47430–47435 were used. Generated browser wheel/Python rendering inputs are
unchanged. Long checks were monitored incrementally without reducing timeouts.

| New check | Result / log |
| --- | --- |
| Focused projection/History/Status Node | **5 passed**, 349 ms; `issue602-status-fixes-focused-node.log` |
| Fast Web CI selection | **577 passed**, 23.54 s; `issue602-status-fixes-fast-web.log` |
| Status/History/Save/Load full browser suites | **18 passed**, 385991 ms, zero skipped/unexpected/flaky; `issue602-status-fixes-browser.log`, report JSON and traces |
| Lazy-import/rollback/divergent-draft Session contracts | **4 passed**, 92643 ms, zero skipped/unexpected/flaky; `issue602-status-fixes-session-browser.log`, report JSON and traces |
| Existing non-slow Python browser gate | **38 passed**, 6258 deselected, 318.89 s; `issue602-status-fixes-python-browser.log` |
| Architecture | **139 passed**, zero failed/skipped, 120441 ms; `issue602-status-fixes-architecture.log` |
| PR smoke | **13 passed**, 128743 ms, zero skipped/unexpected/flaky; `issue602-status-fixes-pr-smoke.log`, report JSON and traces |
| Authority/receipts/ancestry and retained evidence hashes | **PASS**; `issue602-status-fixes-audit.log` |
| Whitespace and tracked generated-artifact diff | **PASS** / no changes |

All logs above are under `/tmp`. Browser saved Sessions, attachments and traces
are under `/tmp/issue602-status-fixes-{browser,session-browser,pr-smoke-browser}`.
The three added round-trip cases are included in the 18, not counted twice.

Reproduction commands (same existing helpers/configurations):

```sh
node --test tests/web/session-request.test.mjs tests/web/history-canonical-owner.test.mjs tests/web/generation-status.test.mjs
node --test tests/web/architecture-contracts.test.mjs
NODE_PATH=/mnt/c/Users/genom/GitHub/gbdraw/node_modules GBDRAW_WEB_TEST_PORT=47432 node /mnt/c/Users/genom/GitHub/gbdraw/node_modules/@playwright/test/cli.js test tests/web/history-generated-authority.playwright.spec.js tests/web/generation-feedback.playwright.spec.js tests/web/session-loading-feedback.playwright.spec.js tests/web/session-save-lifecycle.playwright.spec.js --config=playwright.config.js --workers=1 --trace=on --reporter=list,json
NODE_PATH=/mnt/c/Users/genom/GitHub/gbdraw/node_modules GBDRAW_WEB_TEST_PORT=47433 node /mnt/c/Users/genom/GitHub/gbdraw/node_modules/@playwright/test/cli.js test tests/web/contracts/current-session-lazy-materialization.playwright.spec.js tests/web/contracts/session-regenerate-intent.playwright.spec.js --grep 'synthetic current session restores|preflight and lazy-access|no-draft session preserves|divergent draft and direct' --config=playwright.config.js --workers=1 --trace=on --reporter=list,json
NODE_PATH=/mnt/c/Users/genom/GitHub/gbdraw/node_modules GBDRAW_WEB_TEST_PORT=47435 node /mnt/c/Users/genom/GitHub/gbdraw/node_modules/@playwright/test/cli.js test --config=playwright.pr-smoke.config.js --workers=1 --trace=on --reporter=list,json
PATH=/mnt/c/Users/genom/GitHub/gbdraw/node_modules/.bin:$PATH NODE_PATH=/mnt/c/Users/genom/GitHub/gbdraw/node_modules GBDRAW_WEB_TEST_PORT=47434 python -m pytest tests/ -m 'browser and not slow' --durations=30 -v
python /tmp/issue602-status-fixes-audit.py
node tools/check-web-change-budget.mjs --base HEAD
node tools/check-web-change-budget.mjs --base origin/dev
WEB_ARCHITECTURE_CHANGE=true node tools/check-web-change-budget.mjs --base c922fc38aac78da9be83342c09ac0164ecef6ff6
git diff --check
```

Fast Web uses the unchanged CI selection: top-level `tests/web/*.test.mjs`,
excluding `architecture-contracts.test.mjs` and `gallery-session-publication.test.mjs`.
Browser evidence commands additionally supplied `PLAYWRIGHT_JSON_OUTPUT_NAME` and
`--output` paths as listed in the receipt. The absolute Node CLI avoids selecting
Python's separately installed `playwright` command.

## Reused evidence and its boundary

Every retained log hash was checked against the prior final-review JSON. S07's
**6268 passed / 17 skipped** remains evidence for unchanged Python core, rendering,
read-only reference output comparisons, recipes and Gallery code/inputs. It is
not reused to certify the modified Session/Status/browser paths. Those are rerun.
S07's **25 documentation contracts**, source/figure provenance from the native
raw-input reference journey, and Ruff apply to unchanged public docs, figures,
capture/source fixtures and Python inputs. That historical journey's Status
observations do not certify this modified comparator.
The prior build verifies unchanged packaging mechanics; it is not a newly built
artifact containing these uncommitted JavaScript edits.

S06 Editor **28 passed** and the presentation-only alignment/Reset portions of
S07's **21 passed** remain evidence for unchanged layout, direction, and geometry
owners. They do not replace new comparison/Session/History checks. S06's older
fast-Web **577** is superseded here by a new run. PR smoke is rerun separately from the
Python browser marker selection. No all-suite, all-staging, or release result is inferred.

## Independent diff audits

- Production: two existing files, **49 additions / 13 deletions**, net **+36**.
  Reviewed proof scope, immutable request/resources, unknown-domain behavior,
  effective defaults, subsequent live resolution and no additional side effects.
- Tests: two existing owners, **204 added lines**, reviewed separately. Public
  round trips, negative evidence/selector cases, rollback, sanitizer and side-effect
  assertions complement the implementation; they do not simply mirror it.
- Docs: this internal result record only. S00–S07, Product authority, public docs
  owner `docs/REFERENCE/web-app.md`, public figures and Gallery instructions remain
  unchanged. The PR body stays in `/tmp`, not in the repository.
- Generated artifacts: no tracked generated changes. Reference outputs, Gallery,
  vendored assets, social preview, dist and egg-info are preserved; wheel is excluded.

## Gate, required CI, and maintainer review

Local correction (`--base HEAD`): **Gate PASS / Review REQUIRED**, Session behavior
risk, no blocking violations or architecture-inventory delta. Actual cumulative
Issue #602 scope against merge-base `c922fc38`: **Gate PASS / Review REQUIRED**,
architecture profile, **7 production files**, **612 additions / 87 deletions**,
gross **699**, net **+525**. The net threshold remains exceeded (525 > 400).
Latest-dev direct tree comparison also returns **PASS / REQUIRED**; its 27 files,
gross 2465/net 1065 include inverse upstream Issue #600 differences and must not
be presented as this follow-up's authored scope or as permission to remove them.

Required CI from actual cumulative diff remains `web-change-budget`, `core-pr`,
`recipes-standard`, `gallery`, `lint`, `web-contracts-pr`, `web-pr-smoke`.
`tools/ci-impact-policy.mjs::classifyChanges` / `requiredJobsFor` are unchanged.
Trusted admission remains `Web base policy (trusted base)` and `PR / gate`.
Reports are local evidence, not remote CI successes for an unpublished candidate.

Maintainer **Review REQUIRED** remains separate from Gate PASS:

1. Applied-width evidence scope, effective/unknown semantics and saved Result/draft
   separation across import, History and live commits.
2. Shared request projection/artifact responsibility and cumulative net additions;
   no new owner/compatibility path or scientific identity regression.
3. Joint compact Editor/review, direction/Reset, focus, scrolling and export
   guarantees retained from S00–S07 and Issue #598.
4. Reconciliation with newer dev Issue #600 changes, followed by validation of the
   actual committed candidate and required PR CI. Current local results do not
   claim a tested integration with that newer base.

`architecture-change` routes review and selects its size profile; it waives no
check. This agent's diff audit does not fulfill human maintainer approval.

## Remaining limits and publication conditions

Inherited limits remain: Legend override live-rerender retry binding constraint;
223 px settings width at 195×422; CSS-viewport-equivalent 200% and focused resize
keyboard evidence; sequential page/list recovery on short screens; unexplained
S06 initial context destruction, no recurrence in S07; parent-dev Issue #564
metadata-free Session two cases unresolved. Issue #601 runtime, complete dev
staging and release readiness are outside this completion claim.

Before publication/integration: separately authorize commit/push/PR creation;
reconcile latest dev without losing Issue #600 or this correction; inspect actual
remote state; publish only the same-named work branch; validate the resulting exact
candidate, obtain required CI and maintainer review. No remote approval applies
to these uncommitted local changes.

Prepared English PR files: `/tmp/issue602-pr-title.txt` and
`/tmp/issue602-pr-body.md`.

Title: **Clarify Linear editing and preserve Status across Session loads**.

Opening paragraph:

> This PR explains Linear Auto label visibility and edit timing, starts fresh Web diagrams with a locked Definition column, and docks the existing Editor or alignment review below compact previews. Save/Load keeps applied live block widths distinct from pending generation settings; an unresolved Circular record selection shows Unknown instead of a false Pending change.

The body describes completed behavior, required CI, new and reused validation,
current-dev reconciliation, required human review and material inherited limits.
F1/F2 are no longer listed as unfixed. It follows the existing PR template and
`write-clear-pull-request` skill; final language checker output is saved in
`/tmp/issue602-status-fixes-pr-language.log` (**PASS**). The exact wording was frozen before
running that check once. Source hashes, evidence-log hashes, browser statistics,
PR-file hashes and final tree state are in
`/tmp/issue602-status-fixes-evidence.json`.

Proposed commit title: **Fix applied edit Status after Session loads**.

Summary: Restore proven live block widths through the existing artifact owner,
keep unresolved Circular selection uncertain, and verify Pending, History and
Session recovery without additional Status work.
