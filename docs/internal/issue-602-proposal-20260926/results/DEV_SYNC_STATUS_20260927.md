# Issue #602 dev synchronization and Status verification — 2026-09-27

This result follows `FINAL_REVIEW_STATUS_FIXES_20260927.md`. It preserves the historical S00–S07 records and distinguishes the locally committed fixes, the normal dev merge, and the subsequent integration corrections. No publication is authorized by this record.

## Local history and synchronization

- Dedicated worktree: `/tmp/gbdraw-issue602-s00-results-20260926`.
- Work branch: `fix/issue-602-linear-live-edit-20260926`; upstream: the same-named origin branch.
- Initial local/upstream HEAD: `9a14db3af8fafd6f8374cfe2a8197a8927cf43f7`.
- Prior reviewed merge-base: `c922fc38aac78da9be83342c09ac0164ecef6ff6`.
- Fetched dev and current comparison base: `f5f86634459e0dcd46c1a452e9219fbba635d429`. A second fetch confirmed the same dev and remote work-branch HEAD.
- Status fix commit: `7955c02f3eb9147566eb93ae44aa92b0b38ccbed` (`Fix applied edit Status after Session loads`). It contains exactly the four handed-off source/test files and their new historical result record.
- Normal merge: `11c81b1cea80caa8b27012cdaef19ec7e16cf830` (`Merge origin/dev into Issue 602 Status fixes`); parents are the Status fix and fetched dev. Both the initial branch HEAD and dev are ancestors.
- Integration correction commit: `f6d44139edcdc47fd3c84a967a873a94bfe8d160` (`Preserve caption edits and cold Generate feedback after dev merge`), containing five production and three test files. The final documentation commit is identified by the handoff receipt; this document does not attempt to embed its own commit ID.
- No branch recreation, reset, amend, rebase, push, PR creation/edit, remote merge, deployment or tag was performed. The read-only PR lookup returned no existing PR.

The merge conflict was confined to the end of `tests/web/session-request.test.mjs`. Removing the three conflict markers retained the complete Status/Lock test block and dev's complete optional-pixel-gap block. `git show --remerge-diff 11c81b1c` records this resolution. Production merges retained dev's annotation-warning persistence, numeric-only Circular gap validation and string UI projection, canonical Python color preparation, and atomic legend/History behavior.

## Corrections and owners

The original F1 correction restores only proven applied live stroke-width evidence through the existing artifact/import owner. Automatic width remains request-owned; zero is valid; missing proof remains Unknown for that domain. F2 excludes unresolved selection domains from definite differences while retaining independent real Pending differences. The request projection/comparator remains single, with no discovery/read/hash/clone solely for Status.

Two integration regressions were exposed by broader actual-browser checks:

1. **Complete default-caption recolor.** Public “Apply to all tRNA” rejected the default legend swatch as a caption collision. The existing color action now supplies the proven previous swatch to dev's canonical atomic legend transaction. It does not change caption normalization, collision rejection, the Python evaluator, or checkpoint History. Added 1440 px and 390 px tests failed before the production fix and passed after it, observing mounted and saved SVG colors, independent Pending differences, one History operation, Undo/Redo, Save/fresh Load and Generate. The existing negative collision/rollback case also passes.
2. **Real first-use Generate progress.** Dev's color-candidate helper initialized the shared Worker before the existing render progress callback, leaving “Preparing input files...” visible during actual cold initialization. The helper now forwards an optional progress observer from the same Generate operation. Runtime readiness is still owned by `diagram-generation.js`; `run-analysis.js` owns the existing message mapping and generation-token guard. No readiness state, Worker, render request or progress timer is added. The unchanged cold/warm/cancel/error/retry browser case failed before this correction and passes afterward. New client tests cover helper cancellation, cleanup, one shared warm render Worker and forwarding through both Python rule stages.

`rule-actions.js`, `history.js`, `history-snapshot.js`, `watchers.js`, and the legend layout owner retain dev's implementation. The integration observer is plumbing through existing owners. Semantic owners, canonical execution paths and persisted compatibility paths do not multiply: changed-scope OE/PE/CB do not increase. There is no exception, duplicate owner, hard-invariant waiver or new persisted compatibility namespace. The old inline Generate stage mapper was removed when its existing mapper was shared with preparation.

Product preflight is `IMPLEMENT_EXISTING_AUTHORITY`: truthful operation progress and the existing full-caption recolor continuation are preserved. No materially different product option is selected. Authority revision 22 and `PD-OI-037` / `OIC-024` govern Status; independent signed contributions remain AND requirements. `PD-OI-027/029/031/034` revisions 5/3/5/5, `PD-OI-035` revision 3, `PD-OI-038`, and `PD-OI-039` revision 2 (`A / EXCLUSIVE_DIRECTIONS_WITHOUT_MATCH`) remain intact. PD-OI-046/047 and Issue #601 runtime are untouched.

## Verification and provenance

The final runtime/test byte inventory is `/tmp/issue602-dev-sync-final-runtime-source.json`. Actual commands, exit receipts, JSON browser reports and traces remain under `/tmp/issue602-dev-sync-*`; final hashes are collected in the handoff evidence receipt. No failed or interrupted process is counted as a successful run.

| Final check | Result | Actual evidence |
| --- | --- | --- |
| Root Web Node contracts | 988 passed, exit 0 | `fast-web-progress-final.log`; includes the 848 CI-selected fast contracts, 139 architecture contracts and one Gallery publication contract |
| Architecture contracts | 139 passed, exit 0 | `architecture-progress-final.log`; two mistyped additional file patterns select no extra tests, so this log is not called a 192-case run |
| Architecture/Product ratchet fixtures | 53 passed, exit 0 | `ratchet-final.log`; also included in the root Node run |
| Focused Worker/progress/rule/resource contracts | 31 passed, exit 0 | `progress-focused-final.log` |
| Status, History, lazy Session and regeneration | 25 passed, no skips/retries | `status-session-progress-report.json`, command/exit receipt and traces; excludes the separately recorded failing structural case |
| Caption collision protection, new recolors and Save/Load lifecycle | 6 passed, no skips/retries | `color-annotation-lifecycle-final-report.json`, command/exit receipt and traces |
| Annotation/style integration, Circular/Linear at 1440/390 px | 4 passed, no skips/retries | `annotation-final-report.json`, command/exit receipt and traces |
| Official PR smoke | 13 passed, no skips/retries | `smoke-final-report.json`, `smoke-python-final-exit.json` |
| Python browser selection | 38 passed, 6,632 deselected, exit 0 | `python-browser-final.log`, `smoke-python-final-exit.json` |
| Public scenario/reference/inventory contracts | 25 passed | `doc-contracts-final.log`; unchanged inputs and tests remain valid after JS-only integration fixes |
| Ruff, authority, whitespace and preservation | PASS | `final-ruff.log`, `final-authority-audit.log`, preservation receipts and final handoff audit |

All abbreviated evidence names in this table have the prefix `/tmp/issue602-dev-sync-`. The environment is the same `/home/kawato/micromamba/bin/python`, Python 3.13.3, pytest 9.0.2 and recorded plugin versions as dev's supervised run. Shared Node Playwright and Python Playwright are available; Chromium checks used the required sandbox escalation. No test timeout was relaxed.

The initial post-merge browser run passed all 18 annotation/alignment cases. The later progress plumbing changes only observed first-use preparation, so the four Generate-containing annotation/style cases are run again on final bytes. The unchanged alignment presentation and geometry owners retain the original 14-case evidence; this does not label the entire earlier 18-case run a final-source run.

### Reused evidence

`/tmp/issue602-dev-sync-reuse-fast-audit.json` and `/tmp/issue602-dev-sync-reuse-audit.json` compare 667 native-production/test/fixture/example/dependency blobs with dev's actual supervised Python source `cb2ad027d3f60fde2331a23ad8f9f6274fe863ba`. Its start receipt, exit receipt and log exist and agree: **6,642 passed, 17 skipped, 11 deselected**, exit 0. The supervised node-ID coverage proves 189/189 recipe and 103/103 Gallery cases, including read-only reference comparisons. This is reused evidence, not a fresh 6,642-case run here. Modified Web, Session, History and progress paths do not inherit that success.

Public documentation/capture inputs and assets retain S07's hash-bound reproduction evidence; public scenario/reference/inventory contracts were freshly checked after the merge. The browser wheel was rebuilt from dev's integrated Python code; later changes are JS only. Its preparation is a disposable testing build, not release preparation. No cache-bust token, public figure, reference output, Gallery/vendor asset or owner-maintained social preview was edited.

### Diagnostic failures and scope limits

The broad 29-case Session/Status diagnostic on the recolor fix passed 25 cases, failed two, and did not run two because the Session describe block is serial (`/tmp/issue602-dev-sync-status-session-final-report.json`). The first-use progress failure was corrected and the unchanged case passed on final bytes. The previously skipped divergent-draft and bare-legacy cases also pass when run separately from the failing direct-edit structural case.

The remaining extra case, `session-regenerate-intent.playwright.spec.js:1267` (“loaded current preview supports direct edits before the first Generate”), requires `artifactCheckpointBuilds` not to increase for a caption-group recolor. Actual count rises from 1 to 3. Source attribution is exact: dev introduced `commitSpecificRules` with `history.runUndoableCheckpoint`; both owners are byte-identical to dev. The caption fix adds previous-swatch evidence but does not change checkpoint allocation. The unmodified dev test encounters the earlier recolor rejection before reaching this assertion, so an exact base execution of the checkpoint failure is **not** claimed. The conflict between this existing structural test and dev's retained atomic History path is unresolved and remains a maintainer review item. The test, timeout, History owner and authority are unchanged. It is not relabeled PASS, waived by Product authority, or covered by the successful functional recolor tests.

The earlier 29-case attempt was intentionally interrupted after the concrete caption-rejection alert, with eight passes, one interrupted and twenty unrun; it was not a timeout conclusion. An intermediate Python-browser run was intentionally stopped for the pending production progress correction after seven passes (exit 2); it is not the final Python-browser result. A first browser launch could not start its local server; the next launch used a verified free port with server debug output. The initial failure supplied no test evidence and no test conditions were weakened. Two incorrect local diagnostic commands were corrected; their exits supplied no passing evidence.

## Independent review and preservation

- **Production:** the two initial Status owners and five integration files were reviewed separately from tests. New wiring retains the shared Worker, allowlisted helper operations, canonical request/admission, color evaluator, token/cancel checks and dev's atomic History. No status-only source work or alternate admission path was added.
- **Tests:** initial F1/F2 tests retain real result/draft/request/resource and unknown/Pending assertions. The added recolor cases observe public controls and the existing owning setup fixture; their Status assertion compares all preexisting Gallery differences, including the separate scale difference. All existing test assertions, mapped contracts and timeouts remain unchanged. New client tests observe actual initialization/cancel boundaries rather than a timer.
- **Docs:** historical S00–S07 and the initial Status-fix result are untouched. This new internal result records source changes, failure limits, reused evidence and review obligations. The public workflow owner remains `docs/REFERENCE/web-app.md`.
- **Generated:** wheel/log/report/trace files are disposable ignored or `/tmp` artifacts. Reference outputs, public examples, social preview, Gallery and vendor trees retain tracked bytes.

Authority SHA-256 remains `5451db7da8e4de970f05afd5d19b003b7de77d52eb84c05d71ff85a533d4b1fb`. The authority audit checks equality at HEAD/dev/prior reviewed base, 36 current and 45 historical Issue #602 signed receipt fields, 36 Issue #598 fields, and retained ancestors. The formal authority and signed decisions are not edited.

The initial 27 stash entries retain exact OIDs/messages. Shared checkout remains on `docs/issue-601-bug15-bug19-proposal-20260926` with all original dirty entries preserved. Another user's untracked `docs/internal/issue-619-proposal-20260927/` appeared during verification; it was observed and preserved, never added to these commits. The initial-state and final preservation receipts distinguish external additions from our changes.

## Gate, review and publication

Final runtime/test bytes pass the unchanged policy against `f5f86634459e0dcd46c1a452e9219fbba635d429` (`/tmp/issue602-dev-sync-final-working-budget.log`): **Gate PASS / maintainer Review REQUIRED**. The final commit-bound evaluation is retained in `/tmp/issue602-dev-sync-final-budget.log` and the handoff receipt. The documentation-only commit adds no runtime change.

Review reasons are architecture-bearing inventory changes, registered Session/compatibility owner changes, and production net additions above the architecture review threshold (537 > 400). The cumulative diff has 11 changed Web production files, 644 additions, 107 deletions and 751 lines of gross churn. Authority locations and privilege entries do not grow; static import cycles remain zero. Gate PASS does not waive the separate failing structural contract or grant publication approval.

Required cumulative PR jobs remain `web-change-budget`, `core-pr`, `recipes-standard`, `gallery`, `lint`, `web-contracts-pr`, and `web-pr-smoke`; trusted admission remains `Web base policy (trusted base)` and `PR / gate`. Local evidence does not substitute for hosted CI or human review. English PR files are `/tmp/issue602-pr-title.txt` and `/tmp/issue602-pr-body.md`; the final language checker uses those exact bytes once.

Inherited limits remain: Legend-override live-rerender retry binding; 223 px settings width at 195×422; 200% CSS-viewport and focus-with-resize keyboard evidence rather than physical-device evidence; sequential page/list scrolling on short screens; unproven initial S06 context destruction; two parent-dev Issue #564 metadata-free Session failures. The additional direct-edit checkpoint-contract conflict above remains unresolved. Complete functional/performance staging, Issue #601 runtime and release readiness are outside the completion claim.

Publication requires new explicit authorization for the intended same-named remote branch and PR, successful required hosted checks, and human review of artifact/Session responsibility, cumulative scope and the retained checkpoint-contract conflict. No permission or general Product reapproval is requested by this local handoff.

## Integration source hashes

- `gbdraw/web/js/app/app-setup.js`: `b301a70d0603c362515db260e850b92833bd152431d417916e0468101495c16d`
- `gbdraw/web/js/app/feature-editor/color-actions.js`: `9cb84242c206144d70c06294869cec31cfde98ddc6c1c11cdb396119d2a42082`
- `gbdraw/web/js/app/rule-matching.js`: `f06cae5d535ad4fa5f335c43761eeee5f4d125bfe748984db9b1a0aa3676b1cd`
- `gbdraw/web/js/app/run-analysis.js`: `f2023d1b918e83a1d03cf901428e4b10a16655746ee5569903a29e38483cc16b`
- `gbdraw/web/js/services/diagram-generation.js`: `22370e88543077552141229309e30076d04240b58b57e154406510001ff4dd76`
- `tests/web/color-captions.playwright.spec.js`: `eed53d96f1b3e2a9688e0ba5456e7f308be4ab58abe811bd5cf7acb3dd85e10c`
- `tests/web/generation-feedback.test.mjs`: `6951f96404ceeb8c3a6b37308afde763c81fe582b519a6a63d1952308ab8da7d`
- `tests/web/rule-matching.test.mjs`: `8618753508c740f3e651311ff6f2b9f1d031053a1552c01e6de998dd49ca227e`

Proposed session commit title: **Preserve Session edit Status through dev integration**.

Summary: Restore applied stroke and unresolved-selector Status, merge current dev without losing annotation/gap/History behavior, and repair complete-caption recolor and real first-use Generate feedback with bounded evidence.
