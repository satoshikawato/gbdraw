# Issue #602 — verified caption History and latest dev synchronization

今このローカル作業を終えるために、ユーザーの追加操作は不要です。将来PRをレビューする際は、次の3点が確認対象です。署名済みProduct決定の再承認は求めません。

1. **保存・Loadの結果**：保存したResultと未適用の設定を分けて保持する。適用済み線幅はApplied、未解決Circular selectionはUnknown、独立した設定差はPendingを保持する。
2. **直接編集とUndo/Redo**：caption全体の色変更を1操作として戻せる。凡例の色も復元し、Undo後に別のルールを追加しても元の凡例が消えない。凡例の追加・削除はdevのcheckpoint処理を保持する。
3. **変更の責任範囲**：既存のrequest、artifact、History、Session、Worker所有者を使う。Historyだけの凡例所有情報を保存Sessionに追加しない。累積14 production files / net 563 linesの範囲を確認する。

The previous `DEV_SYNC_STATUS_20260927.md` is historical. Its direct-edit checkpoint conflict is resolved by this follow-up; S00–S07 and every earlier result remain intact. This record grants no publication authority.

## Local history and exact comparison

Dedicated worktree: `/tmp/gbdraw-issue602-s00-results-20260926`. Branch: `fix/issue-602-linear-live-edit-20260926`; upstream: the same-named origin branch. The remote work branch remains `9a14db3af8fafd6f8374cfe2a8197a8927cf43f7`.

| Local commit | Purpose |
| --- | --- |
| `7955c02f3eb9147566eb93ae44aa92b0b38ccbed` | Save the handed-off applied-width and unresolved-selector Status corrections. |
| `11c81b1cea80caa8b27012cdaef19ec7e16cf830` | Normal merge of dev `f5f86634459e0dcd46c1a452e9219fbba635d429`. |
| `f6d44139edcdc47fd3c84a967a873a94bfe8d160` | Preserve complete default-caption recolor and real cold Generate progress. |
| `6daee48a2446c3007e6c52a504bd5e55b7a5dc1d` | Record the then-current integration and unresolved structural case. |
| `e63493ec3fb428c511dda1f160de7079a3093a2a` | Keep caption recolor History bounded and reversible. |
| `952a6293e6165cdad1bed98690c52252204b299b` | Normal merge of dev `88028fd242d263f0fe86aaf9da57b8dc9eb082f6` / PR #618. |
| `9b3f6f8b0b96cd2f6a17a92a2d4b744797418279` | Normal merge of dev `d313b70b9f97c2c1d70f9ae885edbead80b62021` / PR #620. |

Final runtime HEAD is `9b3f6f8b0b96cd2f6a17a92a2d4b744797418279`; comparison base and latest observed dev are `d313b70b9f97c2c1d70f9ae885edbead80b62021`. The documentation-only commit containing this record is identified by `/tmp/issue602-followup-final-evidence.json`, avoiding a self-referential commit ID.

The first merge retained both complete Status/Lock and dev gap test blocks. The PR #618 merge conflicted only in `docs/REFERENCE/web-app.md` and `docs/capture/README.md`: the updated alignment/direction/Reset explanation remains, including the compact Editor/reopening/scroll continuation, and both independent capture workflows remain. Its updated Gallery/reference-only Custom alignment assertions are included in fresh browser verification.

PR #620 changes only `tools/web-change-policy.json`. The accepted upstream bytes are imported unchanged, permitting the independently approved Session import Worker constructor and codec import. No candidate authority edit or Issue #597 runtime implementation is added. The formal Product Contract remains revision 22 with its original SHA-256.

## Correction and owner/path evidence

The unchanged direct-edit case previously observed `artifactCheckpointBuilds` rising from 1 to 3 for a caption recolor. `commitSpecificRules` now consumes the existing canonical legend diff: color-only changes use existing bounded intent History; legend additions/removals retain existing checkpoint History. The legend owner finishes preparation and measurement before the guarded synchronous rules/geometry/Result commit inside the selected transaction. Stale candidates and collision failures still leave state and History unchanged.

Actual browser investigation exposed two additional requirements for bounded History. Vue refs expose `value` through their prototype, but the existing intent read/write helpers required an own property and recorded an empty legend. Those same helpers now recognize inherited `value`. On restoration the legend owner uses the captured swatch instead of overwriting it with the original palette. History also captures small per-target caption ownership records from the same legend owner and restores them, so Undo cannot leave a default entry falsely owned by the specific-color table.

Ownership metadata lives only in History intent. It is not attached to reactive legend entries, exported editor state, or a Session schema/reader. No full artifact clone, catalog build, Result serialization, biological-source work or Status-only side effect is added by this History metadata capture. `history.js`, watchers, legend geometry/layout, renderer and Worker implementations are preserved. Ordinary OE/PE/CB evidence is non-increasing: existing canonical legend, History snapshot and transaction owners remain; no competing owner, persisted adapter or execution entry is introduced.

Product preflight remains `IMPLEMENT_EXISTING_AUTHORITY`. The existing atomic-edit, Undo/Redo, separate draft/Result, and truthful Status continuations are jointly preserved. No materially different product option is selected.

## Verification and source binding

Every follow-up run has a command/start/exit receipt, per-file runtime/test SHA-256 inventory and immutable log under `/tmp/issue602-followup-<name>.*`. Successful runs assert that those source files were unchanged during execution. Final receipt binds the local HEAD, base, policy results, accepted commands and review files.

| Evidence | Result and boundary |
| --- | --- |
| `node-corrected` | 989 root Web Node tests passed: 849 CI fast contracts, 139 architecture contracts and one Gallery publication contract. |
| `architecture-final` | 192 architecture/Product ratchet tests passed after PR #620 policy input changed; this overlaps the root run. |
| `direct-corrected` | 4 browser cases passed: unchanged direct-edit contract, collision rejection, full default-caption recolor at 1440/390 px, including legend Undo/Redo and independent subsequent edits. |
| `canonical-color` | 3 browser cases passed: multicolor captions at 1440/390 px and both Linear legend orientations; native comparison, Session/export and topology changes included. |
| `session-final` | All 4 cases in the unfiltered Session-regeneration file passed after PR #618 merge, with no retry or exclusion, including the former direct-edit failure, divergent draft and bare legacy case. |
| `session` + `session-final` | 36 distinct Session/Status/History/palette/mode cases are covered: 32 other fresh successful cases have unchanged code, inputs and acceptance conditions; the complete 4-case Session file supplies final success. This is not represented as a single 36-pass invocation. |
| `annotation-alignment` | All 18 annotation/style/alignment/Reset cases passed on merged PR #618 assertions, including both modes at 1440/390 px and reference-only Custom alignment. |
| `smoke` | All 13 official PR-smoke cases passed, no retry or skip. |
| Python browser | 36 browser cases passed; both comparison-contract shards then passed with Node Playwright PATH and a free dedicated port. These are 38 distinct pytest cases across commands, not one all-pass invocation. |
| `docs-final` | 56 public-reference/scenario/capture/figure-inventory/link contracts passed after documentation conflict resolution. |
| Audits | Ruff, whitespace, authority/receipts/ancestry, source binding and preservation pass; exact final policy/CI classification is attached to the final receipt. |

PR #618 changes no Web runtime, Node test, Gallery session, example inventory, dependency or comparison fixture bytes. Its changed alignment browser spec is not executed by earlier accepted commands and is freshly verified in `annotation-alignment`. PR #620 changes only the accepted policy input, which is freshly covered by `architecture-final` and the exact-head policy check. These boundaries permit the earlier accepted source-bound successes to apply to the final runtime; no old unmodified-source assertion is used for a changed owner.

Native reuse remains restricted to 667 unchanged native production/test/fixture/example/dependency blobs, compared in batches with actual supervised dev source `cb2ad027d3f60fde2331a23ad8f9f6274fe863ba`. Its existing 6,642-pass / 17-skip / 11-deselected log and receipts include 189 recipe and 103 Gallery cases. This is reused native evidence; modified Web/History and newly merged documentation inputs do not inherit that entire run. Public reproduction evidence is retained for unchanged figures; PR #618's updated figures/capture flows remain exact upstream assets. No public figure or reference output was regenerated here.

### Diagnostic outcomes retained

The before-fix unchanged direct-edit case failed at the checkpoint assertion. The first bounded-History trial passed that case but failed Redo legend-color persistence; the added Undo-swatch assertion then exposed the missing Vue ref capture. Final tests observe restored swatches, valid ownership after Undo and an independent subsequent rule.

The first expanded 36-case browser diagnostic passed 34, failed divergent-draft SVG equivalence and did not run its following serial legacy case. The observed difference was a legend/canvas width change of approximately 0.000824 SVG units on fresh Load. No comparator, geometry owner or timeout was changed. The same two cases passed in isolation, and the complete final four-case Session file also passed. The earlier failure and unproven intermittent measurement cause remain recorded, not erased or counted as success.

Two Python-browser attempts passed 36 but could not execute the comparison shards: first the Python Playwright executable shadowed Node's CLI; then the default server port was occupied. Correct PATH and a free dedicated port produced both shard successes. One invalid docs test selector collected zero cases and was replaced by the existing inventory test names. One intermediate verification was intentionally interrupted to correct a ref-condition syntax error; its exits are not successful evidence. Test-owned timeouts, authorities, or acceptance checks were not weakened.

The PR #620 merge's first automatic approval review timed out before execution. Read-only HEAD/status confirmed no merge started; the tool-authorized single retry completed the ordinary local merge. This was not a safety-policy rejection or publication action.

## Independent review and preservation

Production was reviewed separately from tests: four follow-up production files change existing owners and retain canonical admission, rule evaluation, guarded atomic commit, shared Worker and both existing History mechanisms. Tests separately retain the unchanged direct-edit structural oracle and add behavioral coverage for inherited refs, immutable History metadata, restored legend swatches and subsequent independent edits. Documentation separately preserves both public workflows and new dev explanations. Generated assets separately remain disposable logs/reports/traces/wheel; Gallery/vendor/reference/social-preview trees match latest dev.

Formal authority SHA-256 stays `5451db7da8e4de970f05afd5d19b003b7de77d52eb84c05d71ff85a533d4b1fb`. Current/historical Issue #602 receipt fields, all Issue #598 fields, PD-OI-027/029/031/034 revisions 5/3/5/5, PD-OI-035 revision 3, PD-OI-038, PD-OI-039 revision 2 / `A / EXCLUSIVE_DIRECTIONS_WITHOUT_MATCH`, PD-OI-046/047 and prior ancestors remain intact. Issue #601 runtime is outside scope.

All 27 original stash OIDs/messages and S00–S07 bytes remain unchanged. Another session independently switched the shared checkout from its original docs branch to `fix/issue-619-circular-track-measure-inputs`, updated dev and committed its Issue #619 plan. Read-only reflog/status documents that external change. All seven original untracked directory entries and the added Issue #619 directory remain. This task never switches or edits the shared checkout or incorporates those changes into its commits.

## Gate and remaining review boundary

Exact-head Gate is PASS; Review remains REQUIRED for architecture inventory/owner changes, registered Session responsibility, 14 production files > 12 and net 563 additions > 400. Production additions/deletions are 724/161 (gross 885). The former checkpoint-contract failure is resolved and is not a remaining maintainer implementation task. Privilege changes relative to latest accepted dev are zero; static import cycles remain zero.

Required PR jobs remain `web-change-budget`, `core-pr`, `recipes-standard`, `gallery`, `lint`, `web-contracts-pr`, `web-pr-smoke`; trusted admission remains `Web base policy (trusted base)` and `PR / gate`. Human review concerns the concrete three points at the start of this record and cumulative scope. It is not general reapproval of signed Product outcomes. Hosted CI, any later publication authorization and complete integrated staging remain separate future boundaries.

Inherited limits remain those in the earlier result: Legend-override rerender retry binding, 223 px settings at 195×422, CSS-viewport/focused-resize device evidence, sequential short-screen scrolling, unexplained initial S06 context destruction, and two parent-dev Issue #564 metadata-free Session cases. The intermittent subpixel measurement diagnostic above also remains disclosed. Complete staging, Issue #601 runtime and release readiness are not claimed.

English PR title/body: `/tmp/issue602-pr-title.txt`, `/tmp/issue602-pr-body.md`. The materially revised final wording is checked once against its exact bytes; prior wording-check evidence is historical. Existing PR lookup is read-only. No Push, PR creation/edit, remote merge, deploy, tag, branch recreation, reset, amend, rebase or stash operation is performed.

Proposed session commit title: **Preserve Session edit Status and bounded color History**.

Summary: Restore applied and unresolved-selector Status, merge current dev, and keep complete-caption edits atomic and reversible without adding full History checkpoints for recolor.
