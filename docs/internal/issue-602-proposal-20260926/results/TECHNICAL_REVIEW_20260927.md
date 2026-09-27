# Issue #602 — cumulative technical review and bounded corrections

今回のローカル作業にユーザーの追加操作は不要です。AIによる技術レビューはmaintainer承認ではありません。将来のPRでは、下記の保存結果・色変更History・責任範囲を人間が確認する必要があります。署名済みProduct結果の再承認は求めません。

## Review conclusion and two findings

The cumulative production, tests, documentation and generated changes were reviewed separately against dev `d313b70b9f97c2c1d70f9ae885edbead80b62021`. Existing Status and bounded caption History corrections remain valid. Read-only observation of that dev's failed Tests run exposed additional inherited failures in owners touched by this branch. Two bounded corrections address four reproducible browser failures; no unrelated implementation is restarted.

1. **Settings-only Save was rejected before producing a download.** `exportSession` unconditionally added `annotationWarnings: []` to `runMetadata`. The existing inventory correctly rejects committed render metadata in a settings-only document, so three ordinary source-free Save workflows failed. The writer now emits `{}` for settings-only `runMetadata`. Full Result track geometry and annotation-warning persistence are unchanged; the inventory/reader/schema is not relaxed. Base authority: `docs/SESSION_COMPATIBILITY.md`, “Session 42: settings before the first source.”
2. **A captured mixed-comparison Generate returned stale after future comparison controls changed.** `setLinearComparisonGlobalAction` invalidates `files.linearCanonicalComparisons`, a derived reuse artifact owned by the committed request. The color-preparation snapshot incorrectly treated that cache identity as an input. Only this derived cache is excluded from that snapshot. Biological and color sources, catalogs, Result identity/names/selection, rules, legend state, mode, generation token and cancellation remain guarded. The unchanged mixed-render test now finishes its captured request and reuses the raw cache; a new asynchronous Node case verifies that derived cache invalidation stays current while biological source replacement becomes stale. Base comparison inventory prohibits binding this artifact outside the committed request; its captured comparison resolution remains authoritative. `PD-OI-016` / `OIC-013` failure isolation and `PD-OI-037` / `OIC-024` truthful Status remain jointly required.

Developer preflight: `IMPLEMENT_EXISTING_AUTHORITY`. No supported outcome is retired, no materially different Product option remains unresolved, and no Product response or authority edits are requested. Existing request, artifact, Session, History, legend, color preparation and Worker owners remain. The two fixes add no semantic owner, canonical entry, persisted adapter or compatibility branch; ordinary OE/PE/CB evidence is non-increasing. Full before/after exception sets are not needed.

## Required cumulative review points

- **A — Result versus draft after Save/fresh Load.** Saved editor evidence restores only proven applied Block Stroke Width, including zero. Automatic width stays request-owned; missing proof makes only width Unknown. Unresolved Circular selection stays Unknown without erasing independent Pending differences. A live field commit clears only its proven Unknown. The draft is never adopted wholesale as the applied Result.
- **B — One reversible caption operation.** Pure recolor uses bounded intent History, while legend inventory additions/removals keep checkpoint History. The canonical transaction performs one guarded commit after legend preparation. Vue inherited refs, saved swatches and History-only legend ownership survive Undo/Redo and a subsequent independent rule. Collision and stale candidate checks protect state and History. No legend ownership metadata is exported as Session fields.
- **C — Existing owner/path boundaries.** Shared generation-intent projection owns Status. Status adds no discovery, biological-source read/hash, Worker work or SVG/checkpoint clone. Request/render/import admission, History transaction, legend geometry and the shared Worker retain their owners. PR #618 documentation/figures and PR #620's unchanged accepted policy remain integrated; Issue #601 runtime remains outside scope.

Tests were audited separately: existing numeric/geometry oracles, structural checkpoint metrics, cold Worker assertions, collision/stale behavior and complete-file coverage remain intact. No checker, timeout, comparison tolerance or test-owned expectation was weakened. The new test exercises an async candidate rather than merely mirroring the filter expression.

Public docs/capture recipes were reviewed separately from internal records. The tracked Lambda–DE3 figure was viewed at readable scale; the unchanged image, original input hashes and twelve-checkpoint capture evidence agree. Generated/vendor/Gallery/reference/social-preview trees retain their protected bytes. No public capture or reference regeneration was performed.

## Exact local boundary

Dedicated worktree: `/tmp/gbdraw-issue602-s00-results-20260926`.

- Continued branch: `fix/issue-602-linear-live-edit-20260926`.
- Upstream: `origin/fix/issue-602-linear-live-edit-20260926`.
- Reviewed starting HEAD: `0173dc2fb337dd0ec4ee8f1ea0fbc2514dc945f5`; prior runtime HEAD: `9b3f6f8b0b96cd2f6a17a92a2d4b744797418279`.
- Comparison base / fetched dev: `d313b70b9f97c2c1d70f9ae885edbead80b62021`; remote work branch: `9a14db3af8fafd6f8374cfe2a8197a8927cf43f7`.
- Final HEAD and clean-tree/remote observations are bound by `/tmp/issue602-review-final-evidence-20260927.json`, avoiding a self-referential commit ID.

Fetch and read-only remote/PR/CI checks found no new dev to merge and no existing PR. No branch recreation or history replacement occurred. The shared checkout remains outside this task and is not switched, edited or staged.

## Verification and evidence reuse

`/tmp/issue602-review-20260927-audit.json` verifies the handed-off receipt before edits: all 67 cumulative-source SHA-256 values, ten artifact hashes, ten accepted command receipts and their actual execution sources (including dirty files), logs/reports and preserved partial runs. The original 36-case union is 32 cases outside the diagnostic four-case file plus its final four-case success. It is not a single 36-pass command. Python browser coverage remains 38 cases across commands, including corrected comparison shards. The prior 989-pass root run and overlapping 192 checks are preserved as historical evidence.

Fresh verification after these two runtime corrections is recorded with command/start/exit/source inventories under `/tmp/issue602-review-<name>-20260927.*`:

| Command | Result / acceptance boundary |
| --- | --- |
| `node` | 990 root Web Node tests passed, including the added asynchronous guard regression; zero failures/skips/cancellations. |
| `architecture` | 192 architecture/Product ratchet checks passed; overlap with root contracts remains explicit. |
| `corrected-final` | All four originally failing source-free/mixed-render cases passed, no retries or exclusions. |
| `session` | The complete unfiltered Session-regeneration file passed all four cases, including direct edits, divergent draft and bare legacy import. |
| `colors` | All six caption cases passed, including collisions, both legend orientations, 1440/390 px recolor, Undo/Redo, ownership after an independent edit, Save/fresh Load and export. |
| `smoke` | All thirteen official PR-smoke cases passed. |
| `annotation` | All eighteen annotation/style/alignment/Reset cases passed, retaining full-Result warning persistence and PR #618 assertions. |
| Lint / whitespace | Ruff and cumulative whitespace checks passed. |

Each accepted fresh receipt asserts unchanged execution sources; its actual dirty-file hashes match the committed final runtime/test bytes. The final receipt lists every accepted command and report digest. No single all-repository success or full staging claim is made. The old 32 Session cases and Python-browser executions are retained as historical audited evidence, not recounted as fresh tests of the corrected runtime. Unchanged public documentation inputs retain the 56-pass documentation-contract receipt.

Native reuse stays restricted to 667 unchanged native production/test/fixture/example/dependency blobs versus actual supervised source `cb2ad027d3f60fde2331a23ad8f9f6274fe863ba`. Raw working bytes and retained log digests were checked. Its 6,642-pass / 17-skip / 11-deselected evidence, including 189 recipes and 103 Gallery cases, does not cover changed Web owners. Public figures and reproduction inputs remain unchanged; their valid evidence is reused.

### Diagnostics retained

- Before the new fixes: five focused browser cases produced one pass and four failures. The successful inherited Feature fill case is preserved; all failures and traces remain in `boundaries`.
- The first corrected selection used a broad Session title filter and completed ten browser cases successfully, but the outer wrapper ended with code 143 before its exit receipt. Its report and matching start/current source hashes are retained in `/tmp/issue602-review-corrected-20260927.termination.json`. Termination cause is unproven. It is diagnostic evidence, not an accepted exit-0 invocation; `corrected-final` supplies a clean four-case receipt.
- The historical approximately 0.000824 SVG-unit fresh-Load difference remains unproven. This review's complete four-case Session file passes without changing geometry, comparison tolerances or timeouts; earlier failures remain intact.
- One read-only draft/log inspection was not executed when automatic approval review timed out; the permitted unchanged single retry completed. This was not a safety rejection or publication attempt.
- Original CLI-shadow/occupied-port Python-browser attempts, old checkpoint/swatch failures and S06 context destruction remain historical. An invalid classifier invocation without required environment is a diagnostic only; the final exact-head invocation uses the required inputs.

## Remote CI and remaining limits

Read-only dev Tests observation: [run 36294259022](https://github.com/satoshikawato/gbdraw/actions/runs/36294259022) failed; Gallery and CodeQL passed. The failed run contains nine browser failures and a canceled shard, so latest dev staging is not green. Three source-free cases and mixed-render failure are now corrected locally. Inherited Feature fill passes locally. Two progress failures are covered by the existing Issue #602 progress correction. Two Issue #564 metadata-free Session cases remain outside the claim. These local results do not turn the remote dev run into a success.

Inherited limits are unchanged: Legend-override rerender retry binding; 223 px settings at 195×422; CSS-viewport/focused-resize device evidence; sequential short-screen scrolling; unexplained S06 context destruction; Issue #564's two cases; the intermittent subpixel diagnostic. Full staging, hosted candidate CI, Issue #601 runtime and release readiness are not claimed.

## Gate, human review and preservation

Final exact-head policy/CI receipts: `/tmp/issue602-review-final-budget-20260927.log`, `/tmp/issue602-review-final-ci-impact-20260927.json`. Gate PASS and maintainer Review REQUIRED remain separate. Cumulative production is fourteen files, 729 additions / 163 deletions, gross 892 / net 566. Architecture and registered Session responsibility plus the file/line thresholds route human review; AI technical review does not satisfy that requirement.

Future maintainer confirmation is concrete: inspect the final diffs and source-bound Session/Status and caption traces for A/B above; inspect the source-free metadata condition and derived-cache-only filter alongside the biological-source rejection and token/cancellation guards for C. The reviewer decides whether those existing owner boundaries and cumulative scope are acceptable. No general Product reapproval is needed.

Formal Contract remains revision 22 / SHA-256 `5451db7da8e4de970f05afd5d19b003b7de77d52eb84c05d71ff85a533d4b1fb`; current/historical Issue #602 and #598 records, retained ancestors, PD-OI-027/029/031/034 revisions 5/3/5/5, PD-OI-035 revision 3, PD-OI-038, PD-OI-039 revision 2 / `A / EXCLUSIVE_DIRECTIONS_WITHOUT_MATCH` and PD-OI-046/047 are preserved. Historical Match is not revived. S00–S07 and all 27 original stash OIDs/messages remain unchanged. Previous result records are not overwritten. Protected policy, workflows, renderer/geometry/History owners and generated trees remain unchanged relative to accepted dev where previously required.

English PR wording stays at `/tmp/issue602-pr-title.txt` and `/tmp/issue602-pr-body.md`; the former bytes are archived as `/tmp/issue602-review-pr-title-before-20260927.txt` and `/tmp/issue602-review-pr-body-before-20260927.md`, matching their historical receipt hashes. The materially updated wording is checked once with its exact bytes. No push, PR mutation, remote merge, deploy, tag, stash action, reset, amend or rebase is performed.

Proposed commit title: **Fix settings-only Session saves and comparison cache invalidation**.

Summary: Review the cumulative Issue #602 changes, fix two inherited Session/generation defects within existing owners, and record source-bound verification without changing Product authority.
