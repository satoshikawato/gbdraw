# Decision Packet: Match affordance in compact direction review

Status: RESOLVED; documentary preflight and exact approved receipt. The normative owner remains the base-branch Contract.

- Concern: `web.similarity-alignment.review-presentation`
- Proposed scenario revision: 2 (supersedes PD-OI-039 revision 1 only as explicitly selected)
- Base: `302dfa1136c95ab50ddef606ad835ec4087b95b2` (`origin/dev`)
- Candidate: `fix/issue-598-alignment-direction-reset-20260926`, S02 after `d5edc4bb00b0d793d94361202c31a270f78a4aa0`
- Prepared by: Codex for Product Decision Owner `satoshikawato`
- Related issues: #598, #602

## Original trigger and authority search (resolved)

Replacing the reference-relative Match checkbox with Keep/right/left/Custom
would remove an independently preserved affordance in PD-OI-039. The four
Issue #598 receipts match the actual dev authority in all nine fields at
revisions 5/3/5/5. They authorize exclusive direction intent and retirement of
the old checkbox/flag, but do not explicitly supersede PD-OI-039's complete
receipt. PD-OI-035 revision 3 and PD-OI-039 revision 1 are jointly required.

| Source | Finding |
| --- | --- |
| Static Contract, PD-OI-027/031/034 | Exclusive Keep/right/left/Custom; old Match checkbox/flag may retire |
| Static Contract, PD-OI-035/039 | Independent exact identity, Select/Skip, canvas/focus and compact review requirements; PD-OI-039 expressly retains “1つのMatch reference direction” |
| Product Impact map and BD registry | No mapped/durable decision resolving this intersection; active BD array is empty |
| Domain/integrity | Source strand, exact identity and unknown exclusions remain mandatory; neither establishes whether Match remains discoverable |
| Released compatibility | No persisted review policy is needed or permitted; S00/S01 main/release observations remain applicable |
| Eligible PR-local decision | No PR-local decision; retirement cannot use that route |
| Current code/tests | One Match checkbox and local preview, evidence only; not a supersession rule |

Authority resolution: CONFLICT for retirement of the Match affordance.
Classification: PRODUCT_DECISION_REQUIRED. More rendering evidence cannot
choose whether that affordance remains. No date/order precedence is inferred.
NOT_ALLOWED does not classify the complete task: independent canonical
transaction/Reset work remains authorized. IMPLEMENT_EXISTING_AUTHORITY does
not resolve the Match retirement. EVIDENCE_REQUIRED alone cannot resolve it.

## Journey and independent realization requirements

A Web user opens an exact-reference alignment review, selects an ambiguous or
non-rendered candidate or Skip, inspects direction arrows and the canvas, and
Applies or retries. Checkpoints: local preview, Python final validation,
commit, recovery, keyboard/focus, fresh Session and History.

| Jointly necessary contribution | A | B |
| --- | --- | --- |
| Keep/right/left/Custom is one exclusive intent; reference participates | Preserve | Preserve |
| One discoverable Match action | Retire expressly | Preserve as a one-shot action that selects Custom outcomes; never a parallel flag/policy |
| Python eligibility, exact reference, keyboard Select/Skip, non-rendered candidates | Preserve | Preserve |
| 390×844/740 canvas at full available width and at least 200 px height, scrollable candidates, reachable Apply/Cancel | Preserve | Preserve |
| Desktop non-modal canvas, pan/zoom, wide drag, focus return, overlays excluded from artifacts | Preserve | Preserve |
| Narrow Editor owner closes it, retains tab, disables opening with a reason, explicit reopen after review | Preserve | Preserve |
| Local zero-Worker edits, one final batch, refreshed material preview/another Apply, atomic recovery | Preserve | Preserve |
| Reset receipt, both Reset scopes, Session/History and ordinary Reverse | Preserve | Preserve |

This is AND-of-OR coverage: a direction control cannot substitute for identity,
canvas, keyboard, Editor lifecycle, failure recovery or persistence.

## Choice A — EXCLUSIVE_DIRECTIONS_WITHOUT_MATCH

Keep the already accepted exclusive direction controls and expressly retire the
one Match affordance from PD-OI-039. All its other independent requirements
remain. Users choose right/left or Custom to obtain their intended output.
The old checkbox, flag and reference-relative projection are removed.

Preserved: all four accepted Issue #598 decisions, exact identity and all
non-Match PD-OI-035/039 effects in the matrix. Added: no additional behavior
beyond Issue #598. Retired: only PD-OI-039's discoverable Match action. Risk:
users accustomed to Match need to choose arrows or Custom explicitly.
Route: DURABLE_AUTHORITY_REQUIRED; explicit receipt, authority-only PD-OI-039
revision-2 supersession merged to dev, then dependent runtime integration.

## Choice B — EXCLUSIVE_DIRECTIONS_WITH_MATCH_ACTION

Retain one named Match action. Activating it makes a one-time local selection
of Custom directions for currently selected known-strand targets using the
reference's currently displayed arrow, keeps the reference, and leaves unknown,
skipped/missing/unusable records unchanged. The resulting Custom values are
editable. Subsequent candidate/reference edits do not maintain a hidden Match
policy: users see/reselect explicit outcomes. Keep/right/left/Custom remain one
intent; there is no independent boolean flag or persisted policy.

Preserved: all independent requirements in the matrix, including discoverable
Match. Added: a one-shot local continuation into Custom. Retired: checkbox and
persistent reference-following intent, not the affordance. Risk: users may expect
Match to follow subsequent candidate edits; disclose its one-shot scope.
Route: DURABLE_AUTHORITY_REQUIRED; explicitly define this continuation and its
limited supersession in authority-only PD-OI-039 revision 2, merged to dev
before dependent runtime. This packet does not infer acceptance by analogy.

## Comparison and constraints

Both choices preserve automatic Keep, arrows, candidate evidence, source bytes,
validation/error/retry, prior artifact on failure, absolute record-owned
transforms, receipt/positions, complete Undo/Redo, fresh Session/regeneration,
Export, compatible comparisons, one renderer/admission/History owner and
resource lifecycle. Review choices/overlays never enter exported artifacts.
Neither adds Worker calls to local editing or LOSAT calls to Reset. Both require
browser/390 px/keyboard acceptance; those checks cannot waive missing authority.
A loses Match discoverability; B adds a disclosed explicit Custom-selection
shortcut and requires tests after candidate edits. Neither changes scientific
strand interpretation or comparison provenance.

Keeping the existing checkbox as an additional independent flag, silently
retiring canvas/Editor requirements, guessing unknown directions, changing
source strand or adding a second render/History owner is NOT_ALLOWED.

Engineering recommendation: A matches the approved #598 operation model with
fewer actions. This recommendation is not Product authority.

## Product Decision Owner response

Supply all nine fields. A letter alone does not resolve preservation, retirement,
rationale or risk. This does not reopen the four already accepted decisions.

```text
PRODUCT_DECISION
Concern: web.similarity-alignment.review-presentation
Scenario revision: 2
Choice: <A / EXCLUSIVE_DIRECTIONS_WITHOUT_MATCH or B / EXCLUSIVE_DIRECTIONS_WITH_MATCH_ACTION>
Rationale: <product-level reason>
Must preserve: <all independent PD-OI-035/039 requirements and all four accepted Issue #598 outcomes; for B also define the one-shot Match continuation>
May retire: <explicit Match affordance scope for A, or old checkbox/flag scope for B; no other retirement>
Accepted residual risk: <bounded risk or none>
Owner: satoshikawato
Decision date: <YYYY-MM-DD>
```

After the explicit response, serialize only that outcome for human review in an
authority-only change. Candidate authority never authorizes the same runtime.
The owner selected A and approved the full displayed nine-field draft. Authority-only PR [#609](https://github.com/satoshikawato/gbdraw/pull/609) merged into dev as `2edc00aebc74e01003da643dfc957b513d5dcfe5`; the implementation branch includes that dev merge through `0453380087cbae640cf0c3b9792a1b7d1d023213`. PD-OI-039 scenario revision 2 authorizes the limited Match retirement, retaining the four #598 receipts and all independent PD-OI-035/039 requirements.

## Approved machine representation

The owner explicitly approved the complete draft with UTF-8 SHA-256 `06a9d2fe9b1d1406f6f8e04c23a9ca031133b9fae24b0683403ec2c6cae55270`. This exact JSON matches the merged Contract and is reproduced for review, not as a second authority store.

```json
{
  "concern": "web.similarity-alignment.review-presentation",
  "scenarioRevision": 2,
  "choice": "A / EXCLUSIVE_DIRECTIONS_WITHOUT_MATCH",
  "rationale": "狭いPreviewでもalignment候補をcanvasで確認できるよう、reviewを図の下段に固定し、候補比較へ操作を集中させる。表示方向はIssue #598のKeep/right/left/Customへ統一し、reference相対のMatch操作による結果との混同を避ける。",
  "mustPreserve": "PD-OI-031/034と現行transform/plan/reset/historyのすべての結果。resolvedの通常自動Apply、ambiguousと明示reviewのlocal draft、独立Select/Skip、候補根拠とreference identity、Issue #598で承認済みの排他的Keep/right/left/Custom、referenceを含むselected known-strand anchorsの方向選択と各recordのbefore/after矢印、unknown/skipped/missing/unusableの理由付き不変、canvas操作、local編集でWorkerを呼ばないこと、Applyの共有Python batch validationとatomic Result/History。失敗時draft/error/retry、Cancel/stale/superseded時の以前のResult/orientation/History、Session/regeneration/Export、focus復帰を維持する。390×844/740では利用可能幅全体かつ高さ200px以上のcanvasを確保し、候補listをscroll、Apply/Cancelを到達可能にする。狭いreview開始時はEditorをownerで閉じ、tabを保持し、review中は理由付きでopenをdisable、終了後は明示reopen可能。wideのdragと非モーダルcanvasを維持する。",
  "mayRetire": "旧Match reference direction checkbox・flag・操作affordance。狭いPreviewでreviewを自由にdragする操作、およびreview中にEditorを同時openする継続。これ以外のPD-OI-035/039の独立要求とIssue #598の承認済み4決定は退役しない。",
  "acceptedResidualRisk": "狭いreviewではlist scrollが増え、自由に位置を動かせなくなる。開始時Editorは閉じるがtabは保持し、終了後再openできる。位置変更で候補draftやResultを変えないことをbrowserで確認する。旧Match利用者はright/leftまたはCustomで表示方向を明示的に選ぶ必要がある。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```
