# Accepted Product Decision Pack: PD-OI-029 — ResetとHistory

## Identity and approval

- Author / Product Decision Owner: `satoshikawato`
- Decision date: 2026-09-26
- Concern: `diagram-generation.similarity-alignment.reset-and-history`
- Scenario revision: 3
- Status: ACCEPTED。選択、根拠、保持条件、廃止範囲、残余リスクは下記の完全な承認文で確定。
- Prepared from base SHA: `d457b7189b137185a8dec800819a312c30b969fa`
- Issue: https://github.com/satoshikawato/gbdraw/issues/598
- Authority update: `product/issue-598-decisions-20260926` / `eba01dd518ff6edd4effb4b00d519cbd78ff8623`
- Implementation: `fix/issue-598-alignment-direction-reset-20260926`
- 実装責任: S02 / S03 / S04

このPackは当該concernだけの承認記録であり、別concernの承認や実装を代替しない。正式authorityは [Option Integrity Product Contract](../../OPTION_INTEGRITY_PRODUCT_CONTRACT.md)。この文書はmachine decision storeではない。

## User journey and outcome

利用者は複数ゲノムの線形図を操作する。目標は、最新Alignの位置効果を解除し、必要ならそのAlignが変えた方向も復元する。 選ぶreferenceは正確なfeature identityであり、多数派方向の代表とは限らない。生物学的sourceと図の表示方向は別の情報である。

positions-onlyを既定にし、positions+alignment direction changesを選べる。combinedは実際の反転対象だけ絶対before方向へ戻す。両scopeは同じ最新baselineを使い、planとreceiptを消費する。

checkpointは操作前preview、final validation、成功artifact、Session fresh Load、regeneration、Undo/Redo、failure recoveryである。このPackの受入条件は総合計画書の **A07–A12、A14**。独立要件をすべて満たし、共通choice IDによって不足を代替しない。

## Authority search and non-waivable boundaries

| 調査対象 | 結果と扱い |
| --- | --- |
| Product map / contract | base OIPCの当該concernをscenario supersessionで更新する。mapやcheckerは変更しない |
| Durable BD registry | 調査baseに当該変更を認可するactive BDはなく、新しいBD storeは作らない |
| Source / tests | 現行Match checkboxと位置のみResetの観測根拠。Product authorityそのものとは扱わない |
| Domain / scientific integrity | source identity/strandを保持し、同一record transformと読みやすいtextを守る |
| Persisted compatibility | 公開済みSessionの証拠をS00で確定し、dev-only migrationを追加しない |
| Owner selection | 下記の9フィールドすべてが承認済み。未記入の根拠やriskを補ってserializeしない |

combined対象の後続manual direction editsが置き換わることを予告する。対象外編集とpending formを保つ。古い履歴欠落を推測しない。Save/fresh Load、style/reorder、complete Undo/Redo、failure rollback、追加Reset LOSATゼロを保つ。

architecture ratchet、canonical rendering、SVG sanitizer/admission、Worker資源管理、History readiness/finalization、scientific geometryのfailureはProduct選択で免除されない。BUG-17の新しいlinear crop/rotation代替は対象外。

## Implementation route and evidence

Route: **DURABLE_AUTHORITY_REQUIRED**。ownerの選択は完了しているが、当該authorityを含むdevをruntimeのbaseへ取り込むまではcandidate runtimeを認可しない。決定更新commitはContractだけを変更する。実装は [総合計画書](../01_MASTER_IMPLEMENTATION_PLAN.md) と対応する独立session promptに従う。

本Pack作成時点でruntimeは未実装。必要なevidenceは上記A条件に対応するsource-bound方向・位置・履歴・操作の実証であり、将来の成功を記録済み結果として扱わない。session evidenceは [保存方針](../evidence/README.md) に従う。受け入れ後も関心と責任を分離し、別concernのscopeをこのPackで広げない。

## Exact Product Decision Owner receipt

Supersedes PD-OI-029 revision 2. Receipt UTF-8 SHA-256: `8cd477428bb9ec0230f3c1f315db5ffffdda1a3478e69f42f3a2c25e05581549`.

```text
PRODUCT_DECISION
Concern: diagram-generation.similarity-alignment.reset-and-history
Scenario revision: 3
Choice: A / SELECTABLE_RESET_WITH_ALIGNMENT_DIRECTION_RESTORE
Rationale: Alignment can change both placement and record direction, so users must be able to choose whether Reset removes only positioning or also restores the directions actually changed by the latest Align. The restoration scope must be explicit without turning alignment plans into direction owners or replaying unrelated edits.
Must preserve: An explicit default Reset positions command that restores positions immediately before the latest successful Align, clears its active plan and keeps all current directions; an additional Reset positions and alignment direction changes command that uses the same position baseline and restores absolute before-Align direction only for records whose direction actually changed in that Align; unchanged direction for every record not reversed by that Align, including the reference when unchanged and all later manual edits on unaffected records; include a reference record in restoration only if an independently authorized alignment direction choice actually reversed it; visible target names, count, current and restored directions, and disclosure that later manual direction edits on restoration targets are replaced by the combined command; replacement of the reset receipt by each new successful Align, with no first-Align or source-orientation fallback; an orientation-independent plan and ordinary record-owned current orientation; source-bound validated restoration information captured from actual successful before/after states, retained across style regeneration and stable reorder, saved and freshly loaded with new Sessions, and cleared atomically with plan invalidation or successful Reset. Missing old restoration information leaves positions Reset available and direction restoration unavailable with an explicit reason, never guessed; an empty modern delta means no Align direction changes. Both Reset scopes consume the active plan and receipt, preserve unrelated settings and pending form edits, and use one canonical artifact transaction. Undo/Redo restores or reapplies complete artifacts including directions, plan and receipt; failed, canceled, stale, superseded, preview-readiness or History-finalization work creates no committed history and preserves the prior artifact and receipt. Original sources, biological identities, target-external settings, readable text, record-consistent feature/label/ribbon geometry, existing canonical request, rendering, Worker, sanitizer/admission and History owners, compatible comparison reuse and zero additional Reset LOSAT jobs remain intact. Restoration receipts never select rendering orientation.
May retire: The rule that Reset always retains directions changed by alignment and that their restoration is available only through ordinary Undo or manual Reverse. Do not retire the position-only choice or either existing recovery workflow.
Accepted residual risk: The combined scope intentionally replaces later manual direction edits on records actually reversed by the latest Align. Bound this to a visible target and before-direction preview, an explicit position-only alternative and one-operation Undo. A positions-only Reset consumes the same active plan and receipt; switching to combined afterward requires Undo of that Reset first. Sessions without historical direction evidence cannot restore it and must show that limit while retaining positions Reset. A compact Session/History receipt adds bounded state-maintenance cost; require source/plan binding, atomic lifecycle, old-information/no-op/re-Align/session/failure/geometry/job regression coverage, and keyboard/390 px acceptance. Guessing missing history, reversing unrelated records, partial restoration and new direction/render/History owners are not accepted.
Owner: satoshikawato
Decision date: 2026-09-26
```
