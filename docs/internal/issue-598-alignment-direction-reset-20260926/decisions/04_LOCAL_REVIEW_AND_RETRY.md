# Accepted Product Decision Pack: PD-OI-034 — draft・検証・再試行

## Identity and approval

- Author / Product Decision Owner: `satoshikawato`
- Decision date: 2026-09-26
- Concern: `web.similarity-alignment.choice-and-retry`
- Scenario revision: 5
- Status: ACCEPTED。選択、根拠、保持条件、廃止範囲、残余リスクは下記の完全な承認文で確定。
- Prepared from base SHA: `d457b7189b137185a8dec800819a312c30b969fa`
- Issue: https://github.com/satoshikawato/gbdraw/issues/598
- Authority update: `product/issue-598-decisions-20260926` / `eba01dd518ff6edd4effb4b00d519cbd78ff8623`
- Implementation: `fix/issue-598-alignment-direction-reset-20260926`
- 実装責任: S01 / S02 / S03

このPackは当該concernだけの承認記録であり、別concernの承認や実装を代替しない。正式authorityは [Option Integrity Product Contract](../../OPTION_INTEGRITY_PRODUCT_CONTRACT.md)。この文書はmachine decision storeではない。

## User journey and outcome

利用者は複数ゲノムの線形図を操作する。目標は、candidate/Skip/方向intentをreview内で変更し、validationと描画を経て適用または再試行する。 選ぶreferenceは正確なfeature identityであり、多数派方向の代表とは限らない。生物学的sourceと図の表示方向は別の情報である。

draftは一つのtagged union、Custom row値はCustom時のみ。編集はWorkerなし。Applyごとのfinal batch validationは一度。最終方向/中心補正がpreviewと違えばrefreshして別Applyを待つ。

checkpointは操作前preview、final validation、成功artifact、Session fresh Load、regeneration、Undo/Redo、failure recoveryである。このPackの受入条件は総合計画書の **A06、A10、A11、A13、A14**。独立要件をすべて満たし、共通choice IDによって不足を代替しない。

## Authority search and non-waivable boundaries

| 調査対象 | 結果と扱い |
| --- | --- |
| Product map / contract | base OIPCの当該concernをscenario supersessionで更新する。mapやcheckerは変更しない |
| Durable BD registry | 調査baseに当該変更を認可するactive BDはなく、新しいBD storeは作らない |
| Source / tests | 現行Match checkboxと位置のみResetの観測根拠。Product authorityそのものとは扱わない |
| Domain / scientific integrity | source identity/strandを保持し、同一record transformと読みやすいtextを守る |
| Persisted compatibility | 公開済みSessionの証拠をS00で確定し、dev-only migrationを追加しない |
| Owner selection | 下記の9フィールドすべてが承認済み。未記入の根拠やriskを補ってserializeしない |

prior Resultと方向/位置/Historyをfailure/Cancel/stale/superseded時に保つ。underlying errorとeditable draftを保持する。Pythonのfinal facts owner、単一resolverと単一artifact transaction、orientation-independent planを保つ。

architecture ratchet、canonical rendering、SVG sanitizer/admission、Worker資源管理、History readiness/finalization、scientific geometryのfailureはProduct選択で免除されない。BUG-17の新しいlinear crop/rotation代替は対象外。

## Implementation route and evidence

Route: **DURABLE_AUTHORITY_REQUIRED**。ownerの選択は完了しているが、当該authorityを含むdevをruntimeのbaseへ取り込むまではcandidate runtimeを認可しない。決定更新commitはContractだけを変更する。実装は [総合計画書](../01_MASTER_IMPLEMENTATION_PLAN.md) と対応する独立session promptに従う。

本Pack作成時点でruntimeは未実装。必要なevidenceは上記A条件に対応するsource-bound方向・位置・履歴・操作の実証であり、将来の成功を記録済み結果として扱わない。session evidenceは [保存方針](../evidence/README.md) に従う。受け入れ後も関心と責任を分離し、別concernのscopeをこのPackで広げない。

## Exact Product Decision Owner receipt

Supersedes PD-OI-034 revision 4. Receipt UTF-8 SHA-256: `280894a4c5bab41a1ae9309cf9f30c1bbd5aaff86541d4a0bc14cad61258d489`.

```text
PRODUCT_DECISION
Concern: web.similarity-alignment.choice-and-retry
Scenario revision: 5
Choice: A / LOCAL_EXCLUSIVE_DIRECTION_REVIEW_WITH_RETRY
Rationale: One exclusive direction selection should govern local candidate and direction previews, and retry must retain user intent without combining a global policy with per-row overrides or committing an unseen reference reversal.
Must preserve: Local no-Worker mode/custom/candidate/Skip editing; one direction tagged-union state with custom row values only in Custom; Python ownership of candidate eligibility and final source/display facts; one final batch validation per Apply attempt; one resolver for visible preview and final absolute orientations/reference-center placement; unchanged unknown/skipped/missing/unusable targets with reasons; an updated review and another Apply if final validated output differs; editable mode/custom choices after validation or render failure with the underlying message; prior Result, directions, placement and History after failed, canceled, stale or superseded work; canvas interaction, retained initially-Keep review after automatic render failure, atomic artifact Undo/Redo, Session/regeneration and the independently accepted Reset contract. Policies remain transient, and successful committed direction changes are actual record-state deltas, including the reference when changed.
May retire: The single draft-level Match flag, reference-relative direction as the only bulk operation, and post-Apply-only per-record exceptions. Do not add a second validation, rendering or History path or persisted direction policy.
Accepted residual risk: A final changed direction/placement receipt may require another Apply. Keep prior artifacts and local intent, show the new preview and limit re-review to actual output differences. Unknown directions are never guessed and stale reference/source bindings reject explicitly.
Owner: satoshikawato
Decision date: 2026-09-26
```
