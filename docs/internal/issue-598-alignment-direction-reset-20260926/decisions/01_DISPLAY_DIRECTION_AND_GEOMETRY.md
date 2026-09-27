# Accepted Product Decision Pack: PD-OI-027 — 表示方向と幾何変換

## Identity and approval

- Author / Product Decision Owner: `satoshikawato`
- Decision date: 2026-09-26
- Concern: `diagram-generation.similarity-alignment.transform-semantics`
- Scenario revision: 5
- Status: ACCEPTED。選択、根拠、保持条件、廃止範囲、残余リスクは下記の完全な承認文で確定。
- Prepared from base SHA: `d457b7189b137185a8dec800819a312c30b969fa`
- Issue: https://github.com/satoshikawato/gbdraw/issues/598
- Authority update: `product/issue-598-decisions-20260926` / `eba01dd518ff6edd4effb4b00d519cbd78ff8623`
- Implementation: `fix/issue-598-alignment-direction-reset-20260926`
- 実装責任: S01 / S02

このPackは当該concernだけの承認記録であり、別concernの承認や実装を代替しない。正式authorityは [Option Integrity Product Contract](../../OPTION_INTEGRITY_PRODUCT_CONTRACT.md)。この文書はmachine decision storeではない。

## User journey and outcome

利用者は複数ゲノムの線形図を操作する。目標は、selected referenceとtargetsの表示方向を決め、全アンカー中心を整列する。 選ぶreferenceは正確なfeature identityであり、多数派方向の代表とは限らない。生物学的sourceと図の表示方向は別の情報である。

自動Keepを保ち、明示reviewではKeep/right/left/Customを一択とする。referenceも方向変更の対象となり、そのfeature centerの直前canvas xと全yを保つ。

checkpointは操作前preview、final validation、成功artifact、Session fresh Load、regeneration、Undo/Redo、failure recoveryである。このPackの受入条件は総合計画書の **A01–A05、A14**。独立要件をすべて満たし、共通choice IDによって不足を代替しない。

## Authority search and non-waivable boundaries

| 調査対象 | 結果と扱い |
| --- | --- |
| Product map / contract | base OIPCの当該concernをscenario supersessionで更新する。mapやcheckerは変更しない |
| Durable BD registry | 調査baseに当該変更を認可するactive BDはなく、新しいBD storeは作らない |
| Source / tests | 現行Match checkboxと位置のみResetの観測根拠。Product authorityそのものとは扱わない |
| Domain / scientific integrity | source identity/strandを保持し、同一record transformと読みやすいtextを守る |
| Persisted compatibility | 公開済みSessionの証拠をS00で確定し、dev-only migrationを追加しない |
| Owner selection | 下記の9フィールドすべてが承認済み。未記入の根拠やriskを補ってserializeしない |

source strandと生物学的identityは変更しない。orientationはrecordが所有し、planやreceiptは描画方向のpolicyを持たない。text可読性、feature/label/annotation/ribbonの同一transform、exact anchor identityを保つ。

architecture ratchet、canonical rendering、SVG sanitizer/admission、Worker資源管理、History readiness/finalization、scientific geometryのfailureはProduct選択で免除されない。BUG-17の新しいlinear crop/rotation代替は対象外。

## Implementation route and evidence

Route: **DURABLE_AUTHORITY_REQUIRED**。ownerの選択は完了しているが、当該authorityを含むdevをruntimeのbaseへ取り込むまではcandidate runtimeを認可しない。決定更新commitはContractだけを変更する。実装は [総合計画書](../01_MASTER_IMPLEMENTATION_PLAN.md) と対応する独立session promptに従う。

本Pack作成時点でruntimeは未実装。必要なevidenceは上記A条件に対応するsource-bound方向・位置・履歴・操作の実証であり、将来の成功を記録済み結果として扱わない。session evidenceは [保存方針](../evidence/README.md) に従う。受け入れ後も関心と責任を分離し、別concernのscopeをこのPackで広げない。

## Exact Product Decision Owner receipt

Supersedes PD-OI-027 revision 4. Receipt UTF-8 SHA-256: `c2e9c2e37b8e54fc013737127db9a650acc48e6b324593fa3c162d595f42fc5d`.

```text
PRODUCT_DECISION
Concern: diagram-generation.similarity-alignment.transform-semantics
Scenario revision: 5
Choice: A / EXPLICIT_DISPLAY_DIRECTION_MODES
Rationale: Users should choose the final displayed direction of selected alignment features instead of making all targets follow a potentially minority-direction reference. Reference identity defines positioning, while record direction remains independent record state.
Must preserve: Default Keep current directions; one exclusive Keep, All selected features right-facing, All selected features left-facing or Custom direction mode; Custom per-record Keep/right/left choices; exact reference identity and selected target anchors; bulk direction scope including the reference and only known-strand selected anchors; explicit unchanged reasons for unknown directions and unchanged missing, unusable and skipped records; record-wide absolute orientation updates without editing biological source strands; a fixed pre-Align canvas x of the reference feature center, adjusted record placement as necessary, unchanged vertical placement and exact idempotent selected-anchor center alignment; per-record before/after arrow previews and truthful scope coverage; readable text and feature/label/annotation/ribbon geometry in the same record transform; one atomic validated orientations/placement/plan/Result commit; orientation-independent plans, ordinary Reverse, accepted Reset scope and complete artifact Undo/Redo. Changed final validated directions or reference-center placement require a refreshed preview and another Apply before commitment. Capture actual direction deltas including the reference if a separately authorized reset receipt is enabled.
May retire: The rule that alignment always preserves the reference record's direction and left-edge position; the single reference-relative Match checkbox; the restriction that target exceptions require post-Apply sidebar Reverse. Do not retire exact reference selection, its fixed anchor-center position, the default Keep mode or ordinary Reverse.
Accepted residual risk: Users can confuse display arrows with biological strand annotations or interpret all as every feature in a source. Show selected-anchor scope, reference participation, record names, before/after arrows and unknown exclusions; preserve sources and preview reference left-edge movement while its feature center stays fixed. Custom increases review density and requires keyboard/390 px acceptance. Majority inference, guessed unknown directions, source annotation changes, hidden flips and separate orientation/render owners are not accepted.
Owner: satoshikawato
Decision date: 2026-09-26
```
