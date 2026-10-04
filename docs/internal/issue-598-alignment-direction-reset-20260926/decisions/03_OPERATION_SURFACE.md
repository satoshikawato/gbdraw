# Accepted Product Decision Pack: PD-OI-031 — 操作面と自動Align

## Identity and approval

- Author / Product Decision Owner: `satoshikawato`
- Decision date: 2026-09-26
- Concern: `diagram-generation.similarity-alignment.surface-scope`
- Scenario revision: 5
- Status: ACCEPTED。選択、根拠、保持条件、廃止範囲、残余リスクは下記の完全な承認文で確定。
- Prepared from base SHA: `d457b7189b137185a8dec800819a312c30b969fa`
- Issue: https://github.com/satoshikawato/gbdraw/issues/598
- Authority update: `product/issue-598-decisions-20260926` / `eba01dd518ff6edd4effb4b00d519cbd78ff8623`
- Implementation: `fix/issue-598-alignment-direction-reset-20260926`
- 実装責任: S02 / S03

このPackは当該concernだけの承認記録であり、別concernの承認や実装を代替しない。正式authorityは [Option Integrity Product Contract](../../OPTION_INTEGRITY_PRODUCT_CONTRACT.md)。この文書はmachine decision storeではない。

## User journey and outcome

利用者は複数ゲノムの線形図を操作する。目標は、正確なreferenceを選び、通常の自動Alignまたは明示reviewから適用する。 選ぶreferenceは正確なfeature identityであり、多数派方向の代表とは限らない。生物学的sourceと図の表示方向は別の情報である。

resolvedなdefault AlignはKeepで自動適用し、不要なreviewを強制しない。ambiguityまたは明示reviewでは一択の方向選択、candidate factsとSelect/Skip、結果矢印を提供する。

checkpointは操作前preview、final validation、成功artifact、Session fresh Load、regeneration、Undo/Redo、failure recoveryである。このPackの受入条件は総合計画書の **A01–A03、A13、A14**。独立要件をすべて満たし、共通choice IDによって不足を代替しない。

## Authority search and non-waivable boundaries

| 調査対象 | 結果と扱い |
| --- | --- |
| Product map / contract | base OIPCの当該concernをscenario supersessionで更新する。mapやcheckerは変更しない |
| Durable BD registry | 調査baseに当該変更を認可するactive BDはなく、新しいBD storeは作らない |
| Source / tests | 現行Match checkboxと位置のみResetの観測根拠。Product authorityそのものとは扱わない |
| Domain / scientific integrity | source identity/strandを保持し、同一record transformと読みやすいtextを守る |
| Persisted compatibility | 公開済みSessionの証拠をS00で確定し、dev-only migrationを追加しない |
| Owner selection | 下記の9フィールドすべてが承認済み。未記入の根拠やriskを補ってserializeしない |

exact popup/drawer reference、Python eligibility/final validation、CLI defaultsとambiguity拒否、underlying error、ordinary record controls、keyboardと既存狭幅palette契約を保つ。

architecture ratchet、canonical rendering、SVG sanitizer/admission、Worker資源管理、History readiness/finalization、scientific geometryのfailureはProduct選択で免除されない。BUG-17の新しいlinear crop/rotation代替は対象外。

## Implementation route and evidence

Route: **DURABLE_AUTHORITY_REQUIRED**。ownerの選択は完了しているが、当該authorityを含むdevをruntimeのbaseへ取り込むまではcandidate runtimeを認可しない。決定更新commitはContractだけを変更する。実装は [総合計画書](../01_MASTER_IMPLEMENTATION_PLAN.md) と対応する独立session promptに従う。

本Pack作成時点でruntimeは未実装。必要なevidenceは上記A条件に対応するsource-bound方向・位置・履歴・操作の実証であり、将来の成功を記録済み結果として扱わない。session evidenceは [保存方針](../evidence/README.md) に従う。受け入れ後も関心と責任を分離し、別concernのscopeをこのPackで広げない。

## Exact Product Decision Owner receipt

Supersedes PD-OI-031 revision 4. Receipt UTF-8 SHA-256: `a10810c7c0df8d912d1d14c9200056e42968a07f49e50011f90b4d1161b157f7`.

```text
PRODUCT_DECISION
Concern: diagram-generation.similarity-alignment.surface-scope
Scenario revision: 5
Choice: A / AUTO_APPLY_WITH_EXPLICIT_DIRECTION_REVIEW
Rationale: Keep automatic resolved alignment uncomplicated while opened reviews expose mutually exclusive final display direction outcomes, including reversing only a minority reference through an all-right or all-left choice.
Must preserve: Exact popup/drawer reference selection; automatic Python-resolved default alignment in Keep mode without forced review; an accessible ambiguity-required or explicit review with candidate facts, Select/Skip, resolution summary, one exclusive Keep/right/left/Custom direction selection and per-record resulting arrows; clear known-strand selected-anchor scope including the reference and unknown/skipped exclusion; shared typed validation, actionable underlying errors, strict CLI ambiguity rejection, existing CLI/API defaults, ordinary record controls, keyboard operation and the existing accepted narrow-screen palette coverage limitation. Apply persists absolute record transforms and placements, never review policies in the alignment plan.
May retire: An opened review exposing only a single reference-relative Match checkbox, permanent preservation of reference direction during an explicit direction operation, and mandatory post-Apply correction for per-record exceptions. Keep automatic/default and explicit review entry points.
Accepted residual risk: Explicit direction choices have more outcomes than the default path and can move the reference record's left edge while its selected feature center stays fixed. Use arrow-based labels, a single radio group, truthful target previews and Custom disclosure. This decision does not redesign narrow-screen palette coverage or introduce new CLI direction flags.
Owner: satoshikawato
Decision date: 2026-09-26
```
