# Product Decision Pack 02 — 単位変更の意味

Status: **Choice A signed**。`satoshikawato` が推奨 A の全文を `2026-09-27` に明示署名。[PR #621](https://github.com/satoshikawato/gbdraw/pull/621) で既存 Product Contract revision 23 へ統合済み。正本は[既存契約](https://github.com/satoshikawato/gbdraw/blob/252986d096011fcf1a0f5564e940480d3b92844d/docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md)。本 Pack は判断の原文を保持する参照文書。

## Identity

- Concern key: `tracks.circular-measure-unit-change`
- Scenario revision: `1`
- Discovery lane: developer preflight
- Prepared from base SHA: `88028fd242d263f0fe86aaf9da57b8dc9eb082f6`
- Prepared for: `fix/issue-619-circular-track-measure-inputs`。runtime headは未作成。
- Prepared by: Codex, 2026-09-27
- Related: Issue #619 / [master plan](../MASTER_PLAN.md)

## Trigger and user journey

pxと倍率を別controlにすると、unit changeが「数値を維持する編集」か「geometryを維持する表現変更」かを明示する必要がある。

manual valueを編集 → unit selectorを変更 → Pendingを確認 → Generate → Undo/Redo。R未解決・異なるRのgrid/batch・invaliddraftでも利用可能な次操作を明示する。

Current behavior: native text inputが保存scalar objectを直接受けて[object Object]になる。従来のbare/px/% grammarと保存draft優先はコード・fixtureで確認できる。表示不具合は科学的値の破損を意味しない。

## Non-waivable constraints

既存 px/factor/% の科学的意味、唯一の canonical request/Worker/Result admission、Load の保存 preview と draft 優先、失敗時の旧 Result、既存 History、local-only/privacy を保つ。不正値を auto/0 にしない。canonical scalar に text/% unit を追加しない。Product receipt は architecture/scientific/required-evidence failure を waive しない。

## Authority search

| Source | Result / gap |
| --- | --- |
| base Product Impact map | request/Result owner/path の関心は登録済み。本 control の具体 UX は決めていない |
| base durable BD / current decision | tools/web-product-decisions.json の decisions は空。対象の BD と exact-head decision はない |
| static Product Contract | PD-OI-037: 適用時点・draft/Result・History/Session/recovery。PD-OI-043: Circular width/radius の factor/% 維持。具体的な selector 動作は未規定 |
| base Web contract | services/session-request.js の唯一の projection、既存 History、privacy/Worker constraints。具体的な selector 動作は未規定 |
| released compatibility | 0.13.0 の ScalarSpec: bare factor、px suffix、%→factor。現在 Session 44/request 8 に同じ意味がある |
| code/tests | 原因と既存受理の証拠。独立 Product authority ではない |
| domain/integrity | px と factor は異なる量。既知の数値・unit を表示や切替で黙示変換しない |

Result: UNRESOLVED（本 Pack の詳細）。Procedural classification: PRODUCT_DECISION_REQUIRED。既存実装の表示不具合の存在は、新しい UX 詳細の承認を代替しない。

## Choice A — KEEP_NUMBER_CHANGE_UNIT（推奨） / Choice B — PRESERVE_RESOLVED_GEOMETRY

| Dimension | Choice A | Choice B |
| --- | --- | --- |
| Complete normative outcome | unit changeは現在表示されている数値textを保持してunitだけ変える。1.5×R→1.5px。半径による換算や自動Generateをしない。nonemptyの変更は1回の通常History transaction。不正numerictextも保持する。同じunitはno-op。Autoの次回入力用unitはPack03で扱う。 | 有効なmanual valueは図上の大きさが同じになるよう換算。R=100pxなら1.5×R→150px。正しいRが不明、複数instanceで異なる、staleまたはinvaliddraftならunit切替をdisableし理由を表示。利用者は数値訂正またはClear→unit選択→新数値入力へ進める。換算だけでGenerateしない。 |
| Preserved effects | 既存の値・unit・draft/Result・request・privacy・recovery | 同じ不変条件 |
| Added effects | 数値維持のunit edit。R未解決でも利用可能 | geometry維持のunit換算。変換不可時の理由とClear経由の継続 |
| Lost / retired effects | なし。新selectorの意味を数値維持の編集として定義する。旧scalarの受理や科学的意味は維持。 | R未解決／異なるR／invaliddraftで、値を保持したまま単位を切り替えられる直接操作。Clear/訂正の継続は残す。 |
| Discoverability / accessibility | numeric/unitのaccessible name、help、keyboard操作 | 同左。disabled時は理由と次操作も提示 |
| Canonical state update | owner actionから既存scalarへ1回。Autoはnull | 同左。Bの新preferenceがある場合もrequestとは分離 |
| Undo / Redo | manual numeric/unitの通常History。AutoはPack03 | 同左。Bの追加保存はPack03に従う |
| Session / regeneration | 保存draftとcommitted Resultを保ち、次回Generateで反映 | 換算後のscalarを既存draftへ保存し、次回Generateで反映。conversion policyを別schemaに保存しない |
| Export / scientific output | 現在のResultを出力。canonical pairの意味を維持 | 同左。換算を選ぶ場合のprecisionはevidenceで確認 |
| Validation / error / recovery | invaliddraftを消さずrow error。旧Resultを保ち訂正可能 | 同左。disable時も理由と訂正/Clearの継続を残す |
| Cache / provenance | 新cacheやrender provenanceを増やさない | 換算が必要な場合だけ現在のgeometry identityを検証 |
| Performance / resources | UI-only。追加Worker/renderなし | 必要な補助情報がある場合はevidenceでoperationを測る |
| Compatibility | 旧valid Session/scalarを受理。型とscientific meaningを保持 | 同左。unit-change availability/precisionのevidenceが必要 |
| Architecture | 既存scalar/editor/request/History owners。不要なparallel pathなし | 同じconstraints。追加owner/schemaにはratchet適用 |
| Evidence available / missing | base defectと既存scalar grammarを確認済み。new UIは実装後検証が必要 | new UIに加えBの換算／永続namespaceのevidenceが必要な場合あり |
| Residual risk | 倍率からpxへ切り替えると図上の大きさが変わる。helpとPendingに、数値維持・次回Generate反映を明記する。 | R取得、複数record、stale Result、丸め誤差により換算結果を誤る可能性。代表geometryと変換可能条件の証拠が必要。 |
| Route | DURABLE_AUTHORITY_REQUIRED | EVIDENCE_REQUIRED |
| Next action | 署名→existing static authorityへ正確にserialize→authority-only統合→runtime実装 | DURABLEなら同じauthority手順。EVIDENCEなら先にruntime不変のevidence-only確認 |

Pack01-Bで%も選択可能なら、Aは65%→65×Rのように表示数値を維持する。Bは%/factor間を100倍の表現差で換算し、pxとの換算に正しいRを使う。旧percentのcanonical意味はどちらでも変えない。

## Comparison and realization requirements

| Independent requirement | A | B |
| --- | --- | --- |
| scalar value/unit integrity | 必須 | 必須 |
| separate numeric/unit controls | 必須 | 必須 |
| Auto/invaliddraft/recovery | 必須 | 必須 |
| History/Session/Generate/export continuity | 必須 | 必須 |
| 本Pack固有のoutcome | 上表のAを完全に実現 | 上表のBを完全に実現 |

これらはAND条件。numeric/unit controlsだけを実装して、AutoやHistory/Sessionを欠いた状態を完了扱いにしない。複数実装が同じoutcomeを実現する場合は実装判断として扱い、Productの別選択肢にしない。

| Observable difference | A | B |
| --- | --- | --- |
| 1.5×R→px | 1.5px | R=100pxなら150px |
| R未解決 | 通常操作可能 | manual切替は理由付きdisable、Clear経由可 |
| grid/batch | 同じscalar編集 | 全instanceのR一致など変換可能条件が必要 |

## Evidence-first route and limits

Aの入力・request・Save/Load・Historyに使うWeb draft表現はS00でdeterministic disposable checksを行う。fixtureはtyped px/factor、旧px/%、text valueのtyped Web draft、null、invalid、disabled/inactive、settings-only、CLI-origin。実行command・raw observation・境界の受理結果を記録する。

Bは換算用Rの正本、全instanceの一致、未生成/stale/invalid時のavailability、conversion precisionを先に確認する。evidence-only作業ではdefault/runtime selection、Product authority、request schema、reference baselineを変更しない。証拠が不足したまま署名でwaiveしない。

## Engineering recommendation

単位選択を、入力した数値の意味を明示的に変更する編集として統一する。現在の円半径や古いResultに依存せず、生成前や複数recordでも同じ操作を使える。

**この推奨はProduct authorityではない。** 他Packの署名を本Packの承認として使わない。Aが署名された場合だけAのresponseをserializeする。Bへの変更はBの完全なresponseと必要evidenceを用意する。

## Product Decision Owner response — A

```text
PRODUCT_DECISION
Concern: tracks.circular-measure-unit-change
Scenario revision: 1
Choice: A / KEEP_NUMBER_CHANGE_UNIT
Rationale: 単位選択を、入力した数値の意味を明示的に変更する編集として統一する。現在の円半径や古いResultに依存せず、生成前や複数recordでも同じ操作を使える。
Must preserve: numericdraftとunitの明示、Auto、Generateまで旧Resultを保つこと、History/Session、失敗復旧、既存scalarの科学的意味、local-only。
May retire: なし。
Accepted residual risk: 倍率からpxへ切り替えると図上の大きさが変わる。helpとPendingに、数値維持・次回Generate反映を明記する。
Owner: satoshikawato
Decision date: 2026-09-27
```

Responseを得たら既存static Product Contractへ全文を正確にserializeし、machine representationを提示する。AはPRODUCT_CHANGEのdurable候補なのでPR-local exact-head blockでは承認しない。authority-only内容がorigin/devへ入る前に依存runtimeを実装しない。
