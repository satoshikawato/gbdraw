# Product Decision Pack 03 — Auto の unit preference の寿命

Status: **Choice A signed**。`satoshikawato` が推奨 A の全文を `2026-09-27` に明示署名。[PR #621](https://github.com/satoshikawato/gbdraw/pull/621) で既存 Product Contract revision 23 へ統合済み。正本は[既存契約](https://github.com/satoshikawato/gbdraw/blob/252986d096011fcf1a0f5564e940480d3b92844d/docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md)。本 Pack は判断の原文を保持する参照文書。

## Identity

- Concern key: `tracks.circular-measure-auto-unit-lifecycle`
- Scenario revision: `1`
- Discovery lane: developer preflight
- Prepared from base SHA: `88028fd242d263f0fe86aaf9da57b8dc9eb082f6`
- Prepared for: `fix/issue-619-circular-track-measure-inputs`。runtime headは未作成。
- Prepared by: Codex, 2026-09-27
- Related: Issue #619 / [master plan](../MASTER_PLAN.md)

## Trigger and user journey

空欄はunitを持たないAutoである。一方、利用者はpxを先に選んでから数値を入力できる。空欄時の選択をHistoryやSessionへ保存するかを決める。

空欄でunitを選ぶ → 数値を入力 → Clear → panel再マウント／Load／Reset／Undo。manual値のunitとAuto時の入力preferenceを区別する。

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

## Choice A — TRANSIENT_AUTO_UNIT（推奨） / Choice B — PERSIST_AUTO_UNIT_PREFERENCE

| Dimension | Choice A | Choice B |
| --- | --- | --- |
| Complete normative outcome | 空欄は常にcanonical null。初期unitは×R。空欄でpxを選べば、componentが存続する間は次の数値がpxになる。Autoのunit選択だけはtransient UIで、History/Session/request/Resultを変更しない。panel再マウント、Load、Resetで既定へ戻る。manual stateのHistoryはtext/unitを正確に復元し、Autoへの復元はnullを復元するがunit preferenceの復元は約束しない。 | 空欄はcanonical null。初期unitは×R。Autoのunit preferenceをslot identity＋fieldで保持し、HistoryとSessionに保存する。再マウントとLoadで同じ選択へ復元し、Resetで既定へ戻す。canonical render requestにはpreferenceを入れない。既存Sessionで欠落していれば×Rを補う。 |
| Preserved effects | 既存の値・unit・draft/Result・request・privacy・recovery | 同じ不変条件 |
| Added effects | component存続中の次回入力unit選択 | unit preferenceのHistory/Session/re-mount復元 |
| Lost / retired effects | なし。Autoのunit preferenceに新しい永続化保証を設けない。manual unitと既存Sessionは維持。 | なし。新しい永続UI preferenceを追加する。 |
| Discoverability / accessibility | numeric/unitのaccessible name、help、keyboard操作 | 同左。disabled時は理由と次操作も提示 |
| Canonical state update | Auto unitだけならtransient preferenceのみ。scalarはnullのまま | Auto unitだけならUI preference ownerを更新。scalar/requestはnullのまま |
| Undo / Redo | manual scalarは通常History。Auto unitはHistoryに入れない | manual scalarとAuto unit preferenceをそれぞれHistoryで復元 |
| Session / regeneration | manual scalarは保存。Auto unit preferenceは保存しない。再生成のAuto意味は同じ | manual scalarとAuto unit preferenceを保存。Auto requestはnullのまま |
| Export / scientific output | Auto preferenceは現在Result/downloadを変えない | 同左。保存したpreferenceもgeometryを変えない |
| Validation / error / recovery | invaliddraftを消さずrow error。旧Resultを保ち訂正可能 | 同左。disable時も理由と訂正/Clearの継続を残す |
| Cache / provenance | 新cache/provenanceなし | 新render cache/provenanceなし。preferenceだけをUIに保存 |
| Performance / resources | UI-only。追加Worker/renderなし | 必要な補助情報がある場合はevidenceでoperationを測る |
| Compatibility | 旧valid Session/scalarを受理。型とscientific meaningを保持 | 同左。追加metadataにreader証拠が必要 |
| Architecture | 既存scalar/editor/request/History owners。不要なparallel pathなし | 同じconstraints。追加owner/schemaにはratchet適用 |
| Evidence available / missing | base defectと既存scalar grammarを確認済み。new UIは実装後検証が必要 | new UIに加えBの換算／永続namespaceのevidenceが必要な場合あり |
| Residual risk | 空欄時だけのunit選択はpanel再マウントやLoadで忘れられる。manual scalarのunitは必ず残り、Auto geometryは変わらない。 | optional UI metadata、identity照合、削除／複製／Reset、Historyの所有が増える。readerが受理する保存namespaceと互換性の証拠が必要。 |
| Route | DURABLE_AUTHORITY_REQUIRED | EVIDENCE_REQUIRED |
| Next action | 署名→existing static authorityへ正確にserialize→authority-only統合→runtime実装 | DURABLEなら同じauthority手順。EVIDENCEなら先にruntime不変のevidence-only確認 |

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
| 空欄→px→1.5 | 1.5px | 1.5px |
| AutoでLoad/再マウント | 既定×Rへ戻る | 保存したunitへ戻る |
| Autoのunitだけの変更 | History/Sessionに入れない | History/Sessionへ入れる |

## Evidence-first route and limits

Aの入力・request・Save/Load・Historyに使うWeb draft表現はS00でdeterministic disposable checksを行う。fixtureはtyped px/factor、旧px/%、text valueのtyped Web draft、null、invalid、disabled/inactive、settings-only、CLI-origin。実行command・raw observation・境界の受理結果を記録する。

Bはoptional UI preferenceのreader/History/slot identity、copy/remove/reset、missing-field defaultを先に確認する。evidence-only作業ではdefault/runtime selection、Product authority、request schema、reference baselineを変更しない。証拠が不足したまま署名でwaiveしない。

## Engineering recommendation

Autoにgeometry上のunitはないため、その選択を次回入力用の小さなtransient preferenceとして扱う。manual値の意味と保存は保ち、追加の永続schemaやunit mirrorを避ける。

**この推奨はProduct authorityではない。** 他Packの署名を本Packの承認として使わない。Aが署名された場合だけAのresponseをserializeする。Bへの変更はBの完全なresponseと必要evidenceを用意する。

## Product Decision Owner response — A

```text
PRODUCT_DECISION
Concern: tracks.circular-measure-auto-unit-lifecycle
Scenario revision: 1
Choice: A / TRANSIENT_AUTO_UNIT
Rationale: Autoにgeometry上のunitはないため、その選択を次回入力用の小さなtransient preferenceとして扱う。manual値の意味と保存は保ち、追加の永続schemaやunit mirrorを避ける。
Must preserve: 空欄/Autoのnull意味、unitを先に選ぶ操作、manual値のunitとHistory/Session、既存preview/request、invaliddraftと失敗復旧。
May retire: なし。Autoのunit preferenceのHistory/Session保証は新設しない。
Accepted residual risk: 空欄時だけのunit選択はpanel再マウントやLoadで忘れられる。manual scalarのunitは必ず残り、Auto geometryは変わらない。
Owner: satoshikawato
Decision date: 2026-09-27
```

Responseを得たら既存static Product Contractへ全文を正確にserializeし、machine representationを提示する。AはPRODUCT_CHANGEのdurable候補なのでPR-local exact-head blockでは承認しない。authority-only内容がorigin/devへ入る前に依存runtimeを実装しない。
