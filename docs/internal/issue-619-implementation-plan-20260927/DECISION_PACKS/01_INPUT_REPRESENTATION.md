# Product Decision Pack 01 — 数値・単位入力と percent の表示

Status: **unsigned proposal**。推奨Aのresponseは本文記入済み。Owner欄への署名は全文への明示同意として扱う。記入日が異なる場合はDecision dateも更新する。未署名の内容はProduct authorityではない。

## Identity

- Concern key: `tracks.circular-measure-input-representation`
- Scenario revision: `1`
- Discovery lane: developer preflight
- Prepared from base SHA: `88028fd242d263f0fe86aaf9da57b8dc9eb082f6`
- Prepared for: `fix/issue-619-circular-track-measure-inputs`。runtime headは未作成。
- Prepared by: Codex, 2026-09-27
- Related: Issue #619 / [master plan](../MASTER_PLAN.md)

## Trigger and user journey

利用者が数値1.5の意味を selector から一意に読めるよう、width/radius を二つの controls に分ける。percent を通常の選択肢にもするか、倍率表示にまとめるかを決める。

既存 Session の保存値を読む → 数値／unit を編集または paste → Generate → Save/Load。checkpoints は Load、入力確定、Generate、History、round trip と download。

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

## Choice A — NUMERIC_PX_FACTOR_WITH_LEGACY_INPUT（推奨） / Choice B — NUMERIC_PX_FACTOR_PERCENT

| Dimension | Choice A | Choice B |
| --- | --- | --- |
| Complete normative outcome | width/radius とも数値 text input と px/×R selector。R は基準円半径。typed px/factor と旧 bare/px/% を受理。65% は0.65＋×Rとして読むが、Loadだけでは raw scalar を書き換えない。有効な単位付き入力/pasteは共有adapterが数値＋unitへ取り込む。不正・未完成textは保持し、validityを広げない。空欄はAuto。 | width/radius とも数値 text input と px/×R/% selector。65% は65＋%として読み・編集できる。%はrequestでfactor0.65へ投影。旧bare/px/%とtyped objectを受理し、不正draftを保持する。空欄はAuto。 |
| Preserved effects | 既存の値・unit・draft/Result・request・privacy・recovery | 同じ不変条件 |
| Added effects | 数値とpx/×Rの明示選択、有効suffixの取り込み、percentの倍率表示 | 数値とpx/×R/%の明示選択、literal percentの数値編集 |
| Lost / retired effects | 通常の編集表示における literal percent spelling と、数値欄内の単位の恒常表示。値・percentによる入力・Session受理は退役しない。 | 数値欄内の単位の恒常表示のみ。percentの通常表示は保つ。 |
| Discoverability / accessibility | numeric/unitのaccessible name、help、keyboard操作 | 同左。disabled時は理由と次操作も提示 |
| Canonical state update | owner actionから既存scalarへ1回。Autoはnull | 同左。Bの新preferenceがある場合もrequestとは分離 |
| Undo / Redo | manual numeric/unitの通常History。AutoはPack03 | 同左。Bの追加保存はPack03に従う |
| Session / regeneration | 保存draftとcommitted Resultを保ち、次回Generateで反映 | 同左。Bの追加metadataは既存namespaceのevidenceが必要 |
| Export / scientific output | 現在のResultを出力。canonical pairの意味を維持 | 同左。換算を選ぶ場合のprecisionはevidenceで確認 |
| Validation / error / recovery | invaliddraftを消さずrow error。旧Resultを保ち訂正可能 | 同左。disable時も理由と訂正/Clearの継続を残す |
| Cache / provenance | 新cache/provenanceなし | 新cache/provenanceなし。%とfactorの表現換算にRは不要 |
| Performance / resources | UI-only。追加Worker/renderなし | 必要な補助情報がある場合はevidenceでoperationを測る |
| Compatibility | 旧valid Session/scalarを受理。型とscientific meaningを保持 | 同左。追加metadataにreader証拠が必要 |
| Architecture | 既存scalar/editor/request/History owners。不要なparallel pathなし | 同じconstraints。追加owner/schemaにはratchet適用 |
| Evidence available / missing | base defectと既存scalar grammar確認済み。codec/controls/roundtripは検証が必要 | 同左。三択とpercent/factorの表現往復を追加検証 |
| Residual risk | percent入力を倍率表示へまとめるため、65%が0.65と読めることをhelpで説明する必要がある。suffix入力の途中と確定を区別する。 | 三択となり、factor↔percentにも100倍の表現差がある。Auto、invaliddraft、unit変更の規則との組合せを明示する必要がある。 |
| Route | DURABLE_AUTHORITY_REQUIRED | DURABLE_AUTHORITY_REQUIRED |
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
| numeric1.5 | 選択したpxまたは×R | 選択したpx/×R/% |
| 保存65% | 0.65＋×R、意味維持 | 65＋%、意味維持 |
| 選択肢数 | 2 | 3 |

## Evidence-first route and limits

Aの入力・request・Save/Load・Historyに使うWeb draft表現はS00でdeterministic disposable checksを行う。fixtureはtyped px/factor、旧px/%、text valueのtyped Web draft、null、invalid、disabled/inactive、settings-only、CLI-origin。実行command・raw observation・境界の受理結果を記録する。

Bの%表示はfactorとは100倍の表示差がある。保存値・selector・requestが同じ意味へ戻ることを検証する。evidence-only作業ではdefault/runtime selection、Product authority、request schema、reference baselineを変更しない。証拠が不足したまま署名でwaiveしない。

## Engineering recommendation

数値の意味を明示しつつ、通常操作の選択肢をpxと倍率の二つに絞る。percentによる既存入力と保存値の意味は維持する。

**この推奨はProduct authorityではない。** 他Packの署名を本Packの承認として使わない。Aが署名された場合だけAのresponseをserializeする。Bへの変更はBの完全なresponseと必要evidenceを用意する。

## Product Decision Owner response — A

```text
PRODUCT_DECISION
Concern: tracks.circular-measure-input-representation
Scenario revision: 1
Choice: A / NUMERIC_PX_FACTOR_WITH_LEGACY_INPUT
Rationale: 数値の意味を明示しつつ、通常操作の選択肢をpxと倍率の二つに絞る。percentによる既存入力と保存値の意味は維持する。
Must preserve: 既存px/factor/%の値とunit、precision、typed request/Session、Auto、invaliddraft、draft/Result分離、適用時点、History、privacy、失敗復旧。
May retire: 数値欄内に単位を恒常表示する旧UIと、percentをliteral spellingのまま通常表示することのみ。percent入力や既存Sessionの受理は退役しない。
Accepted residual risk: percent入力を倍率表示へまとめるため、65%が0.65と読めることをhelpで説明する必要がある。suffix入力の途中と確定を区別する。
Owner: ____________________ （maintainer loginによる署名）
Decision date: 2026-09-27
```

Responseを得たら既存static Product Contractへ全文を正確にserializeし、machine representationを提示する。AはPRODUCT_CHANGEのdurable候補なのでPR-local exact-head blockでは承認しない。authority-only内容がorigin/devへ入る前に依存runtimeを実装しない。
