# 承認済み Decision Pack — 装飾配置の継承

Concern: `web.composition-decoration-continuity`。Scenario revision: `1`。
選択: `A / CARRY-MATCHED-DECORATION-DELTAS`。Product Decision Owner: `satoshikawato`。
承認日: 2026-09-26。状態: **製品結果は明示承認済み、恒久的 base authority への記録は S00 の作業**。

## 承認の範囲と根拠

対象は [Issue #599](https://github.com/satoshikawato/gbdraw/issues/599) のこの concern のみ。
Product Decision Owner は、装飾配置の継承、Layout edit の発見、検索・toolbar の配置という3つの独立した推奨 A に対して「すべて推奨案で承認します。」と明示回答した。
この Pack は、そのうち `web.composition-decoration-continuity` の完全な A を記録する。他 concern の承認記録は [一覧](00_APPROVED_PRODUCT_DECISIONS.md) から参照する。選択 B は記録しない。
電子的な明示承認として記録し、架空の手書き署名は作らない。

調査基準: `origin/dev@d457b7189b137185a8dec800819a312c30b969fa`。
この基準では active BD は存在せず、OIPC revision 16 が intent/failure 等を規定するが、この完全な outcome は未登録だった。
developer preflight の Product 判断は上記承認で解決済み。恒久記録の経路は **DURABLE_AUTHORITY_REQUIRED**。
S00 は最新 dev の map/契約/accepted decisions を再調査し、既存 authority と重複しない記録先と番号を確定する。
ここに仮 BD/PD-OI 番号や runtime 自己承認は置かない。

## 選択した完全な製品結果

| Dimension | 承認済み A の outcome |
| --- | --- |
| Complete normative outcome | 同じ図の legend/title/Linear scale の delta を新 automatic 配置へ1回加算。mode/grouping、validated source/region、record identity で照合。prefix/配列順/DOM順を使わず、batchの各出力も別々に対応。非ゼロ target の消失・未知対応は候補を公開せず、旧 Result を保持して対象 Reset または設定修正を案内。zero/fresh は自動 Generate |
| Preserved effects | 同じ request path、record/alignment semantics、drag/Reset、History/Session/current Result export、failure isolation。diagram全体/個別record/padding/legend順の新継承保証は含めない |
| Added effects | 通常 Generate / committed-candidate render / automatic reflow で対象継承。対応不能から明示的に回復 |
| Lost effects | Generate による無条件 reset。対応不能の場合は Reset/修正の追加手順 |
| Retired effects | 通常 Generate が対象の非ゼロ delta を無言で捨てる動作だけ |
| Discoverability/accessibility / immediate feedback | 通常は自動継承。対応不能の理由・対象・Reset/修正を keyboard/touch から到達可能に表示 |
| Canonical state update | Result SVG を差分の正本とし、transaction-local snapshot を candidate に適用。UI refs は同期値。latent/global map なし |
| Undo/Redo | Generate は継承を含む1 replacement。保存 Result を復元し delta を再加算しない |
| Session / regeneration | 保存 Result から次の Generate の差分を取得。通常 load は保存 Result を維持。未知対応は無言破棄しない |
| Export/artifact | 適用済みの current Result を各形式へ出力。raw Python recipeだけで手動位置を再現する保証なし |
| Validation/error | finite delta、一意 target、図の同一性を検査。source/region/mode/grouping/record集合変更や欠落/重複/未知が転用不能なら候補公開前のエラー |
| Failure/recovery / next available action | render/transform/bind失敗、Cancel/staleは旧Result/request/History保持。対応不能は対象 Reset または設定修正→Generate |
| Scientific-output | 装飾位置のみ。recordTranslations/active alignment の record delta を二重加算しない |
| Cache/provenance | 既存 validated identity/digest 使用。raw search/cache key に delta を追加しない。candidate と保存Resultを一致 |
| Performance | 非ゼロ対象のみ。候補の既存parse/serialize共用、batch旧SVGは必要分のみ、zero fast path維持 |
| Compatibility | writer/readerを維持。照合できない保存図は明示回復 |
| Architecture | composition が演算 owner。既存 candidate transform seam に注入。第二 replacement pathなし |
| Evidence available/missing | base観測/source。C01–C08/R01の新動作、batch/failure注入は未実装・未実行 |
| Residual risk | 新 automatic + 同じ delta なので絶対位置は変わり、clipping/overlapが残りうる。自動clampせず padding/Reset で調整。対応不能では明示 Reset/修正が必要 |
| Route | DURABLE_AUTHORITY_REQUIRED |

## 承認本文

以下が承認された回答の全文である。owner/date/rationale/preservation/retirement/risk を補完・拡張せず恒久契約へ転記する。

```text
PRODUCT_DECISION
Concern: web.composition-decoration-continuity
Scenario revision: 1
Choice: A / CARRY-MATCHED-DECORATION-DELTAS
Rationale: 色やfontを直すたびに装飾の配置をやり直す負担をなくし、別の図に位置を誤転用しない。
Must preserve: legend/title/Linear scaleのdrag・適用・Reset・History・Session・current Result export、source/region/record同一性、既存record/alignmentの意味、zero fast path、失敗/Cancel/stale時の旧Resultとcommitted requestを保持する。未知/欠落からの無言削除、別sourceへの誤転用、record deltaの二重加算、新schema/Worker/全History cloneを認めない。通常Generate/committed-candidate/automatic reflowに同じ候補境界を使い、batch全出力を対応identityへだけ適用する。
May retire: 通常Generateがlegend/title/Linear scaleの非ゼロdeltaを無言で捨てる動作だけ。diagram全体、個別record、padding、legend順の新しい継承保証は含めない。
Accepted residual risk: 新automatic配置に同じdeltaを加えるので絶対位置は変わり、clipping/overlapが残りうる。自動clampせずpadding/Resetで調整する。対応不能は候補公開前に止まり、明示Reset/設定修正が必要。
Owner: satoshikawato
Decision date: 2026-09-26
```

承認本文 UTF-8（末尾 newline を含めない）の SHA-256: `c7668d708ac6cd90b4373df0a685e418f21be7e653e99cbe899b82fba25e5e7d`。

## レビュー用の機械表現

これは上記本文のフィールド対応を見せる非実行 JSON であり、checker が読む新 schema/registry ではない。
unmapped concern の実行可能な decision store は追加しない。active authority は既存規約の記録先へ S00 で記録し、dev merge を確認する。

```json
{
  "concern": "web.composition-decoration-continuity",
  "scenarioRevision": 1,
  "choice": "A / CARRY-MATCHED-DECORATION-DELTAS",
  "rationale": "色やfontを直すたびに装飾の配置をやり直す負担をなくし、別の図に位置を誤転用しない。",
  "mustPreserve": "legend/title/Linear scaleのdrag・適用・Reset・History・Session・current Result export、source/region/record同一性、既存record/alignmentの意味、zero fast path、失敗/Cancel/stale時の旧Resultとcommitted requestを保持する。未知/欠落からの無言削除、別sourceへの誤転用、record deltaの二重加算、新schema/Worker/全History cloneを認めない。通常Generate/committed-candidate/automatic reflowに同じ候補境界を使い、batch全出力を対応identityへだけ適用する。",
  "mayRetire": "通常Generateがlegend/title/Linear scaleの非ゼロdeltaを無言で捨てる動作だけ。diagram全体、個別record、padding、legend順の新しい継承保証は含めない。",
  "acceptedResidualRisk": "新automatic配置に同じdeltaを加えるので絶対位置は変わり、clipping/overlapが残りうる。自動clampせずpadding/Resetで調整する。対応不能は候補公開前に止まり、明示Reset/設定修正が必要。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

## 非免除条件と実装前提

既存 typed request、SVG sanitization、科学的な値と record/alignment identity、保存 Result/Session/History/export、failure isolation、lazy Worker、性能 gate を維持する。
承認は gate failure、追加 dependency、未知値の無言破棄、別 schema/互換経路の増加を免除しない。
[総合計画](MASTER_PLAN.md) の受入条件と各 concern の結果を照合する。候補の同じ選択 ID だけで、全 preservation requirement が成立したとはみなさない。

依存する runtime の実装前に恒久 authority が dev に必要である。
S00 の変更は authority-only とし、runtime は含めない。PR 作成・dev merge の権限境界は総合計画に従う。
