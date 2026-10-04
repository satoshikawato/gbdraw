# 承認済み Decision Pack — 検索・toolbar の専用行

Concern: `web.preview-search-placement`。Scenario revision: `1`。
選択: `A / DOCKED-SEARCH-AND-CONTROLS`。Product Decision Owner: `satoshikawato`。
承認日: 2026-09-26。状態: **製品結果は明示承認済み、恒久的 base authority への記録は S00 の作業**。

## 承認の範囲と根拠

対象は [Issue #599](https://github.com/satoshikawato/gbdraw/issues/599) のこの concern のみ。
Product Decision Owner は、装飾配置の継承、Layout edit の発見、検索・toolbar の配置という3つの独立した推奨 A に対して「すべて推奨案で承認します。」と明示回答した。
この Pack は、そのうち `web.preview-search-placement` の完全な A を記録する。他 concern の承認記録は [一覧](00_APPROVED_PRODUCT_DECISIONS.md) から参照する。選択 B は記録しない。
電子的な明示承認として記録し、架空の手書き署名は作らない。

調査基準: `origin/dev@d457b7189b137185a8dec800819a312c30b969fa`。
この基準では active BD は存在せず、OIPC revision 16 が intent/failure 等を規定するが、この完全な outcome は未登録だった。
developer preflight の Product 判断は上記承認で解決済み。恒久記録の経路は **DURABLE_AUTHORITY_REQUIRED**。
S00 は最新 dev の map/契約/accepted decisions を再調査し、既存 authority と重複しない記録先と番号を確定する。
ここに仮 BD/PD-OI 番号や runtime 自己承認は置かない。

## 選択した完全な製品結果

| Dimension | 承認済み A の outcome |
| --- | --- |
| Complete normative outcome | searchを専用top row、toolbarを専用bottom row、canvasと同じeditorを中間workspaceへ配置。drawerがsearch/toolbarを覆わない構造にし、全幅で検索自由dragを退役。wrap/scrollで全機能を維持。通常高さ740px以上の受入条件はworkspace高さ200px以上、short viewport/keyboardでは全操作へscroll到達可能 |
| Preserved effects | 全search field/query/regex/Prev/Next/Open/Enter、active match/focus、toolbar、同じcanvas/editor、drawer/tab/Close/Escape、Session/Export/History、PD-OI-035 |
| Added effects | search/toolbar/drawerの非被覆を構造で保証、drawer/resizeでquery/focus維持 |
| Lost effects | searchを好きな場所へdragする操作。専用rowがcanvasの縦領域を使う |
| Retired effects | search自由drag、固定360px退避、旧absolute translateと競合CSSだけ。editor/reviewの意味は退役しない |
| Discoverability/accessibility / feedback | 安定した順序、keyboard、wrap/scroll、focus維持。short高さでoverflow clipによる隠れなし |
| Canonical state update | query等は既存search owner。geometryはCSSのみ、新座標ref/observerなし |
| Undo/Redo | chromeはartifact Historyに入れず、図の履歴を維持 |
| Session / regeneration | chrome位置は新たに保存しない。既存Session/Result/draftを維持 |
| Export/artifact | HTML chromeは図/SVG/PNG/PDF/保存Resultへ入らない |
| Validation/error | CSSで専用row/workspaceを分離し、短い高さで到達性検査 |
| Failure/recovery / next action | drawer/resizeでquery/active/focus保持。Result消失は既存visibilityで閉じる |
| Scientific-output | chromeのみ、生物学的値/比較/scale不変 |
| Cache/provenance | chromeをWorker/cache/requestに入れない |
| Performance | CSSのみ。位置computed/global drag listenerを削除 |
| Compatibility | 保存形式維持、自由drag退役を明示 |
| Architecture | index.htmlへgeometry ownerを収束しJS判断と旧CSS例外を削除。第二canvas/editorなし |
| Evidence available/missing | 10条件のbase矩形/source。新P01–P03/R01、hit target/keyboard/touch/short高さは未実装 |
| Residual risk | 専用rowがcanvas縦領域を減らし、検索を図の近くへ動かせなくなる。wrap/scroll/最低高さで実操作を確保 |
| Route | DURABLE_AUTHORITY_REQUIRED |

## 承認本文

以下が承認された回答の全文である。owner/date/rationale/preservation/retirement/risk を補完・拡張せず恒久契約へ転記する。

```text
PRODUCT_DECISION
Concern: web.preview-search-placement
Scenario revision: 1
Choice: A / DOCKED-SEARCH-AND-CONTROLS
Rationale: 検索とzoom/resetをdrawer開閉や画面幅にかかわらず操作できることを優先し、検索バーの自由移動に伴う衝突と場外配置をなくす。
Must preserve: searchの全field/query/regex/Prev/Next/Open/Enter、active match/focus、全toolbar操作、drawer tab/Close/Escape、同じsearch/canvas/SVG/editor DOM、既存Session/Export/Historyを維持する。search/toolbarをworkspace外の専用rowへ置きdrawer被覆を防ぐ。通常高さ740px以上の受入条件でworkspace200px以上、short viewport/keyboard/200%zoomで全操作へscroll到達可能にする。JS位置判断と競合CSSを削除する。alignment reviewとeditor上下dockの仕様は変更しない。
May retire: 全幅での検索バー自由drag、固定360px退避、search/toolbarの旧absolute translateと競合CSSだけ。検索/編集機能、保存意味、editor/reviewの製品仕様は退役しない。
Accepted residual risk: 専用rowがcanvasの縦領域を減らし、検索を図の近くへ動かせなくなる。wrap/scroll、canvas最低高さ、実hit targetとkeyboard/touch検証で実操作を確保する。
Owner: satoshikawato
Decision date: 2026-09-26
```

承認本文 UTF-8（末尾 newline を含めない）の SHA-256: `4a215b8ffb32217529703db8c3ae83d0ea104f69b2fe32e6d5153af791e05f63`。

## レビュー用の機械表現

これは上記本文のフィールド対応を見せる非実行 JSON であり、checker が読む新 schema/registry ではない。
unmapped concern の実行可能な decision store は追加しない。active authority は既存規約の記録先へ S00 で記録し、dev merge を確認する。

```json
{
  "concern": "web.preview-search-placement",
  "scenarioRevision": 1,
  "choice": "A / DOCKED-SEARCH-AND-CONTROLS",
  "rationale": "検索とzoom/resetをdrawer開閉や画面幅にかかわらず操作できることを優先し、検索バーの自由移動に伴う衝突と場外配置をなくす。",
  "mustPreserve": "searchの全field/query/regex/Prev/Next/Open/Enter、active match/focus、全toolbar操作、drawer tab/Close/Escape、同じsearch/canvas/SVG/editor DOM、既存Session/Export/Historyを維持する。search/toolbarをworkspace外の専用rowへ置きdrawer被覆を防ぐ。通常高さ740px以上の受入条件でworkspace200px以上、short viewport/keyboard/200%zoomで全操作へscroll到達可能にする。JS位置判断と競合CSSを削除する。alignment reviewとeditor上下dockの仕様は変更しない。",
  "mayRetire": "全幅での検索バー自由drag、固定360px退避、search/toolbarの旧absolute translateと競合CSSだけ。検索/編集機能、保存意味、editor/reviewの製品仕様は退役しない。",
  "acceptedResidualRisk": "専用rowがcanvasの縦領域を減らし、検索を図の近くへ動かせなくなる。wrap/scroll、canvas最低高さ、実hit targetとkeyboard/touch検証で実操作を確保する。",
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
