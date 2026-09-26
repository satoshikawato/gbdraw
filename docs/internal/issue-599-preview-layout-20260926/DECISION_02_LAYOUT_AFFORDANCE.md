# 承認済み Decision Pack — Layout edit の発見

Concern: `web.layout-edit-affordance`。Scenario revision: `1`。
選択: `A / EXPLICIT-MODE-WITH-DISCOVERABLE-TARGETS`。Product Decision Owner: `satoshikawato`。
承認日: 2026-09-26。状態: **製品結果は明示承認済み、恒久的 base authority への記録は S00 の作業**。

## 承認の範囲と根拠

対象は [Issue #599](https://github.com/satoshikawato/gbdraw/issues/599) のこの concern のみ。
Product Decision Owner は、装飾配置の継承、Layout edit の発見、検索・toolbar の配置という3つの独立した推奨 A に対して「すべて推奨案で承認します。」と明示回答した。
この Pack は、そのうち `web.layout-edit-affordance` の完全な A を記録する。他 concern の承認記録は [一覧](00_APPROVED_PRODUCT_DECISIONS.md) から参照する。選択 B は記録しない。
電子的な明示承認として記録し、架空の手書き署名は作らない。

調査基準: `origin/dev@d457b7189b137185a8dec800819a312c30b969fa`。
この基準では active BD は存在せず、OIPC revision 16 が intent/failure 等を規定するが、この完全な outcome は未登録だった。
developer preflight の Product 判断は上記承認で解決済み。恒久記録の経路は **DURABLE_AUTHORITY_REQUIRED**。
S00 は最新 dev の map/契約/accepted decisions を再調査し、既存 authority と重複しない記録先と番号を確定する。
ここに仮 BD/PD-OI 番号や runtime 自己承認は置かない。

## 選択した完全な製品結果

| Dimension | 承認済み A の outcome |
| --- | --- |
| Complete normative outcome | OFF の drag は従来の canvas pan。supported target は help cursor/hover枠と「Turn on Layout edit to move this item」を表示。toolbar の常設説明、keyboard focus、touch でも同じ情報へ到達。ON は grab、drag 中は grabbing。target 上で mode を自動ONにしない |
| Preserved effects | feature click、editable label、legend個別編集、Shift/Ctrl、background pan、record/alignment、History/Session/export |
| Added effects | OFF でも対象と有効化方法を発見。hover に閉じない説明 |
| Lost effects | 機能の喪失なし。有効化の1手順は残る |
| Retired effects | OFF/ON の意味を区別できないcursorと説明不足だけ |
| Discoverability/accessibility / feedback | 常設説明、toggle aria-pressed/説明、focus/touch。大量の SVG tab stop は追加しない |
| Canonical state update | 既存mode ref。hintは派生表示でSVG/History/canonicalを変更しない |
| Undo/Redo | hover/hint は履歴なし。実 drag のみ既存1操作 |
| Session / regeneration | mode/Result復元後に表示をrebind。hintは保存しない。Generate継承はPack01が決める |
| Export/artifact | cursor/outline/hint は Preview 専用。plain/interactive SVG、PNG/PDF、保存Resultへ入れない |
| Validation/error | 既存 composition eligibility で対象限定。未対応targetを動かせると説明しない |
| Failure/recovery / next action | Result/load/Historyの既存bind。hintが出せない場合も常設説明とtoggleを使用可能 |
| Scientific-output | 発見方法だけ。生物学的意味/comparison/alignment不変 |
| Cache/provenance | hintをrequest/cache keyに含めず、Preview transientをclean serializationで除去 |
| Performance | 既存bindで対象限定、hoverで全走査/Worker/History cloneなし |
| Compatibility | Session/modeの意味、writer/readerを維持 |
| Architecture | 既存target ownerとHTML/CSSの派生表示。別mode ownerなし |
| Evidence available/missing | base cursor/source。新hint、keyboard/touch、serialization の after は未実装 |
| Residual risk | 有効化の1手順が残る。常設説明とkeyboard/touch検証で発見を補う |
| Route | DURABLE_AUTHORITY_REQUIRED |

## 承認本文

以下が承認された回答の全文である。owner/date/rationale/preservation/retirement/risk を補完・拡張せず恒久契約へ転記する。

```text
PRODUCT_DECISION
Concern: web.layout-edit-affordance
Scenario revision: 1
Choice: A / EXPLICIT-MODE-WITH-DISCOVERABLE-TARGETS
Rationale: canvas panと配置編集の区別を保ちながら、対象と有効化方法を初めての利用者にも示す。
Must preserve: OFFの従来panと明示toggle、ONのtarget drag、supported targetへのhelp/hover説明とtoolbarの常設説明、keyboard focus/touchで同じ説明へ到達すること、feature/label/legend個別編集とShift/Ctrlの優先順位、record/alignmentの既存動作、実dragの1 History操作、Session/Exportを維持する。hintだけでcanonical値を変えず、成果物へhintを保存しない。modeを自動ONにしない。
May retire: OFF/ONの意味を区別できないcursorと説明不足だけ。gesture/編集機能/保存意味は退役しない。
Accepted residual risk: 移動前に有効化の1手順が残る。常設説明とkeyboard/touch検証で発見可能性を補う。
Owner: satoshikawato
Decision date: 2026-09-26
```

承認本文 UTF-8（末尾 newline を含めない）の SHA-256: `0365e0057b61e6b606c6714dc17f0e5c89d755ce651c1a8eec66d0139e0c65ec`。

## レビュー用の機械表現

これは上記本文のフィールド対応を見せる非実行 JSON であり、checker が読む新 schema/registry ではない。
unmapped concern の実行可能な decision store は追加しない。active authority は既存規約の記録先へ S00 で記録し、dev merge を確認する。

```json
{
  "concern": "web.layout-edit-affordance",
  "scenarioRevision": 1,
  "choice": "A / EXPLICIT-MODE-WITH-DISCOVERABLE-TARGETS",
  "rationale": "canvas panと配置編集の区別を保ちながら、対象と有効化方法を初めての利用者にも示す。",
  "mustPreserve": "OFFの従来panと明示toggle、ONのtarget drag、supported targetへのhelp/hover説明とtoolbarの常設説明、keyboard focus/touchで同じ説明へ到達すること、feature/label/legend個別編集とShift/Ctrlの優先順位、record/alignmentの既存動作、実dragの1 History操作、Session/Exportを維持する。hintだけでcanonical値を変えず、成果物へhintを保存しない。modeを自動ONにしない。",
  "mayRetire": "OFF/ONの意味を区別できないcursorと説明不足だけ。gesture/編集機能/保存意味は退役しない。",
  "acceptedResidualRisk": "移動前に有効化の1手順が残る。常設説明とkeyboard/touch検証で発見可能性を補う。",
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
