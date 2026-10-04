# Issue #602 署名済みProduct Decisions

状態: **SIGNED**。Product Decision Owner: **satoshikawato**。署名日: **2026-09-26**。

5つのConcernを独立して採択した記録である。01はB、02〜05はA。各回答のRationale、Must preserve、May retire、Accepted residual risk、Owner、Decision dateは採択された本文をそのまま保存する。

本書は人間の決定記録とserializationのレビュー資料であり、active Product authorityでも新しいdecision storeでもない。正式authorityは既存 `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` にS00で反映し、依存runtimeはそれがorigin/devへmergeされたbaseから実装する。mapped concernsへ変更されたbaseでは既存のmap/BD routeに従う。署名はpush、PR、merge、deploy、tag等の外部操作を許可しない。

[総合計画](00_MASTER_PLAN.md) / [実装セッション一覧](INSTRUCTION_PROMPTS/README.md)

## 採択一覧

| Pack | Concern | Revision | Signed choice |
| --- | --- | --- | --- |
| [01](01_METADATA_VISIBILITY_DECISION_PACK.md) | `linear.record-label-auto-visibility` | 2 | B / AUTO-FRESH-RESET-WITH-DISCLOSURE |
| [02](02_DEFINITION_DEFAULT_DECISION_PACK.md) | `linear.definition-display` | 2 | A / D1-B-WEB-LOCKED-FRESH; retain D2-P and D3-A |
| [03](03_EDIT_FEEDBACK_DECISION_PACK.md) | `web.edit-application-feedback` | 1 | A / DERIVED-APPLICATION-STATUS |
| [04](04_COMPACT_EDITOR_DECISION_PACK.md) | `web.editor.compact-presentation` | 1 | A / DOCKED-COMPACT-EDITOR |
| [05](05_COMPACT_ALIGNMENT_REVIEW_DECISION_PACK.md) | `web.similarity-alignment.review-presentation` | 1 | A / DOCKED-COMPACT-ALIGNMENT-REVIEW |

## 署名済み回答全文

### 01 — linear.record-label-auto-visibility

```text
PRODUCT_DECISION
Concern: linear.record-label-auto-visibility
Scenario revision: 2
Choice: B / AUTO-FRESH-RESET-WITH-DISCLOSURE
Rationale: 共有行の図では簡潔な既定表示を維持し、情報が非表示になる理由とShowへの変更先を配置操作の場所で明示する。
Must preserve: fresh/resetの独立Auto、Show/Hideの明示値、diagram-wide Auto解決、休眠行除外、既存Sessionと保存Result、Undo/Redo、GenerateとExportの区別。Auto非表示時はLayoutにも理由・対象field・図全体の範囲・次回Generateの効果・Record Labelsへの変更先を表示する。
May retire: なし。
Accepted residual risk: 共有行でAccession/Lengthが非表示になる結果自体は残る。説明を見落とす可能性があるため、layout操作場所とLabelsの両方で実効値と変更先を示す。
Owner: satoshikawato
Decision date: 2026-09-26
```

署名: **satoshikawato** — **2026-09-26**。上記回答全文を承認。

### 02 — linear.definition-display

```text
PRODUCT_DECISION
Concern: linear.definition-display
Scenario revision: 2
Choice: A / D1-B-WEB-LOCKED-FRESH; retain D2-P and D3-A
Rationale: Webで新しく作るLinear比較図ではDefinitionを共通左列にそろえ、行のalignやoffset後も名前を比較しやすくする。
Must preserve: Lock=trueの共通左端とconfigured gap、明示Lock=falseの共通幅中央とrow追従、単一/共有/混在行、既存text_anchorの受入範囲、保存Sessionの明示値と対応済み旧省略意味、読込時の保存Result、CLI/Python省略default。PD-OI-024のD2-Pの保存/手入力Subtitleと継承・ラベル区別、およびD3-AのReplicon/Organelle選択順・独立制御・既定falseをすべて維持する。Linear LayoutでON/OFFの違いとGenerate適用を常時説明する。
May retire: D1-Aのうち、Web fresh/resetがLock=falseを初期値として選ぶ部分だけ。OFFの明示操作とCLI/Python既定は退役しない。
Accepted residual risk: 新しいWeb図のDefinition外観が従来のfresh図と変わり、Definitionがrowに追従しなくなる。利用者はOFFを選べ、既存Sessionの値と保存Resultは勝手に変更しない。
Owner: satoshikawato
Decision date: 2026-09-26
```

署名: **satoshikawato** — **2026-09-26**。上記回答全文を承認。

### 03 — web.edit-application-feedback

```text
PRODUCT_DECISION
Concern: web.edit-application-feedback
Scenario revision: 1
Choice: A / DERIVED-APPLICATION-STATUS
Rationale: 利用者が編集後の図と次回Generateの変更を区別できるよう、操作の適用タイミングと生成設定の未適用状態を事実から表示する。
Must preserve: 操作単位のLive edit、Applies on Generate、Apply requiredを区別し、左Palette Instant Previewと右Alignment reviewを例外なく正しく分類する。canonical即時commitと必要時自動rerender、reviewのlocal draft、target-only操作とpending設定の分離、対応済みoverride継承、atomic Generate、失敗/Cancel/stale時の旧ResultとHistory、Undo/Redo、SessionのResult/draft分離、Exportの現在Result出力を維持する。Pendingとlive applying/errorを独立に示し、invalid/unknownをAppliedとしない。Generateによる配置再計算とzoom reset、Save/Exportの意味を事前に説明する。
May retire: 製品の適用タイミング・保存・復旧・編集機能は退役しない。全DOM座標や手動位置を再生成後にも無条件に保持する保証は新設しない。
Accepted residual risk: 生成intent比較の漏れや比較基準の誤更新は誤表示を生みうる。生成・live commit・履歴・Sessionの一致テストを必須とし、根拠不足はunknownとして表示する。Status目的のWorker呼出、genome byte読取り/hash、SVG/checkpoint cloneは受け入れない。
Owner: satoshikawato
Decision date: 2026-09-26
```

署名: **satoshikawato** — **2026-09-26**。上記回答全文を承認。

### 04 — web.editor.compact-presentation

```text
PRODUCT_DECISION
Concern: web.editor.compact-presentation
Scenario revision: 1
Choice: A / DOCKED-COMPACT-EDITOR
Rationale: 狭いPreviewでも即時編集の変化を図で確認できるよう、図とEditorを上下の領域へ配置する。
Must preserve: 同じSVGとEditor、全tabと同期可用性、canonical live commitと必要時rerender、既存History/Session/Export、camera操作、keyboard、Close/Escapeのvisibility-only意味、選択tab、Result置換/失敗復旧。390×844/740では利用可能幅全体かつ高さ200px以上のcanvasを確保し、Editor内容を独立scrollさせ、Close/headerとtoolbarを操作可能にする。短いviewport/soft keyboardでは全操作へscrollで到達できる。wideのside drawerを維持する。
May retire: 狭いPreviewでEditorが横から全面高さを覆う表示配置だけ。編集機能や保存意味は退役しない。
Accepted residual risk: 上下分割で図とEditor listの縦領域が短くなり、list scrollが増える。実操作のpointer/keyboard/browser検証を必須とし、複製Preview・SVG clone・第二editorによる回避は受け入れない。
Owner: satoshikawato
Decision date: 2026-09-26
```

署名: **satoshikawato** — **2026-09-26**。上記回答全文を承認。

### 05 — web.similarity-alignment.review-presentation

```text
PRODUCT_DECISION
Concern: web.similarity-alignment.review-presentation
Scenario revision: 1
Choice: A / DOCKED-COMPACT-ALIGNMENT-REVIEW
Rationale: 狭いPreviewでもalignment候補をcanvasで確認できるよう、reviewを図の下段に固定し、候補比較へ操作を集中させる。
Must preserve: PD-OI-031/034と現行transform/plan/reset/historyのすべての結果。resolvedの通常自動Apply、ambiguousと明示reviewのlocal draft、独立Select/Skip、候補根拠とreference identity、1つのMatch reference directionと各targetの結果方向、canvas操作、local編集でWorkerを呼ばないこと、Applyの共有Python batch validationとatomic Result/History。失敗時draft/error/retry、Cancel/stale/superseded時の以前のResult/orientation/History、Session/regeneration/Export、focus復帰を維持する。390×844/740では利用可能幅全体かつ高さ200px以上のcanvasを確保し、候補listをscroll、Apply/Cancelを到達可能にする。狭いreview開始時はEditorをownerで閉じ、tabを保持し、review中は理由付きでopenをdisable、終了後は明示reopen可能。wideのdragと非モーダルcanvasを維持する。
May retire: 狭いPreviewでreviewを自由にdragする操作、およびreview中にEditorを同時openする継続だけ。候補や方向の選択、Apply前draft、failure/retryは退役しない。
Accepted residual risk: 狭いreviewではlist scrollが増え、自由に位置を動かせなくなる。開始時Editorは閉じるがtabは保持し、終了後再openできる。位置変更で候補draftやResultを変えないことをbrowserで確認する。
Owner: satoshikawato
Decision date: 2026-09-26
```

署名: **satoshikawato** — **2026-09-26**。上記回答全文を承認。

## 機械表現

以下は署名済み回答のfieldを機械表現へ変換したもの。`Scenario revision` のみinteger、他のfieldは原文のstring。ID割当、authority merge SHA、runtime head、evidence referencesはProduct回答に追加していない。正式recordへの配置・supersession・Contract revisionはS00の作業であり、この配列を新しいCI registryとして読み込ませない。

```json
[
  {
    "concern": "linear.record-label-auto-visibility",
    "scenarioRevision": 2,
    "choice": "B / AUTO-FRESH-RESET-WITH-DISCLOSURE",
    "rationale": "共有行の図では簡潔な既定表示を維持し、情報が非表示になる理由とShowへの変更先を配置操作の場所で明示する。",
    "mustPreserve": "fresh/resetの独立Auto、Show/Hideの明示値、diagram-wide Auto解決、休眠行除外、既存Sessionと保存Result、Undo/Redo、GenerateとExportの区別。Auto非表示時はLayoutにも理由・対象field・図全体の範囲・次回Generateの効果・Record Labelsへの変更先を表示する。",
    "mayRetire": "なし。",
    "acceptedResidualRisk": "共有行でAccession/Lengthが非表示になる結果自体は残る。説明を見落とす可能性があるため、layout操作場所とLabelsの両方で実効値と変更先を示す。",
    "owner": "satoshikawato",
    "decisionDate": "2026-09-26"
  },
  {
    "concern": "linear.definition-display",
    "scenarioRevision": 2,
    "choice": "A / D1-B-WEB-LOCKED-FRESH; retain D2-P and D3-A",
    "rationale": "Webで新しく作るLinear比較図ではDefinitionを共通左列にそろえ、行のalignやoffset後も名前を比較しやすくする。",
    "mustPreserve": "Lock=trueの共通左端とconfigured gap、明示Lock=falseの共通幅中央とrow追従、単一/共有/混在行、既存text_anchorの受入範囲、保存Sessionの明示値と対応済み旧省略意味、読込時の保存Result、CLI/Python省略default。PD-OI-024のD2-Pの保存/手入力Subtitleと継承・ラベル区別、およびD3-AのReplicon/Organelle選択順・独立制御・既定falseをすべて維持する。Linear LayoutでON/OFFの違いとGenerate適用を常時説明する。",
    "mayRetire": "D1-Aのうち、Web fresh/resetがLock=falseを初期値として選ぶ部分だけ。OFFの明示操作とCLI/Python既定は退役しない。",
    "acceptedResidualRisk": "新しいWeb図のDefinition外観が従来のfresh図と変わり、Definitionがrowに追従しなくなる。利用者はOFFを選べ、既存Sessionの値と保存Resultは勝手に変更しない。",
    "owner": "satoshikawato",
    "decisionDate": "2026-09-26"
  },
  {
    "concern": "web.edit-application-feedback",
    "scenarioRevision": 1,
    "choice": "A / DERIVED-APPLICATION-STATUS",
    "rationale": "利用者が編集後の図と次回Generateの変更を区別できるよう、操作の適用タイミングと生成設定の未適用状態を事実から表示する。",
    "mustPreserve": "操作単位のLive edit、Applies on Generate、Apply requiredを区別し、左Palette Instant Previewと右Alignment reviewを例外なく正しく分類する。canonical即時commitと必要時自動rerender、reviewのlocal draft、target-only操作とpending設定の分離、対応済みoverride継承、atomic Generate、失敗/Cancel/stale時の旧ResultとHistory、Undo/Redo、SessionのResult/draft分離、Exportの現在Result出力を維持する。Pendingとlive applying/errorを独立に示し、invalid/unknownをAppliedとしない。Generateによる配置再計算とzoom reset、Save/Exportの意味を事前に説明する。",
    "mayRetire": "製品の適用タイミング・保存・復旧・編集機能は退役しない。全DOM座標や手動位置を再生成後にも無条件に保持する保証は新設しない。",
    "acceptedResidualRisk": "生成intent比較の漏れや比較基準の誤更新は誤表示を生みうる。生成・live commit・履歴・Sessionの一致テストを必須とし、根拠不足はunknownとして表示する。Status目的のWorker呼出、genome byte読取り/hash、SVG/checkpoint cloneは受け入れない。",
    "owner": "satoshikawato",
    "decisionDate": "2026-09-26"
  },
  {
    "concern": "web.editor.compact-presentation",
    "scenarioRevision": 1,
    "choice": "A / DOCKED-COMPACT-EDITOR",
    "rationale": "狭いPreviewでも即時編集の変化を図で確認できるよう、図とEditorを上下の領域へ配置する。",
    "mustPreserve": "同じSVGとEditor、全tabと同期可用性、canonical live commitと必要時rerender、既存History/Session/Export、camera操作、keyboard、Close/Escapeのvisibility-only意味、選択tab、Result置換/失敗復旧。390×844/740では利用可能幅全体かつ高さ200px以上のcanvasを確保し、Editor内容を独立scrollさせ、Close/headerとtoolbarを操作可能にする。短いviewport/soft keyboardでは全操作へscrollで到達できる。wideのside drawerを維持する。",
    "mayRetire": "狭いPreviewでEditorが横から全面高さを覆う表示配置だけ。編集機能や保存意味は退役しない。",
    "acceptedResidualRisk": "上下分割で図とEditor listの縦領域が短くなり、list scrollが増える。実操作のpointer/keyboard/browser検証を必須とし、複製Preview・SVG clone・第二editorによる回避は受け入れない。",
    "owner": "satoshikawato",
    "decisionDate": "2026-09-26"
  },
  {
    "concern": "web.similarity-alignment.review-presentation",
    "scenarioRevision": 1,
    "choice": "A / DOCKED-COMPACT-ALIGNMENT-REVIEW",
    "rationale": "狭いPreviewでもalignment候補をcanvasで確認できるよう、reviewを図の下段に固定し、候補比較へ操作を集中させる。",
    "mustPreserve": "PD-OI-031/034と現行transform/plan/reset/historyのすべての結果。resolvedの通常自動Apply、ambiguousと明示reviewのlocal draft、独立Select/Skip、候補根拠とreference identity、1つのMatch reference directionと各targetの結果方向、canvas操作、local編集でWorkerを呼ばないこと、Applyの共有Python batch validationとatomic Result/History。失敗時draft/error/retry、Cancel/stale/superseded時の以前のResult/orientation/History、Session/regeneration/Export、focus復帰を維持する。390×844/740では利用可能幅全体かつ高さ200px以上のcanvasを確保し、候補listをscroll、Apply/Cancelを到達可能にする。狭いreview開始時はEditorをownerで閉じ、tabを保持し、review中は理由付きでopenをdisable、終了後は明示reopen可能。wideのdragと非モーダルcanvasを維持する。",
    "mayRetire": "狭いPreviewでreviewを自由にdragする操作、およびreview中にEditorを同時openする継続だけ。候補や方向の選択、Apply前draft、failure/retryは退役しない。",
    "acceptedResidualRisk": "狭いreviewではlist scrollが増え、自由に位置を動かせなくなる。開始時Editorは閉じるがtabは保持し、終了後再openできる。位置変更で候補draftやResultを変えないことをbrowserで確認する。",
    "owner": "satoshikawato",
    "decisionDate": "2026-09-26"
  }
]
```

## Authority反映時の限定範囲

- 01-B: fresh/resetのAutoと現在のdiagram-wide解決を維持し、説明と変更先を正式な保証にする。01-Aのdefault Showは採択していない。
- 02-A: PD-OI-024のD1のうちWeb fresh/resetだけLock ON。D2-P/D3-Aと明示OFF、保存、省略、CLI/Pythonを維持する。
- 03-A: 新しい表示は既存の適用境界を観測する。全設定自動Generateや新たな保存形式を許可しない。
- 04-A: 同じEditor/SVGの上下配置。live編集の意味は維持する。
- 05-A: 同時canvas操作と可視領域を新しい正式結果にする。PD-OI-035の旧390px遮蔽許可はその範囲でsupersedeし、identity、keyboard/Skip、非描画候補、desktop、focus、overlay非保存の寄与は維持する。narrow dragとreview中Editor同時openの限定退役は署名本文どおり。

各レコードの意味を拡張する修正はこの署名から推測しない。最新baseとの実質的な競合がある場合は、対応するConcernだけを再評価する。選択済み結果の一般的な再承認は不要。
