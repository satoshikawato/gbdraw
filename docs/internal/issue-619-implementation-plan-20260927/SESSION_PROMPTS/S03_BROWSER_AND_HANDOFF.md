# INSTRUCTION PROMPT S03 — 統合browser検証と最終handoffを完成させる

## 開始条件と隔離

[共通開始・終了手順](README.md)を実行する。remoteの`fix/issue-619-circular-track-measure-inputs`を取得し、対象branchの既存checkoutを前sessionから引き継ぐ。同branchのwriterは直列。cloneは不要で、worktreeは別作業との隔離が必要な場合だけ使う。他sessionのcheckout・server・browser・artifactへ干渉しない。dirty/unrelated changesをstage/revertしない。

[総合計画](../MASTER_PLAN.md)、AGENTS/CLAUDE/Web CLAUDE、architecture/Product/Web change policiesを読む。scalarの科学的意味、唯一のtyped request/Worker、draft/Result分離、Auto/invaliddraft、History/Session/privacyを守る。実装選択はSOLID/KISS/DRY/YAGNIに照らし、少数owner/pathとsuperseded pathの削除を優先する。

## 目的

数値＋unit controlsの実アプリ、保存復元、History、Generate失敗復旧、download/native replay、狭幅accessibilityを検証し、ユーザーdocumentationと最終reviewを完成させる。

## Ownership / 前提

所有: `tests/web/circular-track-measure-input.playwright.spec.js`、必要fixture/対応tests、既存public technical pages、S03結果。見つかったin-scope defectはそのownerで修正し、production/testを再reviewする。checker/workflow/reference/social previewを変更しない。

S00/S01/S02のremote commits、authority、commands/input/environment、未完項目を確認する。同じ実装への有効な証拠は再利用し、changed/failing partsのみ再検証する。

## 作業

1. browser specを既存functional config/discoveryへ接続する。自分のcheckout/専用port/serverを使い、server reuseで別branchのcodeを試さない。実際のsource/wheelを記録する。
2. tobacco Gallery Sessionの1440px/390pxで四scalarのnumeric/unit表示を確認。no-objecttext、keyboard/tab、help、label/error association、no-overlapを確認する。
3. typedpx/factor、bare/px/%、1.5＋selector、decimal/exponent、tiny/largeprecisionを実編集する。Auto選択→入力/Clear/panel remount/Reset/Loadが承認済み規則であることを確認する。
4. focus/blur/no-opと表示操作のscalar/request/Result/History/Worker不変を確認する。numeric/uniteditをUndo/Redo。IME/trim/invalidtext保持を確認する。
5. valid Generate、invalid Generateのrowerrorと旧Result/request維持、訂正後retryを確認する。value/unitはexact equality、same-value geometryは一致。pixel↔factorの意図したgeometry変更を確認する。
6. Save Sessionの実ファイルをLoad。draft≠committed、disabled/inactive、settings-only、CLI-originを代表fixtureで確認する。手でschema/hashを改変しない。
7. 実SVGをdownloadして現在Resultとgeometryを比較し、保存Sessionをnativeでreplayする。Linear smoke/sharedprojection regressionを行う。trackedreferenceを再生成しない。
8. 既存public technical ownerへnumeric/unit/%/Auto/Generate/Historyのexact semanticsを記載する。新publicpageを増やさない。controlsを説明するGallerycaptureが実際に影響を受ける場合だけ、該当skillを読みowner recipeで再生成して視覚確認する。minimal testfigureを公開showcaseに使わない。
9. production/test/docs/generated diffを独立にreviewし、重複owner/parser/formatter/直接bindingが残っていないことを確認する。全C619結果とrequired gatesをまとめる。

## Commands / 受入

masterのfocused Node checks、architecture-contracts、通常policy gate。newbrowser specは`GBDRAW_WEB_TEST_PORT=<own port> npx --no-install playwright test --config=playwright.functional.config.js --project=chromium --workers=1 tests/web/circular-track-measure-input.playwright.spec.js`。Node/Python両Playwrightを確認し、runner欠落時は等価checkを実施して限界を記録する。

native replayはexisting Session CLI/API pathを用い、保存file・command・outputsを記録する。wheelはcurrentbranchからprepareする。全C619-01〜10が必須。未実行・失敗・未署名を残して全完了と呼ばない。

PR/deploymentは別の公開境界。PR wordingが依頼された場合だけapplicable PR skill/language checkを行う。dev staging/release promotionは既存workflowに従い、第二pipelineを作らない。

next action: 最終handoff。追加sessionが必要な実failureならownerとrequiredcheckを明示する。

## 完了と公開

所有範囲のproduction/test/docs diffを別々にレビューする。`SESSION_RESULTS/S03.md`にcommands/exits、基準SHA・environment/input・authority、C619結果、未実行check、残るboundary、次sessionの前提を記録する。実測していないgate/browser checkをpassと呼ばない。

**共通終了手順に従い、成果と結果文書を対象branchへcommitし、同名remote branchへpushする。** commit/pushを省略して終了しない。push直前にremote state、branch/upstreamを確認し、force push/main/devへの直接pushをしない。別writerが先へ進んでいればその内容を確認する。authority boundaryが未成立でも独立文書・evidenceはcommit/pushし、依存runtimeは変更しない。自分が起動したprocessだけを停止し、対象checkoutを次sessionへ引き継ぐ。

handoffにはremote branchと実際のcommit SHA、English commit title/summary、通ったchecks、未成立の境界、次に使うpromptを含める。
