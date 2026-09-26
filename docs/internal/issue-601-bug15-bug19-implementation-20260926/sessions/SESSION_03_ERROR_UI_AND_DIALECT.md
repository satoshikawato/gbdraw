# INSTRUCTION PROMPT — S03 — 原因の返却・回復案内・Details/Copy・regex方言

あなたはoperation outcomeと利用者向けエラー表示の担当です。
構造化原因を失わずにGenerate/Align/import/exportの次操作へつなぎ、承認済みの限定Details/Copyと正しいregex方言を表示してください。

## 必ず当該ブランチを取得する

SESSION_CODEはs03です。
[総合計画書](../MASTER_PLAN.md)の共通取得手順で**origin/fix/issue-601-bug15-bug19を取得し、
専用cloneのfix/issue-601-bug15-bug19で作業してください**。
共有treeや他sessionのbranch/index/環境を使わず、同じremote branchのwriterを一つにします。
AGENTS.md、CLAUDE.md、Web CLAUDE、総合計画、
[診断公開の承認](../decisions/DECISION_01_ERROR_DISCLOSURE.md)、
[field回復の承認](../decisions/DECISION_02_REGEX_EDIT_RECOVERY.md)、先行SESSION_XX_RESULT.mdを読んでください。

## 所有範囲

- web/js/app/run-analysis.js、similarity-alignment.jsのoutcomeと既存rollback/retry。
- rule/color/label import、export caller、app-setup.jsの通知/composition配線。
- web/index.htmlのoperation alert、Details/Copy、方言説明。
- web/js/app/feature-search/search-core.js、services/standalone-interactivity-assets.jsの方言案内。
- 関連tests、SESSION_03_RESULT.md。
- PDF font/asset、regex evaluator、field draft寿命は変更しない。

## 作業

1. S02のcontractとfailure一覧を受け取り、run-analysisの全status:error returnを監査する。
   callerが原因を必要とする全pathは{status:'error',error}を返す。
   generatedArtifactCandidateの除去でerrorを捨てず、global errorLogから原因を推測しない。
2. 実run-analysisのengine-error/catchをAlign Applyまで接続する。
   原因を同じreviewへ表示してdraftとApply retryを残し、Result/request/orientation/Historyを保護する。
   errorなし応答は原因情報欠落として案内し、入力ミスやnetwork failureと断定しない。
3. cancel/stale/supersededをerror UIに変換せず、旧operationで新しいalertを上書きしない。
   初回・rollback成功・rollback failureの案内はtransactionの実状態から導く。
4. 全移行callerを同じnormalizerへ収束する。準備失敗をInvalid ruleと呼ばず、
   ExportをGeneration Errorと名付けない。source側で個別prefix/classifierを増やさない。
5. Detailsを初期collapseの任意表示にし、bounded安全情報だけを描画する。
   Copy diagnosticsは表示中の安全modelだけを手動コピーし、隠れたexception fieldを再取得しない。
   Clipboard不可時もDetailsの選択コピーと普通のrecoveryを維持し、元failureを置換しない。
6. keyboard/focus/role/name、copy成否通知、390 pxの到達性を確認する。
   retry/edit/saveは各既存ownerの実在するactionへつなぐ。汎用dispatcherは作らない。
7. Color/LabelにPython regex・case-insensitiveと短い(?i)/(?P<name>)の例を案内する。
   実ownerがPythonでない入口へ誤った説明を付けない。
8. Feature SearchとstandaloneにはJavaScript regex・case-insensitive、
   Invalid JavaScript regular expressionと単語検索へ戻れる既存操作を示す。
   JS→Python翻訳、runtime埋込み、検索target変更をしない。
9. 移行したraw summary/alert prefix/Circular特別extractorを除去し、known修正情報を保つ。
   S04用のfield error表示はS02 modelを再利用できるようにするが、別draft管理は追加しない。

## 必須検証

実orchestrationのfailure→caller→UIでknown causeと次actionをassertする。
stubbed runAnalysisがerrorを返すだけの既存unitを、実経路failureの代用にしない。
Alignでcause保持、draft/Result/orientation/History/request保持、retry成功を確認する。

synthetic sentinelを原exception/cause/stdout/cleanupに入れ、summary/Details/Copy/consoleへ出ないこと、
上限と再正規化、cancel/staleの非error、初回/旧Resultの真実な表示を検査する。
Circular/Linear、desktop/390 px、keyboard/focus、Details/Copy/clipboard不可をbrowserで確認する。
Python Color/Labelのmanual/TSV/preset/History/Session/Generate parityと
JS Search/standalone download bytesの方言・targetを保護する。原監査未特定は限界として残す。

SESSION_03_RESULT.mdへcaller移行表、全error return結果、actual integration evidence、
公開情報/Copy allowlist、UX/復旧、S04へのfield model契約を保存する。

English commit title: Show actionable operation errors and clarify regex dialects
Summary: Retain caller causes, safe diagnostics, recovery actions, and existing search semantics.

## セッション終了時のコミット・プッシュ

総合計画書の共通終了手順に従い、SESSION_03_RESULT.mdと対象変更を
**検証後にコミットし、当該同名remote branchへプッシュしてください**。
通常のtargetはfix/issue-601-bug15-bug19です。
branch/upstream、staged scope、remote実状態を確認し、force-pushやmain/dev直接pushは行いません。
remote/local SHA一致、result file、完了/未完了の開始条件を次sessionへ引き継ぎます。
