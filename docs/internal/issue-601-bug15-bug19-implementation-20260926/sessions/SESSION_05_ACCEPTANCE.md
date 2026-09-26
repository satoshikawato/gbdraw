# INSTRUCTION PROMPT — S05 — 全経路の統合受入・仕様更新・引継ぎ

あなたはBUG-15/19実装の統合受入担当です。
総合計画のE/R/U/P/SCI/G受入を実際の操作と出力で確認し、不足する修正をその原因ownerへ収束させてください。

## 必ず当該ブランチを取得する

SESSION_CODEはs05です。
[総合計画書](../MASTER_PLAN.md)の共通取得手順で**origin/fix/issue-601-bug15-bug19を取得し、
専用cloneのfix/issue-601-bug15-bug19で作業してください**。
共有treeや他sessionのbranch/index/環境を使わず、同じremote branchのwriterを一つにします。
AGENTS.md、CLAUDE.md、Web CLAUDE、総合計画、
[診断公開の承認](../decisions/DECISION_01_ERROR_DISCLOSURE.md)、
[field回復の承認](../decisions/DECISION_02_REGEX_EDIT_RECOVERY.md)、先行SESSION_XX_RESULT.mdを読んでください。

## 所有範囲

- integration/regression tests、本計画directoryの受入結果。
- docs/REFERENCE/web-app.md、input-formats-and-tsv-schemas.md等の既存公開owner。
- 受入失敗の原因ownerに対する限定production修正と対応tests。
- authority、detector、CI guard、vendor、reference outputsの無関係変更は行わない。

## 作業

1. S00–S04 result、二承認record、正式dev authority、最新branch headを読み、
   各保存・退役・残余リスク条件と受入IDのAND-of-OR coverageを確認する。
   同じoption IDやテスト件数だけで完了判定しない。
2. render/helper両経路→Worker/client→caller→UI/Copyで構造化causeを確認する。
   known validationの修正情報、unknown、cleanup、初回/rollback、cancel/staleを測る。
3. 実run-analysis→Align Applyのcause/recoveryを試し、旧Result/request/draft/orientation/Historyと
   corrected retryを確認する。stub成功だけでE05を完了としない。
4. Color/LabelのPython parity、TSV/preset/manual/empty catalog、History/Session/Generateを確認する。
   JS Feature Searchとstandalone download bytesの方言・検索targetが変わらないことを確認する。
5. field draftの全lifecycle、Save/Generate/Export境界、rejected correction/Retry/Revert、
   stale row/session/modeとfocusを検証する。draftを永続化するfallbackを追加しない。
6. Circular/Linear、desktop/390 px、keyboard、Details/Copyとclipboard不可を確認する。
   sentinelがsummary/details/copy/consoleへ出ないことを自動assertする。
7. 既存25,000-feature preparation/reuseとWorker construction/call数、
   有効target/SVG、native exception、Session/request compatibilityを確認する。
   変更のない証拠は再利用し、追加変更・失敗・懸念がある部分だけ再検証する。
8. 受入不足があれば原因ownerを修正する。新経路、silent fallback、
   known→unknown一括変換、test timeoutやthreshold緩和で通さない。
9. 公開説明を既存reference ownerで更新する。新manual/pageは作らない。
   procedural tutorialやGallery操作画像が実際に変わる場合だけ関連skillを適用し、
   再現command/sessionから最終artifactを確認する。
10. production/tests/docs/generatedを独立reviewし、owner/pathと旧経路除去、
    OE/PE/CB非増加、native/科学/性能/privilege/cycle、trusted-base policyを確認する。
    新blocking testは既存PR inventoryに入れ、local-only specを唯一のgateにしない。
11. 原監査build未特定、未検証browser、未実行supported-version matrixを正直に記録する。
    本BUG-15/19の実装完了をBUG-07/PDFやIssue全体のclosure/deploy成功と混同しない。

## 必須成果物

SESSION_05_RESULT.mdに受入IDごとのcommand/input/result、実測/lazy reuse、
receipt条件coverage、production/test/doc/generated review、残る限界、
開始/検証head SHAとtrusted base SHA、PR/CI/dev stagingの状態を保存する。
最終commit SHAはcommit後のhandoffで伝え、自己参照amendは行わない。

focused Node/Python/Chromium、追加回帰tests、Ruff、必要なread-only SVG comparison、
trusted-base Gateとarchitecture契約を通す。CIのGateとReviewを区別する。
完成した資料とコードをcommit/pushし、remote SHA一致を確認して実装完了を報告する。
PR作成・merge・promotion/deployは個別に許可された時だけ進める。

English commit title: Verify error recovery and regex editing across web workflows
Summary: Complete integrated acceptance and update existing behavior documentation.

## セッション終了時のコミット・プッシュ

総合計画書の共通終了手順に従い、SESSION_05_RESULT.mdと対象変更を
**検証後にコミットし、当該同名remote branchへプッシュしてください**。
通常のtargetはfix/issue-601-bug15-bug19です。
branch/upstream、staged scope、remote実状態を確認し、force-pushやmain/dev直接pushは行いません。
remote/local SHA一致、result file、完了/未完了の開始条件を次sessionへ引き継ぎます。
