# INSTRUCTION PROMPT — S04 — 既存Color ruleの拒否された入力を保持

あなたは既存Color rule pattern fieldの未確定入力と回復操作の担当です。
承認されたKEEP_REJECTED_PATTERN_DRAFTだけを実装し、Python評価とcanonical commitのownerを維持してください。

## 必ず当該ブランチを取得する

SESSION_CODEはs04です。
[総合計画書](../MASTER_PLAN.md)の共通取得手順で**origin/fix/issue-601-bug15-bug19を取得し、
専用cloneのfix/issue-601-bug15-bug19で作業してください**。
共有treeや他sessionのbranch/index/環境を使わず、同じremote branchのwriterを一つにします。
AGENTS.md、CLAUDE.md、Web CLAUDE、総合計画、
[診断公開の承認](../decisions/DECISION_01_ERROR_DISCLOSURE.md)、
[field回復の承認](../decisions/DECISION_02_REGEX_EDIT_RECOVERY.md)、先行SESSION_XX_RESULT.mdを読んでください。

## 所有範囲

- web/js/app/feature-editor/rule-actions.jsと既存rule composition。
- web/index.htmlの既存Color rule pattern field、必要なfocused helper。
- Session/reset/History/row replacement/mode/drawerの既存ownerへの限定transition呼出配線。
- 関連Node/Playwright tests、SESSION_04_RESULT.md。
- Label TSV、新規rule、preset、whitelist、Search、Session schemaへdraft機構を拡張しない。

## 作業

1. 対象fieldの表示textとaccepted canonical ruleを明確に分ける。
   target row identity/revision/current documentに結び付く最小transient stateを既存rule ownerに持たせる。
   新しいpersisted rule IDや汎用draft manager、canonicalの複製ownerを作らない。
2. 通常changeで既存preparationの一回評価を使い、valid current responseだけ既存atomic live commitへ通す。
   pendingはvalid/non-match/Appliedと扱わない。JS事前compile・validate-only・keystroke連続callは追加しない。
3. failureでtextを同じfieldに残し、S02 normalizerのsyntax/runtime区別とfield errorを関連付ける。
   Not applied、Save/Generateはlast accepted rule、Exportは現在Resultを使うことを示す。
   位置が不明なら架空の位置を表示しない。
4. Retryは同じdraftの一回評価、Revertは現在accepted値へ戻す明示owner transitionとする。
   field draft/failure/Retry/Revertはartifact Historyを増やさず、成功editだけ従来の一entry。
5. drawer close/reopenと同documentの一時mode切替はdraftを保持し、旧pending responseを無効化する。
   draftを保持したまま不可視modeへcanonical mutationをcommitしない。
6. row削除、対象ruleのUndo/Redo置換、成功document/session置換、resetでownerがdraftを解放する。
   Session置換が失敗したら旧sourceを復元してdraftも維持する。
   unrelated Historyや単なるdrawer closeをdraft削除とみなさない。
7. 古いresponseはrow/revision/catalog/Result/mode/documentの現在性で除外する。
   keyの取り違え、row reorder後の別行commit、Session load後の旧commitを防ぐ。
   watcherに正しさを依存せず、bulk ownerが明示transitionを呼ぶ。
8. canonical Session保存/Generate projection/Export/UndoRedoの意味を維持し、draftは非永続とする。
   元patternをerror details/Copy/consoleへ自動公開しない。
9. 無条件input復帰のfinallyを置換し、同じfieldのdraft更新をtemplateとownerで二重に所有しない。
   長くなる単目的helperは既存feature-editor配下へ分け、単一entry ownerを維持する。

## 必須検証

不正[、Python専用有効式、runtime initialization failure、empty/unrelated catalogで
display text、accepted rule、Result、Historyを別々に確認する。
修正/Retry成功のatomic commitと一History entry、Revertの無mutationをassertする。

遅い旧edit、新edit、row remove/reorder、対象rule Undo/Redo、unrelated History、
成功/失敗Session置換、reset、drawer close/reopen、一時mode切替を検査する。
Save/fresh Load/Generate/Exportでdraft非永続とaccepted/current Resultの一致を確認する。
desktop/390 px、keyboard、errorのfield関連付け、focus、操作到達性をbrowserで確認する。
25,000-feature prepared reuseとWorker呼出数に余分な評価がないことを確認する。

SESSION_04_RESULT.mdにdraftの全transition、canonical境界、race/History/Sessionの結果、
owner/path削除、S05で残る受入を保存する。

English commit title: Keep rejected color-rule edits available for correction
Summary: Add focused transient draft recovery while preserving Python validation and canonical state.

## セッション終了時のコミット・プッシュ

総合計画書の共通終了手順に従い、SESSION_04_RESULT.mdと対象変更を
**検証後にコミットし、当該同名remote branchへプッシュしてください**。
通常のtargetはfix/issue-601-bug15-bug19です。
branch/upstream、staged scope、remote実状態を確認し、force-pushやmain/dev直接pushは行いません。
remote/local SHA一致、result file、完了/未完了の開始条件を次sessionへ引き継ぎます。
