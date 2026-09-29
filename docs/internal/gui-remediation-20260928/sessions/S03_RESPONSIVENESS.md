# S03 INSTRUCTION PROMPT: 操作応答、Generate受付、二重計算の退行防止

gbdraw Webの比較切替・入力・Generate受付を軽くし、不要なPending表示と表示専用計算を除去してください。同じ重複計算が再導入されたときに失敗する検証も実装します。

## 前提と場所

**fix/gui-feedback-remediation-20260928** と /mnt/c/Users/genom/GitHub/gbdraw/.worktrees/gui-feedback-remediation-20260928を継続使用。
[総合計画](../00_MASTER_PLAN.md)、[SESSION_LOG](../SESSION_LOG.md)、results/S00.md/S02.md、AGENTS/CLAUDE/Web CLAUDEを読む。
S02がorigin/devへ反映済みであることを確認する。S00の固定runner/入力/予算を用いる。

## 変更責任

app/app-setup.js、app/run-analysis.js、app/generation-status.js、services/config.js、services/session-request.js、app/feature-editor/label-override-table.jsと対象tests。
index.htmlはResult/Generate上の指定Pendingブロックのみ担当。検索・Editor・tooltip・Align/ReviewのmarkupはS06担当。

1. Product Contract PD-OI-037と受入記述を、指定された常時Pending説明の削除へ同時更新する。static Contractのlifecycle文も改正済み同時更新ルールに合わせる。既存の明示要求を再質問せず、不足する正式receipt事項だけ扱う。
2. consumerを追跡し、generationApplicationFeedback/getGenerationApplicationStatusと表示専用のintent baseline/bookkeepingを除去する。非表示化だけで終了しない。
3. canonical request builder、compareCanonicalRenderRequests/assertCanonicalRenderRequestsEquivalent、committedCanonicalSession、artifact capture/restore、History/Save/Exportの意味を維持する。
4. buildLabelOverrideRowsはlabel/visibility overrideが空ならmetadata/index構築前に正しい空結果を返す。visibilityだけある場合などの負例も検証する。
5. Generateの既存operationが受付→前処理→検証→render→settlementを所有するよう整理する。別の並行busy ownerを作らず、開始時に可視状態を公開してpaint機会を与える。
6. 処理中・取消・実エラー・live edit失敗の通知は既存operation/error refsから提供する。重いcomputedをsr-only用に残さない。
7. CW-01〜06の初期自動検証を生成label表とselector metadata/index、指定UI操作に限定し、既存runtime-test-hooksとcontract/performance testsへ実装する。計測coverageとpositive controlがなければ0仕事と判定しない。既存test/helperがbase mapのhard contractかを先に確認し、mapped coverageを候補変更だけへ置き換えない。

## テスト

- 空override、labelのみ、visibilityのみ、bulk、非空編集後のGenerate/Save/Load。
- 比較切替はselection/historyを正しく更新し、不要なPython/LOSAT dispatch・全feature表構築を起こさない。
- Generate直前のinput blur/Historyを含めて即時受付を計測。
- 遅延helperを用い、二重クリック1操作、前処理Cancel、遅着結果拒否、validation review、error、retryでbusy/Result/Historyが正しいこと。
- 空早期終了の除去、status全feature走査の再導入を一時変異としてテストが失敗すること。変異は復元する。
- 実在する複数consumerが残る場合だけ同一構築結果の共有を検証する。残らなければ実在経路への重複構築変異で感度を確認し、仮想consumerを作らない。独立したbefore/after観測、新Generate、retryの正当な再計算は通す。
- generation-feedback等のテストは表示assertionだけを更新し、Result/draft、applied authority、Save/Export、Undo/Redoのbehavior coverageを維持する。

性能はS00と同じ条件で測り、可視反映とsettlementと完了時間を分ける。初期対応ではWeakMap、toRaw一括適用、新Worker、新cache、Vue build変更をしない。残ったボトルネックが実測された場合だけ範囲を再検討する。

## 終了

既存generation/session/history unit testsと対象Playwrightを実行し、必要なfast/architecture gateを確認する。変更したruntime、tests、契約を別々にレビューする。

results/S03.mdにCWごとの正例/負例、変異感度、性能raw data、consumer削除境界、残余リスク、browser configとtest選択を含む実行commandを記録する。SESSION_LOGを更新して担当差分をcommitし、英語title/summaryを残す。外部公開は承認範囲に従う。
