# Session 04 — Pendingと操作ごとの適用タイミングを表示する

## Missionと製品結果

署名済み `web.edit-application-feedback` revision 1 / **A / DERIVED-APPLICATION-STATUS** を画面へ接続する。生成設定draftと現在Resultを区別し、操作がLive edit、Applies on Generate、Apply requiredのどれかを説明する。

左PaletteのInstant Previewは即時反映、右Alignment reviewはApply前draftである。パネルの左右を適用規則にしてはいけない。Pendingとlive applying/errorを同時表示可能にし、失敗・Cancel・staleでPendingを消さない。Generateは配置再計算とzoom resetを伴い、SaveはResult/draft保存、Exportは現在Result出力である。

開始条件: S03の共用projection/適用基準と検証が完了。

## 共通の作業規則

- `AGENTS.md`、`CLAUDE.md`、`gbdraw/web/CLAUDE.md`、`docs/internal/PRODUCT_IMPACT_RATCHET.md`、`docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md` を読む。総合計画は `docs/internal/issue-602-proposal-20260926/00_MASTER_PLAN.md`、署名全文と機械表現は同directoryの `06_SIGNED_PRODUCT_DECISIONS.md`。別worktreeの場合は計画checkoutを参照元に保持する。
- 最初にworking tree、branch/upstream、前段commitと証拠を確認し、originをfetchする。他者の変更をpreserveし、main/devへcommitしない。実行する指示書の範囲でlocal編集・検証・修正を完了する。
- 正式authorityがbaseへmerge済みか、署名されたconcern/revision/全寄与が一致するかを確認する。署名済み選択を再度選び直したり、一般的な再承認を求めたりしない。未merge authorityや候補文書はruntimeを自己承認しない。
- 既存owner/pathを再利用し、canonical request/Worker/Result admission/Historyの第二経路を作らない。純粋な処理と副作用、適用の意味と表示を分離する。重複処理は同じ変更で除去する。新schema/互換reader/依存/汎用frameworkは追加しない。
- browser確認前にNodeとPythonのPlaywrightを確認する。Nodeがなければ同等のPython確認を実行する。Chromium sandboxだけの失敗は同じ確認を適切なescalationで再実行する。port競合は空きportを使い、既存serviceを停止しない。
- focused checksと必要gatesを実行し、失敗を修正する。同一コード/入力/環境の検証済み証拠は再利用する。長いtestは最低30分を許容して増分monitorし、test所有の短いtimeoutは弱めない。
- production、tests、docs、generated artifactsの差分を分けてレビューする。参考SVGは通常read-only、生成wheelはgitignored、social previewはowner-maintained。新たなProduct結果やarchitecture例外が必要なら該当範囲だけを停止し、独立作業を完了する。
- 検証済み変更を一つのlogical commitとして扱い、branch/upstreamを確認してEnglishのproposed commit titleとsummaryを提示する。push/PR/merge/deploy/tagは個別に明示された許可の範囲だけ。PRを作成・編集する場合は指定PR skillを読む。

## Branchと開始確認

実装には計画commitを含む既存の `fix/issue-602-linear-live-edit-20260926` を必ず使用する。originをfetchし、S00の正式authorityがorigin/devへmerge済みであることを確認する。S01ではそのorigin/devをこの実装ブランチへmergeし、S02–S07は前段の変更と証拠を持つ同じブランチを継続する。計画commitを保持し、ブランチの再作成やresetによる置換を行わない。正しいbaseの祖先関係と前段handoffを確認する。upstreamとpush先は同名の `origin/fix/issue-602-linear-live-edit-20260926` のみとする。

## Ownership / 非対象

- `app/app-setup.js` のwiring、必要ならfocused `app/generation-status.js`、`index.html` の操作/Result/Generate表示。
- 既存Palette、feature/label/legend、review、Generate、History/restore ownersからoperation factsを受け取る。
- S03の正規化をUIに複製しない。第二dirty flag/History/Worker、全設定自動Generate、confirmation modalを追加しない。
- 対応済みoverride継承の実際のregressionが見つかった場合だけ、`candidate-render.js` 等の現在ownerでfocused correctionを行う。

## 実装

1. Result/Generate付近に有効生成設定のPending/invalid/unknown/未生成を説明する。常に「現在のResult」と「次の成功Generate」の意味を分ける。S03の比較結果から派生し、入力eventごとのflagsを立てない。
2. 各操作へ適用分類を置く。Scale/crop/slots等はGenerate、単一live editingとInstant PreviewはLive、reviewはApply required。広いパネルheadingだけで例外を隠さない。
3. live rerender中/失敗を既存ownerのstateから表示する。直接適用済み編集は失敗でも保持する契約を説明し、まだ反映できていないgeometryをAppliedとしない。S03基準更新を利用して、scale Pending+色Applied等の組合せを成立させる。
4. record-display/Paletteの局所説明は共用比較の担当部分へ接続するか、その範囲と明示する。globalに同じ意味を決めるownerを並立させない。Statusの表示actionはcanonical stateやHistoryを変更しない。
5. Generate前に配置再計算、zoom reset、Undo復帰を短く説明する。全手動positionを保持すると約束しない。現在対応する色/ラベル/visibility/record layoutの継承は確認する。
6. Save/Export付近でResult/draft保存と現在Result出力を説明する。操作をGenerate triggerにしない。role=statusは重要遷移だけpoliteに読み上げ、全keystrokeの反復を避ける。

## 検証と受入

E01–E05の表示とend-to-end継続。生成→scale Pending→live色/label→Save/Load→Export→Generate成功→Undo/Redoを実操作で確認する。Generate失敗/Cancel/stale、live rerender error、review Cancel/Apply、値を元へ戻すscenarioも含む。Export図にpending設定が混ざらず、LoadだけでResultが変わらないことを比較する。

```bash
node --test tests/web/session-draft-authority.test.mjs tests/web/palette-history-preservation.test.mjs tests/web/candidate-render.test.mjs tests/web/history.test.mjs
npx playwright test tests/web/palette-history-preservation.playwright.spec.js tests/web/history-generated-authority.playwright.spec.js --workers=1
node tools/check-web-change-budget.mjs --base origin/dev
git diff --check
```

必要な表示regressionは既存fixtureへ追加し、Circular/Linear両modeで確認する。Status表示のためにWorkerやfull serializationが増えていないことをS03のstructural evidenceとwiring diffで示す。

## 完了handoff

`docs/internal/issue-602-proposal-20260926/results/S04.md` にbase/head、branch/upstream、正式authority、変更対象、削除した旧処理、command/result、受入ID、visual observations、残る制約を記録する。未実行checkを合格扱いしない。最終回答は実現結果と確認済み範囲を説明し、Englishのproposed commit titleとsummaryを示す。次は `SESSION_05_COMPACT_EDITOR.md`。開始条件と引き継ぐcommit/証拠を明記する。
