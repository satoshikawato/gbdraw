# Session 05 — 狭いPreviewのEditorを下段へ配置する

## Missionと製品結果

署名済み `web.editor.compact-presentation` revision 1 / **A / DOCKED-COMPACT-EDITOR** は、狭いPreviewに図とEditorの非重複領域を作る。390×844のbaselineでは幅374pxのPreviewを334pxのdrawerが覆い、toolbar6操作のcenter hitを遮る。既存テストはその状態を許容しているため、見える図と操作可能性の新しい受入を追加する。

同じResult、Editor内容、visibility/tab owner、live commit、History/Session/Exportを維持する。wideのside drawerは維持。Close/Escapeはvisibilityだけを変更する。

開始条件: S04完了、04-Aがbase authorityにある。

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

- `index.html` のPreview/Editor responsive CSSとmarkup、`app/right-drawer.js` の既存action呼出、必要最小限のpresentation helper/wiring。
- `tests/web/right-drawer*.js` の実geometry/pointer/keyboard checks。
- 第二visibility/tab ref、mini-preview、複製SVG、mobile専用Editor、generic overlay manager、外部sheet dependencyは追加しない。review draft/ApplyはS06が担当。

## 実装

1. Preview container基準の既存40rem境界を使い、狭い場合は上段canvas、下段Editorにする。grid/flexで実際の領域を分け、見かけのz-index変更だけで解決したとしない。
2. same markup/contentを使い、Editorのlistだけを独立scroll。header/Closeとcanvas toolbarを操作可能にする。drawn SVGは不変で、view fit/cameraは既存presentation意味を利用する。
3. wideのside lane、toggle、tab可用性と同期reconcile、reset/restore/失敗復旧を維持する。responsive切替でHistory、Result、canonical overridesを変更しない。
4. CSS/dvh/safe-areaで可視高さを扱い、不足が実証された場合だけviewport計測を既存lifecycleへ加える。JSとCSSに別の狭い判定を作らない。
5. narrow open状態でのtoolbar遮蔽を許容する既存testは、選択済みの新受入へ修正する。toleranceやhit判定を弱めず、図と操作領域を実際に確保する。

## 検証と受入

M01–M03/M05。390×844/740はsticky header/Generate barを差し引き、canvasを利用可能幅全体・高さ200px以上にする。390×500、320px幅、844×390、200% zoom、入力focusとsoft-keyboardに相当するviewport縮小で全操作に到達しscrollで回復できることを示す。実機keyboardと自動resizeの差が残る場合は証拠の限界を明記し、未観測を観測済みとしない。

閉/openでbounding rect、elementFromPoint、actual zoom/pan、Close/Escape、tab移動を確認する。色/label/visibilityを一つ編集して図に反映されることも確認する。幅変更やopen/closeでSVG/canonical state/Historyが変わらないことを測る。Circular/Linearとdesktopの既存操作を含む。

```bash
node --test tests/web/right-drawer.test.mjs
npx playwright test tests/web/right-drawer.playwright.spec.js --workers=1
node tools/check-web-change-budget.mjs --base origin/dev
git diff --check
```

既存Gallery Sessionを実用的なbrowser fixtureに使い、最終desktop/narrow画像をreadable scaleで目視する。internal evidenceは基準SHA、viewport、Session、capture commandを記録する。public Gallery画像を変更する作業はこのsessionに含めない。

## 完了handoff

`docs/internal/issue-602-proposal-20260926/results/S05.md` にbase/head、branch/upstream、正式authority、変更対象、削除した旧処理、command/result、受入ID、visual observations、残る制約を記録する。未実行checkを合格扱いしない。最終回答は実現結果と確認済み範囲を説明し、Englishのproposed commit titleとsummaryを示す。次は `SESSION_06_COMPACT_ALIGNMENT_REVIEW.md`。開始条件と引き継ぐcommit/証拠を明記する。
