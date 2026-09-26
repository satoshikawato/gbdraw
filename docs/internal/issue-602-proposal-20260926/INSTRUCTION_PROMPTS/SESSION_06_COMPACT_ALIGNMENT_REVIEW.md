# Session 06 — Alignment reviewを図の下段へ配置し、draftを維持する

## Missionと製品結果

署名済み `web.similarity-alignment.review-presentation` revision 1 / **A / DOCKED-COMPACT-ALIGNMENT-REVIEW** は、狭いPreviewでも候補とcanvasを同時に確認・操作できる結果を選択している。

既存PD-OI-035 scenario 2は390px遮蔽を許可していた。この旧mobile例外をS00で限定supersedeした正式authorityがbaseに必要である。resolvedの通常自動Apply、ambiguous/明示reviewのlocal draft、exact identity、Select/Skip、1つのMatch reference direction、Python validation、atomic Result/History、失敗retry、Cancel/stale保持はすべて維持する。

開始条件: S05完了、05-Aと旧390px例外の限定supersessionがbaseにある。

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

- `index.html` の同一review表示、`app/app-setup.js` のpalette位置/focus wiring、必要なら既存top-level create entry下のfocused presentation helper。
- `app/similarity-alignment.js` のlocal draft/Apply/Cancel/retry ownerは維持。`app/right-drawer.js` のclose actionを使う。
- 候補resolver、Python strand/eligibility、plan/transform/Session schema、LOSAT cache、別review controller、第二SVGは変更しない。

## 実装

1. 実装前にresolved明示review、複数target ambiguity、automatic render失敗からのretryを390pxでcaptureし、遮蔽率・Apply/Cancel到達・canvas pointerを記録する。Editorのbaselineをreviewの測定へ転用しない。
2. 狭いPreviewでは同じreviewを下段へdockし、上段canvasを確保する。candidate listは独立scroll、Apply/Cancelは到達可能。Teleport/DOM配置を切替える場合は同じcontroller/draftを保持し、候補やフォームを初期化しない。
3. narrowの自由dragを退役し、wideは既存dragとnon-modal paletteを維持する。位置/clamp処理をfocused presentation責務へ寄せるなら旧inline処理を同じ変更で除去する。幅変更でPython helper/render/Historyを起動しない。
4. review開始とwideからnarrowへのresizeではEditorをownerでcloseし、tabを保持する。review中はopen toggleを理由付きでdisable。review終了後は明示reopenできる。hidden open状態や自動復元用mirrorを作らない。
5. local Select/Skip/方向変更はWorkerを呼ばずcanvas markerと結果方向を更新。Applyは既存batch validation/atomic path。failed Apply/automatic renderはdraftとunderlying errorを保ちretryできる。source/Resultが変わるstale draftは既存ownerのguardで処理する。
6. 非モーダルrole=dialogとfocus/復帰を維持する。aria-modal/focus trapを加えずcanvasへkeyboardで到達できるようにする。narrowで候補marker、pan/zoom、toolbarが実際に操作できることを確認する。

## 検証と受入

M01–M05に加え、PD-OI-031/034と既存transform/plan/reset/historyの全継続。390×844/740でcanvas幅全体・高さ200px以上、短い高さ/landscape/zoom/keyboardで全操作に到達できることを確認する。

opened reviewのままresize、Select/Skip、Match reference direction、Apply成功/失敗/retry、Cancel、stale/superseded、Undo/Redo、Save/Load/regenerationを確認する。drawer close/tab保持/open disable/reopenも確認する。guides/numbers/draftがResult/download/Sessionへ混ざらないことを比較する。resolved normal Alignが不要なreviewを開かないことも保護する。

```bash
node --test tests/web/similarity-alignment-actions.test.mjs tests/web/right-drawer.test.mjs
npx playwright test tests/web/similarity-alignment-ui.playwright.spec.js tests/web/right-drawer.playwright.spec.js --workers=1
node --test tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs --base origin/dev
git diff --check
```

画像の目視とactual pointer/keyboard操作を両方記録する。仕様に合わせて旧遮蔽許容assertionを更新するが、候補選択/方向/復旧のassertionを消して品質を下げない。

## 完了handoff

`docs/internal/issue-602-proposal-20260926/results/S06.md` にbase/head、branch/upstream、正式authority、変更対象、削除した旧処理、command/result、受入ID、visual observations、残る制約を記録する。未実行checkを合格扱いしない。最終回答は実現結果と確認済み範囲を説明し、Englishのproposed commit titleとsummaryを示す。次は `SESSION_07_ACCEPTANCE_AND_DOCS.md`。開始条件と引き継ぐcommit/証拠を明記する。
