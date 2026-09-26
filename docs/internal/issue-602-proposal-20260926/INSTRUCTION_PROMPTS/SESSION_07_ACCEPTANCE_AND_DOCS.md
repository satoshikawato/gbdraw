# Session 07 — 統合受入、ユーザー文書、最終handoffを完成する

## Missionと採択済み結果

gbdraw Issue #602の完成を検証する。Product Decision Owner `satoshikawato` / 2026-09-26の採択は**01-B、02-A、03-A、04-A、05-A**である。

- Accession/Lengthはfresh/reset Autoを維持し、共有行によるdiagram-wide非表示を説明する。
- Web fresh/resetのDefinition列はLock ON、明示OFF/保存/CLI/Python defaultは維持する。
- PendingとLive applying/error、Apply前reviewを区別し、Generate/Save/Exportの意味を表示する。
- 狭いPreviewは図と同一Editor/reviewの上下配置。review中はEditorをcloseしてtab保持、理由付きopen disable、終了後明示reopen。
- 狭いreviewのdragは退役するが、候補draft、方向設定、canvas、Apply/Cancel、失敗retryと既存History/Sessionは維持する。

開始条件: S01–S06完了、全対応authorityと旧PD-OI-035 mobile例外の限定supersessionがruntime baseにある。未merge authority、前段未検証、受入漏れを完了扱いしない。

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

## Ownership

- 総合計画のV01–V03、D01–D03、E01–E05、M01–M05、G01と現行ci-impactが要求するgates。
- ユーザー文書は `docs/REFERENCE/web-app.md`。既存の文書所有者を使い、新しいpublic guideを増やさない。
- 統合で見つかった不具合は対応する既存ownerで修正し、その範囲の証拠だけ更新する。authority/checkerをruntime候補に合わせて緩めない。

## 統合確認

1. latest origin/devとruntime ancestry、5件の全signed outcome、PD-OI-024のD2/D3保持、PD-OI-035旧mobile例外の終了を確認する。
2. 実用的なLinear SessionでAuto説明→共有行変更→Lock ON/OFF→Generate→Live編集→別設定Pending→Save/Load→Export→Generate→Undo/Redoを実操作する。Current Resultとdraftを各checkpointで比較する。Circularでは共通Status/Editorのnon-regressionも確認する。
3. resolved normal Align、明示resolved review、ambiguity、automatic failureからreview retryをnarrow/wideで実行する。Select/Skip/方向/marker、Apply/Cancel/失敗/stale、drawer調停とresizeを確認する。
4. failure/cancel/staleのResult/request/Historyと、直接live適用済み編集の保持を確認する。保存draftを適用済み基準にしたりStatus purposeのWorker/byte/hash/cloneを加えたりしていないことをtraceする。
5. geometryはactual SVG共通左端/gap、可視canvasのrect、elementFromPoint、zoom/pan、focusとscrollで確認する。snapshotやpanel幅だけで品質を判定しない。
6. desktop、390×844/740、390×500、320px幅、landscape、200% zoom、keyboard時の到達を確認する。自動viewport縮小と実機keyboardの証拠を区別し、残るlimitationsを明示する。
7. Web referenceへAuto理由と導線、Lock ON default、3種の適用分類、Pending/live状態、Generate再配置/zoom、Save/Export、compact操作を記載する。文書は実際のUI名と完了した結果を説明し、会話や選択過程を読者へ要求しない。
8. Gallery tutorialやpublic screenshotが直接affectedならweb-gallery-screenshot-maintenance skillを読む。手順文書の再生成証拠が必要ならlove-me-love-my-docs skillを適用する。実際のsource recipe/sessionとgeneratorで更新し、public図をminimal fixtureへ置換しない。

## 最終checks

前段の同一HEAD/入力/環境のfocused証拠を再利用し、materialに変わったscopeと未解決点だけ再実行する。次は最終gatesの基本commandであり、actual diffに対するcurrent CIのrequired setも確認する。

```bash
ruff check gbdraw/
node --test tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs --base origin/dev
python -m pytest tests/ -v -m "not slow"
python -m build
git diff --check
```

wheel依存browserチェックは先に `python tools/prepare_browser_wheel.py`。cache-bust refreshはdeployable bundleを準備する場合だけ。Python描画を変更した場合、またはgeometryに未解決の懸念があれば `python -m pytest tests/test_output_comparison.py::TestOutputComparison -v` をread-onlyで実行する。新しいintentional geometryが必要な場合だけreference再生成をreviewする。

changed focused Node/browser specsと既存CIのWeb/Python contractsを確認する。依存/privacy/Worker lifecycleが不変なら追加offline auditは不要。最後にproduction、tests、docs、generated diffsを独立reviewし、旧projection/位置処理/重複global Pendingが残っていないことを確認する。

## 完了の判定と報告

全受入条件が成立し、required gatesが合格、authorityがruntime baseに存在するときに完了とする。Review REQUIREDはGate PASSと別に扱い、必要なmanual review/例外を省略しない。既存scientific content/identityとsupported persistenceは維持、新依存/第二renderer/schema/互換pathはゼロであることを示す。

`docs/internal/issue-602-proposal-20260926/results/S07.md` にexact HEAD/base、authority、受入IDとcommand/artifact、visual observations、limitations、変更のない証拠の再利用根拠、owner/path/CB非増加を記録する。例外が必要ならcomplete OE/PE/CBと別maintainer判断を用意する。

最終回答はユーザーに見える完成結果、tests、material limitationsとEnglish proposed commit title/summaryを示す。公開権限がある場合だけ指定work branchへpush/PRを行い、retry前にremote stateを確認する。許可がなければlocal成果をレビュー可能な状態で渡し、push/mergeの未実施を機能実装の未完了と混同しない。
