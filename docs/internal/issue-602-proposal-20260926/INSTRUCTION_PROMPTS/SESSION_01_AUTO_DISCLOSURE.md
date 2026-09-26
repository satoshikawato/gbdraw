# Session 01 — Autoの既定を維持し、非表示の理由と変更先を示す

## Missionと製品結果

gbdraw Linearは複数recordをshared rowへ置くと、Autoを選んだAccession/Lengthを図全体で隠す。署名済み `linear.record-label-auto-visibility` revision 2 / **B / AUTO-FRESH-RESET-WITH-DISCLOSURE** はこの既定と解決規則を維持し、理由を配置操作の場所にも示す。

fresh/resetは各Auto、単一record rowsはShown、実際のshared rowが一つでもあればAutoのfieldだけdiagram-wide Hidden。Show/Hide、disabled layoutの休眠行除外、Session選択、保存Resultは維持する。説明は次の成功Generateのdraft効果であり、現在Resultを推測しない。

開始条件: S00の正式authorityがorigin/devへmerge済み。計画commitを持つ既存の実装ブランチを継続使用する。

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

- `app/linear-label-visibility.js` の既存resolver、`app/linear-record-layout.js` のrendered shared-row判定、`app/app-setup.js` の派生binding、`index.html` のLinear Layout/Record Labels表示。
- 必要なfocused Node/browser tests。default factoryは変更しない。
- Definition既定、global Status、mobile layout、Python rendererは後続の責務。per-row Auto、新collision engine、default Showは追加しない。

## 実装

1. resolverのAuto/Show/Hideとshared-row判定がrequest projection/表示で同じ意味かtraceする。別のshared-row検出やvisibility ownerを作らない。
2. LayoutとLabelsにAutoが隠すfieldを特定した説明を出す。理由はshared rendered row、範囲は図全体、適用は次回Generate。AccessionだけAuto等の混合でも対象を間違えない。
3. Record Labelsへ移動する明示buttonを置く。click/keyboardでnative disclosureを開き、対象selectへscroll/focusする。値をShowへ変えたりGenerateしたりしない。
4. Autoの選択と実効結果を分け、共有行解除やlayout disableで説明を再評価する。説明やnavigationにHistoryを作らない。未知modeの既存validationを維持する。
5. polite statusは重要な効果変更だけを伝え、繰返しreadoutやDOM textからの状態逆推定を避ける。

## 検証と受入

総合計画のV01–V03。single/shared/mixed rows、Auto/Show/Hide全組合せ、共有行解除、disabled休眠layout、Load済み旧Resultと変更後draftが異なるscenarioを含める。navigationの前後で選択・SVG・History件数が不変、focusが対象selectへ到達することをbrowserで示す。Generate後には説明したbooleanと実SVG textが一致する。

```bash
node --test tests/web/linear-label-visibility.test.mjs tests/web/linear-typography.test.mjs tests/web/session-draft-authority.test.mjs
npx playwright test tests/web/linear-typography.playwright.spec.js --workers=1
git diff --check
```

既存browser specに回帰を加え、実際の説明・keyboard移動を確認する。Python実装を触っていなければ現行geometryの証拠を再利用する。required Web checksはactual diffから実行する。

## 完了handoff

`docs/internal/issue-602-proposal-20260926/results/S01.md` にbase/head、branch/upstream、正式authority、変更対象、削除した旧処理、command/result、受入ID、visual observations、残る制約を記録する。未実行checkを合格扱いしない。最終回答は実現結果と確認済み範囲を説明し、Englishのproposed commit titleとsummaryを示す。次は `SESSION_02_DEFINITION_DEFAULT.md`。その開始条件を示し、利用者が次セッションを指定したらこの指示書をそのまま使える状態で渡す。
