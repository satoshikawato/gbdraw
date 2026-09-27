# Session 02 — WebのDefinition列を既定Lock ONにする

## Missionと製品結果

Linear比較図のrowをalign/offsetするとLock OFFのDefinitionも追従する。署名済み `linear.definition-display` revision 2 / **A / D1-B-WEB-LOCKED-FRESH; retain D2-P and D3-A** は、Web fresh/resetをLock ONにする。

保存false/true、対応済み旧省略意味と保存Result、CLI/Python省略default、明示OFFの共通幅中央とrow追従は維持する。PD-OI-024のSubtitle保持/継承/ラベル区別と、Replicon/Organelleの選択順・独立制御・既定falseをすべて残す。Lock ON座標計算は基準SHAで修正済みであり、新たな配置engineを作らない。

開始条件: S01完了、対応authorityがbaseにある。

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

- `services/session-active-config-contract.js` のdefault factory、`services/reset.js` の利用、既存reader/request projectionのdefault/省略境界、Linear Layoutの短い説明。
- existing focused defaults/restore/request/browser tests。
- Autoのdefault/解決規則、CLI/Python default、新schema/reader、alignment biology、global Status、mobile layoutは変更しない。

## 実装

1. factoryがfresh/reset以外のlegacy restoreでも使われるかtraceし、保存false/trueとsupported old omissionのpositive fixturesを先に確認する。new defaultが古いSessionへ漏れる経路を特定する。
2. Web fresh/reset値を既存factoryでtrueにする。既存reader境界で旧省略意味を保持し、別互換pipelineを作らない。current malformed inputをdefaultで受理しない。
3. LayoutにONは共通左列、OFFはrow追従、変更はGenerateで適用、という短い常時説明を置く。
4. 既存共通placement、final translations後の列origin、collision bandsとpaintの共有を保持する。failが既存契約違反ならこのownerでfocused correctionを行い、その実SVG影響を確認する。tests/reference_outputsを成功のために書き換えない。

## 検証と受入

総合計画D01–D03。fresh/reset ON、保存OFF/ON、supported省略、読込Result不変、save/load/regenerationを確認する。CLI/Pythonの省略はfalseのまま。Subtitle、Replicon、record-local labels、text_anchor受入範囲も確認する。

正/負の不均等translation、center/Similarity alignment、単一/共有/混在行で、actual SVGのDefinition左端とnearest sequenceへのgapをChromiumで測る。text-anchor属性だけでは証明しない。既存nonbrowserテストの変更のない証拠は再利用できるが、新default/restoreは新しいassertionが必要。

```bash
node --test tests/web/session-active-config-contract.test.mjs tests/web/session-request.test.mjs
python -m pytest tests/test_linear_definition_alignment.py -q
npx playwright test tests/web/linear-multi-record.playwright.spec.js tests/web/linear-typography.playwright.spec.js --workers=1
git diff --check
```

browser wheelが古ければ通常prepare commandで準備する。Pythonを変更した場合はRuffとread-only OutputComparisonも必要。

## 完了handoff

`docs/internal/issue-602-proposal-20260926/results/S02.md` にbase/head、branch/upstream、正式authority、変更対象、削除した旧処理、command/result、受入ID、visual observations、残る制約を記録する。未実行checkを合格扱いしない。最終回答は実現結果と確認済み範囲を説明し、Englishのproposed commit titleとsummaryを示す。次は `SESSION_03_CANONICAL_STATUS_PROJECTION.md`。その開始条件を示し、利用者が次セッションを指定したらこの指示書をそのまま使える状態で渡す。
