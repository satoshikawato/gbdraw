# Session 00 — 署名済みProduct Decisionsを正式authorityへ反映する

## Missionと採択済み結果

gbdraw Issue #602はLinearのmetadata非表示説明、Definition列既定、編集適用の理解、狭いPreviewのEditor/review占有を修正する。このセッションはProduct authorityだけを扱う。

Product Decision Owner `satoshikawato` が2026-09-26に次の**全文**へ署名済みである。rationale、must preserve、may retire、accepted residual risk、owner/dateを含む正確な回答は `docs/internal/issue-602-proposal-20260926/06_SIGNED_PRODUCT_DECISIONS.md`。

| Concern / revision | 採択 |
| --- | --- |
| linear.record-label-auto-visibility / 2 | B / AUTO-FRESH-RESET-WITH-DISCLOSURE |
| linear.definition-display / 2 | A / D1-B-WEB-LOCKED-FRESH; retain D2-P and D3-A |
| web.edit-application-feedback / 1 | A / DERIVED-APPLICATION-STATUS |
| web.editor.compact-presentation / 1 | A / DOCKED-COMPACT-EDITOR |
| web.similarity-alignment.review-presentation / 1 | A / DOCKED-COMPACT-ALIGNMENT-REVIEW |

Accession/Lengthのfresh/resetをShowへ変えない。02はWeb fresh/resetのみLock ON、他surfaceと保存値を維持。03は操作単位の説明と派生Pending、04/05は同一図/同一panelの上下配置。05は狭いreviewのdragと同時Editor openを限定退役し、canvas、local draft、Apply/Cancel/retryを維持する。

## 共通の作業規則

- `AGENTS.md`、`CLAUDE.md`、`gbdraw/web/CLAUDE.md`、`docs/internal/PRODUCT_IMPACT_RATCHET.md`、`docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md` を読む。総合計画は `docs/internal/issue-602-proposal-20260926/00_MASTER_PLAN.md`、署名全文と機械表現は同directoryの `06_SIGNED_PRODUCT_DECISIONS.md`。別worktreeの場合は計画checkoutを参照元に保持する。
- 最初にworking tree、branch/upstream、前段commitと証拠を確認し、originをfetchする。他者の変更をpreserveし、main/devへcommitしない。実行する指示書の範囲でlocal編集・検証・修正を完了する。
- S00では署名されたconcern/revision/全寄与を正式authorityへ反映する。runtimeを扱うS01以降は、そのauthorityがbaseへmerge済みであることを確認する。署名済み選択を再度選び直したり、一般的な再承認を求めたりしない。未merge authorityや候補文書はruntimeを自己承認しない。
- 既存owner/pathを再利用し、canonical request/Worker/Result admission/Historyの第二経路を作らない。純粋な処理と副作用、適用の意味と表示を分離する。重複処理は同じ変更で除去する。新schema/互換reader/依存/汎用frameworkは追加しない。
- browser確認前にNodeとPythonのPlaywrightを確認する。Nodeがなければ同等のPython確認を実行する。Chromium sandboxだけの失敗は同じ確認を適切なescalationで再実行する。port競合は空きportを使い、既存serviceを停止しない。
- focused checksと必要gatesを実行し、失敗を修正する。同一コード/入力/環境の検証済み証拠は再利用する。長いtestは最低30分を許容して増分monitorし、test所有の短いtimeoutは弱めない。
- production、tests、docs、generated artifactsの差分を分けてレビューする。参考SVGは通常read-only、生成wheelはgitignored、social previewはowner-maintained。新たなProduct結果やarchitecture例外が必要なら該当範囲だけを停止し、独立作業を完了する。
- 検証済み変更を一つのlogical commitとして扱い、branch/upstreamを確認してEnglishのproposed commit titleとsummaryを提示する。push/PR/merge/deploy/tagは個別に明示された許可の範囲だけ。PRを作成・編集する場合は指定PR skillを読む。

## Ownership / 非対象

唯一の実装対象は `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`。既存の広い/unmapped concerns用static authorityである。runtime、tests、checker、registry、workflow、計画書を同じauthority-only候補へ含めない。新decision storeやBD番号を捏造しない。正式Contract以外のhandoffは元の計画checkoutへ保存する。

## 手順

1. 基準SHA `d457b7189b137185a8dec800819a312c30b969fa` と最新origin/devを比較し、allowlistで署名者、現在のconcern/revision、すでに同内容がmergeされていないかを確認する。5回答と機械表現をfield単位で照合する。一般的な再署名は求めない。
2. 最新origin/devから `product/issue-602-decisions-20260926` をupstreamなしで作る。実装ブランチ `fix/issue-602-linear-live-edit-20260926` は計画commitを保持したまま残し、別worktreeを使用する。計画文書をauthority-only候補へ混ぜない。実装ブランチのcheckoutにある計画/署名ファイルを参照する。
3. 現在のProduct Impact mapが該当結果をfaithfully所有しているか再照合する。基準SHAではunmappedのため既存static Contractを使う。新baseでmappedへ変わっている場合だけ、正しい既存authority routeへ調整する。候補の自己承認や別storeを作らない。
4. 5件の署名本文を正確にserializeする。新record番号はlatest baseで重複を確認して割り当てる。PD-OI-024のD1はWeb fresh/reset部分だけ改訂し、D2-P/D3-Aの独立寄与を残す。
5. PD-OI-035 scenario 2の390px遮蔽許可と「図を確認するにはreviewを閉じる」continuationは、署名済み05-Aの可視canvas/同時操作へ限定的にsupersedeする。新review-presentation recordへの適用関係を明示し、旧mobile例外をactiveな許可として残さない。identity、keyboard/Skip、非描画候補、desktop canvas、focus、overlay非保存はそのまま維持する。新しいrationale/riskを推測せず、05-Aの本文を根拠にする。
6. Contract revision、選択/署名metadata、supersessionと必要なacceptance referencesを整える。old Autoの意味や他のalignment結果を変更しない。署名本文と新machine JSONを元checkoutのS00handoffで並べてレビュー可能にする。
7. authority-only diffを独立reviewし、既存checkerを実行する。全差分が正式Contractだけであることを確認する。base checkerをcandidate checkerへ置換しない。

## 検証

```bash
node --test tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs --base origin/dev
git diff --check
```

commandと実際のbase/headを記録する。Contract checkerが追加されていればcurrent baseのrequired checkも実行する。署名回答の9項目が完全一致し、未選択01-Aを採用していないことを確認する。checkerだけの成功はmerge済みauthorityではない。

## 終了境界

検証済み候補をauthority branchの一つのlogical commitにする。push/PR/mergeはその具体的な外部操作の明示許可がある場合だけ。リモート操作をretryする前に実状態を確認する。S01の開始条件は**5結果と旧mobile例外の限定supersessionがorigin/devへmerge済み**であること。local候補があるだけなら独立baseline収集を除いてruntimeを開始しない。

## 完了handoff

`docs/internal/issue-602-proposal-20260926/results/S00.md` にbase/head、branch/upstream、正式authority、変更対象、削除した旧処理、command/result、受入ID、visual observations、残る制約を記録する。未実行checkを合格扱いしない。最終回答は実現結果と確認済み範囲を説明し、Englishのproposed commit titleとsummaryを示す。次は `SESSION_01_AUTO_DISCLOSURE.md`。その開始条件を示し、利用者が次セッションを指定したらこの指示書をそのまま使える状態で渡す。
