# Session 03 — requestとStatusで共用する生成意味と比較基準を整える

## Missionと製品結果

gbdraw WebのResultは最後に成功した図、生成設定draftは次のGenerate用の設定である。署名済み `web.edit-application-feedback` revision 1 / **A / DERIVED-APPLICATION-STATUS** は、この差を事実から表示する。生成設定PendingとLive applying/errorは独立であり、右の色編集がscaleのPendingを消してはいけない。

このセッションはcanonical projection、artifactの適用済み基準、比較の整合性を担当する。UI表示はS04。03-Aの適用タイミングを変更せず、右のlive actionは即時canonical commit、reviewはApply前local draft、Generateはatomic replacementのままにする。

開始条件: S02完了、03-Aが正式authorityとしてbaseにある。

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

- `services/session-request.js` のcanonical正規化/record/track/config projection、既存resource backing/token、既存artifact/History/session ownerの適用済み表現。
- `record-display-options.js` とPaletteのpartial Pendingはtrace対象。target-only操作の可用性など別責務を一括削除しない。
- request組立とStatusの共通意味を一つのownerへ集約。global UIや各inputのイベントへ個別dirty flagを追加しない。
- 新request/schema、SVG clone、genome hash、File I/O、Worker呼出、generic event bus、第二artifact checkpointは追加しない。

## 実装

1. current request builderと保存Result/editor overrides/active draftの表現をtraceする。実際にrequestへ影響する全domainを既存projectionに対応付ける。normalized current writer contractと既存typed inputを根拠にし、別render-field inventoryをStatusだけに作らない。
2. representative inputでbuilderのrequest意味をcharacterizeする。正規化を抽出した後も同じrequestを得ることを検証する。record topology、tracks、comparison、label/render options、unmanaged有効設定を漏らさない。
3. DOM/Worker/File内容不要の軽いprojectionを同じcanonical ownerへ置き、builderとStatus比較で共用する。UIはhelperを利用するだけでbuildCanonicalRenderRequestの新privileged callerにならない。
4. 比較するのは現在modeの有効描画意味とresource binding。tab/scroll/focus/休眠slot/他mode profile/保存メタデータは除く。Auto/Show/Hideの選択は保存状態に残し、現在のeffective booleanが同じならそれだけでPendingにしない。shared row変更時は同じresolverで再評価する。
5. Fileは既存backing/token/bindingで区別し、同名・同サイズの別入力をcleanとしない。Fileを再読取り/hash/base64化しない。resource identityが確定していない場合はunknownとして扱い、Statusのためにdiscovery Workerを開始しない。
6. 比較基準は既存artifact ownerに所属させ、保存requestと適用済みeditor状態を使う。現在draftで基準を更新しない。Generate成功、Session restore、Undo/Redoでは対応artifactへ切替。Live commitは実際に適用したfieldだけ反映し、別のPendingを残す。
7. 保存Resultとactive draftが違うSessionをLoadしても、そのdraftをapplied基準にしない。既存saved overrideから復元できる部分は復元し、不明はunknown、invalidはinvalid。未生成は未生成として表す。falseの安心表示を作らない。
8. full checkpoint/signature/cloneを増やさず、既存pre-stateを唯一のrollback authorityとして保つ。obsolete projectionや同じglobal Pending ownerを同じ変更で除去する。partial scopeの表示はS04へ担当範囲を引き継ぐ。

## 検証と受入

E01/E03/E04のデータ層。scale/crop/slots、active comparison、metadata、source交換、inactive mode/slot、pending Session、Live palette、Undo/Redo、Generate失敗/Cancel/staleを含む。元の値へ戻して差が解消すること、同名別Fileを識別することを確認する。invalid/unknownをcleanにしない。field coverageと適用済み基準の更新範囲を実装レビューで説明する。

```bash
node --test tests/web/session-request.test.mjs tests/web/session-draft-authority.test.mjs tests/web/session-active-config-contract.test.mjs tests/web/candidate-render.test.mjs tests/web/history-canonical-owner.test.mjs tests/web/record-display-options.test.mjs
node --test tests/web/architecture-contracts.test.mjs
node tools/check-web-change-budget.mjs --base origin/dev
git diff --check
```

focused testsへmeaningful regressionsを加える。Status操作のWorker/File-byte/hash/SVG-clone呼出がゼロであることをspy/既存structural hooksで確認する。browserの保存/履歴fixtureも一つ以上でartifact/draftの継続を確認し、S04の表示wiring用に比較結果とoperation factsの小さなinterfaceを渡す。

## 完了handoff

`docs/internal/issue-602-proposal-20260926/results/S03.md` にbase/head、branch/upstream、正式authority、変更対象、削除した旧処理、command/result、受入ID、visual observations、残る制約を記録する。未実行checkを合格扱いしない。最終回答は実現結果と確認済み範囲を説明し、Englishのproposed commit titleとsummaryを示す。次は `SESSION_04_APPLICATION_FEEDBACK.md`。開始条件と引き継ぐcommit/証拠を明記する。
