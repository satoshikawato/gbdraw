# INSTRUCTION PROMPT — S00: Authority intake and approved decision serialization

あなたは gbdraw の Issue #597 修正担当者です。このpromptとrepository内の資料だけで作業してください。
実装対象は `fix/issue-597-input-session-20260926`。他sessionのcheckout/branchを変更してはいけません。

## ブランチ取得と必読資料

最初に [SESSION_WORKFLOW.md](../SESSION_WORKFLOW.md#1-独立-checkout-に取得する) の手順で、
remoteの最新 `fix/issue-597-input-session-20260926` をこのsession専用の独立checkoutへclone/fetchして使ってください。
shared checkoutをswitch/reset/cleanせず、別の新規実装branchを使わないでください。
`AGENTS.md`、`CLAUDE.md`、`gbdraw/web/CLAUDE.md`、[MASTER_PLAN.md](../MASTER_PLAN.md)、
共通作業規約、[承認済みdecisions](../README.md)を読んでください。S00には前提sessionのresultはありません。
Product outcomesは承認済みで、再選択を求めません。authority-before-runtime等の別条件は維持してください。

## 目的

2件のapproved receiptsを既存authorityに統合する候補を完成させる。Product再承認を要求せず、runtimeには着手しない。

## 開始条件

計画branchがremoteに存在する。2件のdecisions Markdownのreceipt/machine fieldsが一致する。

## 所有範囲

`authority-candidates/` と `results/S00_RESULT.md`。2026-09-26の追加指示により、この計画ディレクトリ内の BUG-01 計画・依存の削除も担当する。残る2件の承認本文は変更しない。既存policy/map/contract/first-parent historyはread-only。active guard/checker/runtimeをこのbranchで編集しない。

他sessionが変更したcodeをrevertせず、最新の同名 branchに合わせてscope内を実装してください。
SOLID/KISS/DRY/YAGNIをstate/owners/paths/compatibility/deliveryにも適用し、
不要なframework、重複owner、並行 fallback、将来用 capabilityを追加しないでください。

## 実行手順

1. Product Impact map/durable store/static Product Contract/privileged detectorを読み、各concernの正しいauthority targetを確定する。unmapped concernを無理にBDへ登録しない。
2. selected choice、rationale/mustPreserve/mayRetire/risk/owner/dateをそのままauthority-only patchへserialize。静的contractのrevision/IDは最新baseから割り当て、既存scopeを変更しない。
3. concrete import-client/Worker pathsのprivileged subjectsとimportersを調査する。初期候補clientは`services/session-import-client.js`、workerは`workers/session-import-worker.js`。constructorの一ownerと実際に必要な最小permissionsを別patchにする。S01でtransportが確定してからfinal pathsをreviewする。
4. mapped behavior contractが変わる場合のevidence-only→authority-ref→runtime順序を整理。checker改修が必要ならchecker-only laneを分ける。
5. separate authority-only PRのexact target、patch application/validation commands、merge後に記録するSHAをresultsへ書く。候補patchの統合はmaintainerが最新devの別PRで行う。implementation branchをauthority-only PRのbaseにしない。
6. 候補をtemporary clean dev worktreeに適用し、schema/field-preservationとdiff pathsを検証してdiscard。active file変更をimplementation branchへ持ち込まない。独立worktree以外を操作しない。

## 検証と完了条件

2 receiptsとcandidateのfield equality、patch checkとallowable paths、既存authorityへのconflictなし、Worker privileged scope根拠。authority統合は未完了ならpendingと記録する。

共通作業規約のverification commandsからscopeに必要なものを実行し、必要な追加targetを結果に記録してください。
失敗を修正し、未測定をpassと報告しないでください。unchanged evidenceは条件一致時だけ再利用してください。

## 保存・コミット・プッシュ

inert candidate patches、receipt comparison、external merge targetsとreadiness表、S00 result。完成したpreparationをcommit/pushし、S01へ進める。

`results/S00_RESULT.md`にcommands/results、source/input/environment SHA、owners/paths、
authority/base references、acceptance ID、remaining boundaryと次sessionの具体的入口を保存してください。
このsessionの実装・tests・docs・結果だけを一commitにし、終了時に必ず同名remoteへpushしてください。
English commit title: `Document approved issue 597 authority integration`。

```bash
git branch --show-current
git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}'
git push origin HEAD:refs/heads/fix/issue-597-input-session-20260926
git ls-remote --heads origin refs/heads/fix/issue-597-input-session-20260926
git rev-parse HEAD
```

stage/commit前チェックとnon-fast-forward回復は共通作業規約に従います。force pushはしません。
local/remote SHA一致を確認してからbranch/commit/checks/残件を報告してください。
このpromptの次sessionまで連続実行せず、S00の範囲で終了してください。
