# INSTRUCTION PROMPT — S04: Session operation consistency

あなたは gbdraw の Issue #597 修正担当者です。このpromptとrepository内の資料だけで作業してください。
実装対象は `fix/issue-597-input-session-20260926`。他sessionのcheckout/branchを変更してはいけません。

## ブランチ取得と必読資料

最初に [SESSION_WORKFLOW.md](../SESSION_WORKFLOW.md#1-独立-checkout-に取得する) の手順で、
remoteの最新 `fix/issue-597-input-session-20260926` をこのsession専用の独立checkoutへclone/fetchして使ってください。
shared checkoutをswitch/reset/cleanせず、別の新規実装branchを使わないでください。
`AGENTS.md`、`CLAUDE.md`、`gbdraw/web/CLAUDE.md`、[MASTER_PLAN.md](../MASTER_PLAN.md)、
共通作業規約、[承認済みdecisions](../README.md)、前提sessionの`results/Sxx_RESULT.md`を読んでください。
Product outcomesは承認済みで、再選択を求めません。authority-before-runtime等の別条件は維持してください。

## 目的

Save/Load中のsemantic mutationを一つのavailabilityで止め、read-only browsingとsame-document consistencyを保つ。

## 開始条件

S03 result、session-operation Product authorityがmerged。Generate/reflowとsession lifecycleを監査する。

## 所有範囲

app-setup/session actions、config.js、Generate/source/editor/History/Reset/cache/file-import ownersとtemplates、focused tests、S04 result。immutable snapshot/lock frameworkは追加しない。

他sessionが変更したcodeをrevertせず、最新の同名 branchに合わせてscope内を実装してください。
SOLID/KISS/DRY/YAGNIをstate/owners/paths/compatibility/deliveryにも適用し、
不要なframework、重複owner、並行 fallback、将来用 capabilityを追加しないでください。

## 実行手順

1. semantic mutation inventoryをaction/DOM/programmatic/async completion別に列挙。pan/zoom/scroll/search等read-onlyを区別。
2. existing save/import pending flagsからavailabilityを導出し、owning actionsとUIに同じpredicateを使う。第二lock refを持たない。
3. same Save join、一download、異種busy、Generate/reflow実行中のSave/Load busy理由とretryを実装。
4. title確定後のpublish/paint/startからsettlementまで一document。遅延file import/cache fill/editor reflowもcheckpointをまたいでpartial commitしない。
5. programmatic mutationsにもexplicit busy result。error/cancel/crash/teardownでpendingを解除し、source first rollback→transient reconcile。
6. private candidateとlive artifactを分離し、failed adoption前後でrequest/resources/Result/Historyの同一性を保護する。

## 検証と完了条件

S-02。session-save-lifecycle/loading-feedback、Generate/reflow busy、duplicate/cross operations、all mutation owners、late completions、read-only pointer/keyboard、failure rollback、settings-only/draft mismatch。

共通作業規約のverification commandsからscopeに必要なものを実行し、必要な追加targetを結果に記録してください。
失敗を修正し、未測定をpassと報告しないでください。unchanged evidenceは条件一致時だけ再利用してください。

## 保存・コミット・プッシュ

mutation inventory/checkpoint evidence、production/tests、S04 resultをcommit/push。disabled HTMLだけで合格にしない。

`results/S04_RESULT.md`にcommands/results、source/input/environment SHA、owners/paths、
authority/base references、acceptance ID、remaining boundaryと次sessionの具体的入口を保存してください。
このsessionの実装・tests・docs・結果だけを一commitにし、終了時に必ず同名remoteへpushしてください。
English commit title: `Keep Save and Load consistent across browser mutations`。

```bash
git branch --show-current
git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}'
git push origin HEAD:refs/heads/fix/issue-597-input-session-20260926
git ls-remote --heads origin refs/heads/fix/issue-597-input-session-20260926
git rev-parse HEAD
```

stage/commit前チェックとnon-fast-forward回復は共通作業規約に従います。force pushはしません。
local/remote SHA一致を確認してからbranch/commit/checks/残件を報告してください。
このpromptの次sessionまで連続実行せず、S04の範囲で終了してください。
