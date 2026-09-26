# INSTRUCTION PROMPT — S05: Measured Session import Worker

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

S01で成立した一transportをproductionに導入し、main-thread JSON parseを削除する。

## 開始条件

S04 result、S01 transport/metrics、constructor/importer permissionsと必要mapped evidenceがtrusted baseにmerged。pathsがpreauthorizedと一致する。

## 所有範囲

session-file boundary、focused session-import-client/worker modules、config import coordination、lifecycle tests/required browser contract instrumentation、S05 result。guard/checker/authority変更を混ぜない。

他sessionが変更したcodeをrevertせず、最新の同名 branchに合わせてscope内を実装してください。
SOLID/KISS/DRY/YAGNIをstate/owners/paths/compatibility/deliveryにも適用し、
不要なframework、重複owner、並行 fallback、将来用 capabilityを追加しないでください。

## 実行手順

1. File/BlobをclientからJS import Workerへ渡し、read/gzip/fatal UTF8/JSON parseをWorkerで実行。main全量text生成を先行しない。
2. S01選択のwhole-objectまたはbounded sectionsの一方式。Vue/DOM/backing internalsを転送せずnew persisted protocolにしない。
3. operation IDs、structured errors、stale/cancel/teardown、settlement時terminate、transfer ownershipをclientの一ownerへ集約。
4. config.jsのdirect JSON.parseを削除。existing migrations/authority/preflight/adoptionを維持し、Worker返却をtrusted documentと見なさない。unsafe keysをcandidate使用前に検証。
5. limits/JSON-gzip/legacy/settings-only/unsupported browserの意味を維持。Worker failureをmain parse fallbackで隠さない。
6. Worker probesをSession importとPython diagramで区別し、saved-preview Python Worker 0 acceptanceは維持。
7. repeat importsでworker/objectURL/resourcesのleakを確認し、transfer/clone/heap metricsをS01と比較。

## 検証と完了条件

S-01/S-02とS-03/S-04 transport範囲。codec/unit、malformed/unsafe/oversize/fatal UTF8、crash/late/duplicate、small+real large JSON/gzip、current+released versions、fresh preview/replay。

共通作業規約のverification commandsからscopeに必要なものを実行し、必要な追加targetを結果に記録してください。
失敗を修正し、未測定をpassと報告しないでください。unchanged evidenceは条件一致時だけ再利用してください。

## 保存・コミット・プッシュ

選択transportだけのproduction/tests、removed main parse evidence、lifecycle/memory data、S05 resultをcommit/push。runtime方式未成立なら原因を直してfocused checksを通す。

`results/S05_RESULT.md`にcommands/results、source/input/environment SHA、owners/paths、
authority/base references、acceptance ID、remaining boundaryと次sessionの具体的入口を保存してください。
このsessionの実装・tests・docs・結果だけを一commitにし、終了時に必ず同名remoteへpushしてください。
English commit title: `Move Session import parsing to a bounded browser Worker lifecycle`。

```bash
git branch --show-current
git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}'
git push origin HEAD:refs/heads/fix/issue-597-input-session-20260926
git ls-remote --heads origin refs/heads/fix/issue-597-input-session-20260926
git rev-parse HEAD
```

stage/commit前チェックとnon-fast-forward回復は共通作業規約に従います。force pushはしません。
local/remote SHA一致を確認してからbranch/commit/checks/残件を報告してください。
このpromptの次sessionまで連続実行せず、S05の範囲で終了してください。
