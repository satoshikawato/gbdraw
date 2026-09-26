# INSTRUCTION PROMPT — S02: Circular sources and canonical Session bindings

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

独立ordered source collectionをfresh input、request、writer/reader、Historyで一度に整合させる。

## 開始条件

S00/S01 results、source Product authorityと必要permissionsがorigin/devにmerged。merge SHAを確認してbranchへ統合する。

## 所有範囲

state.js、components.js/index.html input cards、app/circular-sources.js、record-display-options.js、run-analysis.js/watchers.js/setup wiring、session-request/resources/resource-backing/config/active-contract/history、必要なtyped Session読取、focused tests、S02 result。guard/checker/static authorityは変更しない。

他sessionが変更したcodeをrevertせず、最新の同名 branchに合わせてscope内を実装してください。
SOLID/KISS/DRY/YAGNIをstate/owners/paths/compatibility/deliveryにも適用し、
不要なframework、重複owner、並行 fallback、将来用 capabilityを追加しないでください。

## 実行手順

1. `{uid,file}` collectionとper-source discoveryを一ownerにする。add末尾、replace同位置新UID、remove対象UID。counts/flat catalogはderived。
2. existing parserとuploaderを再利用。run-analysisの旧独立discovery lifecycle/scalar runtime経路を削除する。restore/reset/Historyはowner transitionsを明示呼出。
3. UID/local selector→stable recordKey、request-local index→resourceというmappingを統一。same-name/same-bytes/duplicate accessionsを別instanceにする。
4. fresh/restoreを既存circularRecords投影へ収束。single/grid/batch、saved opt-out、single crop、transforms/positions/depth/annotations/prefixes/endpointsを同じordered universeから投影。
5. bindings3 writerとreleased1/2 readerのdirect normalization。current c_gb writerを除去、legacy compositeはtruthful mappingだけ分離、曖昧combinedは一source/backingとして維持。branch-only intermediatesのreaderを作らない。
6. before-first-Generate、inactive-mode sources、draft/committed mismatch、settings-onlyを保存。source交換対象だけintentをinvalidate、imported comparisonをunresolved ownerに返す。
7. add/replace/removeは一History operation。async discoveryは別Historyを作らず、stale/Undoをsource versionでguard。sourceエラーでcommitted Resultを消さない。
8. schema namespace evidence、complete owner/path/CB setsとremoved pathsをS02 resultへ記録。positive CBならfinal-head architecture reviewの必要性を残す。
9. new current writerで拒否されるbranch-owned Sessionがrequired testsの入力なら、このsessionでowner generatorによりcurrent formatへ再生成してcommitへ含める。S07へ延期するためのtemporary reader/skipは作らない。S07には再生成evidenceと残るpublic docs/capture対象を渡す。

## 検証と完了条件

C-01〜C-05とA-01。session-request/resources/backing/active-files/record-display/typed-session tests、browser multi-file fresh/restored/undo、CLI/Python replay。別sourceのsparse depth/selector/rotationをassert。

共通作業規約のverification commandsからscopeに必要なものを実行し、必要な追加targetを結果に記録してください。
失敗を修正し、未測定をpassと報告しないでください。unchanged evidenceは条件一致時だけ再利用してください。

## 保存・コミット・プッシュ

一貫したproduction/testsとS02 resultをcommit/push。partial new writerや並行scalar pipelineの状態で終了しない。required inputsのcurrent artifact再生成もこのcommitに含める。S07へ残る範囲とreuse可能evidenceを列挙する。

`results/S02_RESULT.md`にcommands/results、source/input/environment SHA、owners/paths、
authority/base references、acceptance ID、remaining boundaryと次sessionの具体的入口を保存してください。
このsessionの実装・tests・docs・結果だけを一commitにし、終了時に必ず同名remoteへpushしてください。
English commit title: `Add source-bound Circular multi-file input and session bindings`。

```bash
git branch --show-current
git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}'
git push origin HEAD:refs/heads/fix/issue-597-input-session-20260926
git ls-remote --heads origin refs/heads/fix/issue-597-input-session-20260926
git rev-parse HEAD
```

stage/commit前チェックとnon-fast-forward回復は共通作業規約に従います。force pushはしません。
local/remote SHA一致を確認してからbranch/commit/checks/残件を報告してください。
このpromptの次sessionまで連続実行せず、S02の範囲で終了してください。
