# INSTRUCTION PROMPT — S01: Baseline and import transport evidence

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

現状のstage bottlenecksと一つのimport transportの実現性を測る。production code/default/authorityを変更しない。

## 開始条件

S00 result。必要authority未mergeでも独立measurementは進められる。基準runtimeのsource fingerprintを固定する。

## 所有範囲

focused characterization tests/measurement recipes、`evidence/`、`results/S01_RESULT.md`。isolated probeはtools/testsの既存ownerに置き、application routingには接続しない。mapped contract変更が必要なら別evidence-only integration対象として準備する。

他sessionが変更したcodeをrevertせず、最新の同名 branchに合わせてscope内を実装してください。
SOLID/KISS/DRY/YAGNIをstate/owners/paths/compatibility/deliveryにも適用し、
不要なframework、重複owner、並行 fallback、将来用 capabilityを追加しないでください。

## 実行手順

1. GenBank/DDBJ/GFF pair、duplicate IDs、rapid replace、saved artifactとdraftの差、settings-onlyの基準を固定する。既存有効evidenceは再実行せず不足だけ追加。
2. main first-parent/tagsでSession/request/bindings/catalog/cacheのnamespaceを確認し、released positive fixturesとSHAを記録する。
3. pinned biological sourcesからreal 10+ records/25,000+ features/full pairwise Sessionの再生成recipeを用意。既存Vibrio fixtureも使用。synthetic duplicationをreal genomesの代用と呼ばない。
4. Save/Loadのread/decompress/parse/transfer/preflight/projection/encode/validation/restore/sanitize/DOM mountを分け、100ms heartbeat、p95/max、longtasks、wall、main+Worker/process memory、bytes/copy/reencodeを測る。
5. disposable File/Blob→JS Worker parse→whole-object replyを比較。不適合ならbounded known sections/feature batches/huge strings transferを比較し、実装用に一方式を選ぶ。main全量text、proxy転送、fallback併存を採用しない。
6. 必要なprivileged pathsをS00候補に反映し、baseline/environmentに対するbudgets、accepted transport、missing evidenceをresultsに記す。目標p95≤250ms/max≤500msを達成済みと捏造しない。

## 検証と完了条件

fixture generation correctness、stage/heap metrics、structured clone停止の有無、payload/request/catalog/cache/SVG equivalence。Vibrio既存budgetsを緩めない。full pipelineとcodec-onlyを区別する。

共通作業規約のverification commandsからscopeに必要なものを実行し、必要な追加targetを結果に記録してください。
失敗を修正し、未測定をpassと報告しないでください。unchanged evidenceは条件一致時だけ再利用してください。

## 保存・コミット・プッシュ

reproducible fixtures/recipes/metrics、chosen transportと理由、namespace evidence、S01 result。方式未成立ならevidenceをcommit/pushしS05依存部分を止める。

`results/S01_RESULT.md`にcommands/results、source/input/environment SHA、owners/paths、
authority/base references、acceptance ID、remaining boundaryと次sessionの具体的入口を保存してください。
このsessionの実装・tests・docs・結果だけを一commitにし、終了時に必ず同名remoteへpushしてください。
English commit title: `Measure issue 597 session responsiveness and transfer costs`。

```bash
git branch --show-current
git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}'
git push origin HEAD:refs/heads/fix/issue-597-input-session-20260926
git ls-remote --heads origin refs/heads/fix/issue-597-input-session-20260926
git rev-parse HEAD
```

stage/commit前チェックとnon-fast-forward回復は共通作業規約に従います。force pushはしません。
local/remote SHA一致を確認してからbranch/commit/checks/残件を報告してください。
このpromptの次sessionまで連続実行せず、S01の範囲で終了してください。
