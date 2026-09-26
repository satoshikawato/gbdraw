# INSTRUCTION PROMPT — S06: Time-bounded projection and restoration

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

S01/S05で測った残るlong tasksを除き、大規模Save/Load全pipelineのbudgetsとcontent equivalenceを満たす。

## 開始条件

S05 result、same-document gate、real fixture/measurement recipe、既存Vibrio acceptance。

## 所有範囲

session-file streaming scheduler、config/projection/validation/candidate preparation/restoreのprivate loops、必要な既存owner helpers、focused/performance tests、S06 result。second validator/renderer/transport/exportWorkerは追加しない。

他sessionが変更したcodeをrevertせず、最新の同名 branchに合わせてscope内を実装してください。
SOLID/KISS/DRY/YAGNIをstate/owners/paths/compatibility/deliveryにも適用し、
不要なframework、重複owner、並行 fallback、将来用 capabilityを追加しないでください。

## 実行手順

1. 保存projection、resource encoding、catalog validation、preflight、restore、sanitization/DOM mountを再測定し、unchanged証拠はreuse。
2. heavy loopsを時間予算とactual task yieldで分割。microtask/nextTickやbyte countだけでは合格にしない。
3. immutable adopted resources/catalog/manifestのreuseを保ち、全graph clone/sign/hash/base64 reencodeを増やさない。小mutable stateだけ必要範囲で固定。
4. candidate preparationのyieldでlive stateのpartial adoptionを公開しない。rollback authorityを一つに保ち、transaction boundaryを早く解除しない。
5. streaming exportを維持しJSON/gzip意味を比較。巨大strings/large object enumeration一つが時間予算を超える場合はそのownerの処理を分割し、silentfallbackを足さない。
6. real 10+ records / 25,000+ features/full pairwiseとVibrioでbefore/after、main+Worker/process peak、stage/wall、p95/max heartbeat、copy/decode/reencodeを比較。
7. preview DOM停止も独立報告し、ユーザーの閲覧・scroll/検索が実際に進むことを確認。未解消なら全体完了にしない。

## 検証と完了条件

S-01〜S-05。既存Vibrio performanceをthreshold変更なしで実行、new fixture/reference environment固定。同じ payload/request/catalog/cache/overrides/SVG、fresh Load→Generate→CLI/Python replayとfailure recovery。

共通作業規約のverification commandsからscopeに必要なものを実行し、必要な追加targetを結果に記録してください。
失敗を修正し、未測定をpassと報告しないでください。unchanged evidenceは条件一致時だけ再利用してください。

## 保存・コミット・プッシュ

required performance/equivalence gates、raw metrics/recipe、owner/path evidence、S06 resultをcommit/push。測定した限界と残件を定量記録する。

`results/S06_RESULT.md`にcommands/results、source/input/environment SHA、owners/paths、
authority/base references、acceptance ID、remaining boundaryと次sessionの具体的入口を保存してください。
このsessionの実装・tests・docs・結果だけを一commitにし、終了時に必ず同名remoteへpushしてください。
English commit title: `Keep large Session projection and restoration responsive`。

```bash
git branch --show-current
git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}'
git push origin HEAD:refs/heads/fix/issue-597-input-session-20260926
git ls-remote --heads origin refs/heads/fix/issue-597-input-session-20260926
git rev-parse HEAD
```

stage/commit前チェックとnon-fast-forward回復は共通作業規約に従います。force pushはしません。
local/remote SHA一致を確認してからbranch/commit/checks/残件を報告してください。
このpromptの次sessionまで連続実行せず、S06の範囲で終了してください。
