# INSTRUCTION PROMPT — S08: Integrated acceptance and final head handoff

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

全acceptanceを統合検証し、必要なarchitecture例外review用packetをfinal exact headに固定する。

## 開始条件

S00・S01・S03〜S07 results、required Product/privileged authorityがdevにmerged、missing independent checksなし。

## 所有範囲

in-scope final fixes、統合tests/evidenceとresults/S08_RESULT.md。policy/checker/authorityをgate通過のため変更しない。final-head packetはrepo外Markdownにする。

他sessionが変更したcodeをrevertせず、最新の同名 branchに合わせてscope内を実装してください。
SOLID/KISS/DRY/YAGNIをstate/owners/paths/compatibility/deliveryにも適用し、
不要なframework、重複owner、並行 fallback、将来用 capabilityを追加しないでください。

## 実行手順

1. D/S/A/W acceptance matrixをactual tests/resultsへ結び付け、未実施・既存reuse・failedを区別。不足はscope内で実装/検証して解消。
2. Circular/Linear の両方とsingle/grid/batch、native/restored/inactive/settings-only、draft/committed mismatch、History/rollback、full replayを確認。
3. production/tests/docs/generated diffを別々にreview。ordinary scopeはowner/path証拠、positive CB等の例外は完全な OE/PE/CB before/after sets、arithmetic、released fixtures、expiry/removal条件。
4. 必須の focused suites、通常 PR の browser contracts、real large/Vibrio performance、pytest の not slow 対象、ruff、必要build、architecture-contractsとtrusted base の policy checkを実行。既存有効な evidenceはcode/input/env/condition一致時だけreuse。
5. testtimeout/budgets/size/privacy/correctnessを緩めず、hard failureをProduct approvalでwaiveしない。必要別authorityが未mergedなら影響runtimeの完了を宣言しない。
6. final source/evidence/resultsを一commitにして同名 branchへpush、remoteSHA一致を確認。
7. 実際の変更でarchitecture例外条件が生じた場合のみ、そのfinalSHAを使い、repo外のMarkdownに全項目を記入した Architecture Ratchet review packetを生成。CI/results/setsとexact headを記載し、maintainerが手動承認する。agentはapproval commentを投稿しない。
8. approval後にdocs commitでheadを変えない。source/fixtures/conditionsが変わればmaterialchecksを再実行し新headのpacketを作る。PR/merge/deployは対象の別authorizationが必要。

## 検証と完了条件

全D-01〜D-04/S-01〜S-05/A-01/W-01。ordinary policy の Gate/Review結果とhumanarchitecture例外を別に記録。final remote SHA、cleancheckout、leaks/privileged paths/cycles、同じ内容/replayを確認。

共通作業規約のverification commandsからscopeに必要なものを実行し、必要な追加targetを結果に記録してください。
失敗を修正し、未測定をpassと報告しないでください。unchanged evidenceは条件一致時だけ再利用してください。

## 保存・コミット・プッシュ

S08resultとfinalfixesのcommit/push、例外条件がある場合のfinalSHAに束縛したrepo外packet、English summaryとremaining external review boundary。未承認例外をapprovedと書かない。

`results/S08_RESULT.md`にcommands/results、source/input/environment SHA、owners/paths、
authority/base references、acceptance ID、remaining boundaryと次sessionの具体的入口を保存してください。
このsessionの実装・tests・docs・結果だけを一commitにし、終了時に必ず同名remoteへpushしてください。
English commit title: `Complete input and Session regression acceptance for issue 597`。

```bash
git branch --show-current
git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}'
git push origin HEAD:refs/heads/fix/issue-597-input-session-20260926
git ls-remote --heads origin refs/heads/fix/issue-597-input-session-20260926
git rev-parse HEAD
```

stage/commit前チェックとnon-fast-forward回復は共通作業規約に従います。force pushはしません。
local/remote SHA一致を確認してからbranch/commit/checks/残件を報告してください。
このpromptの次sessionまで連続実行せず、S08の範囲で終了してください。
