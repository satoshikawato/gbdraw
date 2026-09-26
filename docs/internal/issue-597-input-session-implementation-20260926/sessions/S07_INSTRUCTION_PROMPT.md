# INSTRUCTION PROMPT — S07: Current artifacts and reproducible public documentation

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

current writerへのartifact更新と既存public workflow説明を完成し、branch-only intermediateの残存を除く。

## 開始条件

S06 result、current bindings/Session grammar、generated artifact inventory。main/tag historyを再確認。

## 所有範囲

既存Gallery generator/refreshが所有するSession/artifact projections、既存public technical/input/session pages、必要なtutorial screenshots/text、compatibility docs、S07 result。social_previewやdist/egg-infoを手編集しない。

他sessionが変更したcodeをrevertせず、最新の同名 branchに合わせてscope内を実装してください。
SOLID/KISS/DRY/YAGNIをstate/owners/paths/compatibility/deliveryにも適用し、
不要なframework、重複owner、並行 fallback、将来用 capabilityを追加しないでください。

## 実行手順

1. reader/writer namespaceとreleased fixture inventoryを確認。未公開Session44/bindings2はcurrentに再生成し、その組の不要reader/migrator/test/docsを除去。
2. Gallery canonical EXAMPLES inventoryとrefresh ownerを使用。S02でcurrent formatへ再生成済みのartifact evidenceはcode/input/environment/conditionsが同じならreuseし、不足したartifactsだけowner generatorで更新する。手編集しない。source bytes、metadata/labels/legend/tracks/comparison contextを保持。
3. 手順docsを変更する場合love-me-love-my-docs、Gallery/tutorial cropsならweb-gallery-screenshot-maintenanceを読んで適用。internal proseだけにcapture workflowを適用しない。
4. existing public page ownersへmulti-file、truthful discovery/single crop、exclusive session semantics、limits/errorsを記述。capabilityごとの新ページを増やさない。
5. documented commands/session/GUI actionsをclean checkoutから実行。readable scale/keyboard/390pxで最終public artifactsを視認しregeneration evidenceを残す。
6. generated wheelは必要時にprepareするだけ。deployable bundle準備の場合だけcache-bust更新。offline依存/lifecycleのruntime変更の検証は該当skillを使う。

## 検証と完了条件

generator/admission/Session replay/current writer、public literal commands/stepsとscreenshots state identity、Gallery/public asset差分、references read-only比較。geometry変更意図がある場合だけ正式reference更新手順。

共通作業規約のverification commandsからscopeに必要なものを実行し、必要な追加targetを結果に記録してください。
失敗を修正し、未測定をpassと報告しないでください。unchanged evidenceは条件一致時だけ再利用してください。

## 保存・コミット・プッシュ

current artifacts、実行可能docs、再生成/視認evidence、compatibility cleanup、S07 resultをcommit/push。public smoke diagramで代用しない。

`results/S07_RESULT.md`にcommands/results、source/input/environment SHA、owners/paths、
authority/base references、acceptance ID、remaining boundaryと次sessionの具体的入口を保存してください。
このsessionの実装・tests・docs・結果だけを一commitにし、終了時に必ず同名remoteへpushしてください。
English commit title: `Regenerate current source bindings and document input workflows`。

```bash
git branch --show-current
git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}'
git push origin HEAD:refs/heads/fix/issue-597-input-session-20260926
git ls-remote --heads origin refs/heads/fix/issue-597-input-session-20260926
git rev-parse HEAD
```

stage/commit前チェックとnon-fast-forward回復は共通作業規約に従います。force pushはしません。
local/remote SHA一致を確認してからbranch/commit/checks/残件を報告してください。
このpromptの次sessionまで連続実行せず、S07の範囲で終了してください。
