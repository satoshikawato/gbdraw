# INSTRUCTION PROMPT — S03: Discovery state and transform disclosure

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

native automatic discoveryを保護し、deferred/loading/errorとsingle crop applicabilityを正しく見せる。

## 開始条件

S02 resultとdiscovery Product authorityがmerged。source collectionがfresh/restoreで一貫している。

## 所有範囲

circular-sources/record-discovery/record-displayとpresentation owner、index.html/setup/watchers、focused browser/state tests、S03 result。new parser/schema/rendererは追加しない。

他sessionが変更したcodeをrevertせず、最新の同名 branchに合わせてscope内を実装してください。
SOLID/KISS/DRY/YAGNIをstate/owners/paths/compatibility/deliveryにも適用し、
不要なframework、重複owner、並行 fallback、将来用 capabilityを追加しないでください。

## 実行手順

1. idle/deferred/loading/ready/errorをsource operationと一致させる。valid native uploadは自動探索、rare formatは既存helper。completionはUID/version/mode/typeを再確認。
2. saved previewはsourceを再探索せずPython Worker0。Records not inspected/Inspect source records、実際のloading、Retry/Replace/Removeを示す。
3. artifactとdraftのsource同一性を確認せずcatalog metadataを流用しない。新source交換は自動探索に戻す。
4. Single-record crop, orientation and titlesと外側の適用不可理由/一件選択導線。applicableになるupload/selector transitionで自動展開し、manual close後は無関係更新で開かない。
5. focus/scroll anchorとkeyboard/390pxを保ち、grouping/先頭record/crop条件を自動変更しない。disclosureをpersisted state/History/Generate triggerにしない。
6. enabled renderingとaction resolution、availability reconciliationを同じowner/predicateに収束。

## 検証と完了条件

D-01〜D-04。record-display-discovery/circular-record-presentation、one/two/duplicates/GFF incomplete/invalid/rapid mutations/History、native upload の Python Worker 0 / helper settlement。upload直後をassertしmanual refreshで代用しない。

共通作業規約のverification commandsからscopeに必要なものを実行し、必要な追加targetを結果に記録してください。
失敗を修正し、未測定をpassと報告しないでください。unchanged evidenceは条件一致時だけ再利用してください。

## 保存・コミット・プッシュ

UI/state testsとS03 resultをcommit/push。完成artifactを読みやすいscaleで確認し、public screenshotsの更新対象をS07に渡す。

`results/S03_RESULT.md`にcommands/results、source/input/environment SHA、owners/paths、
authority/base references、acceptance ID、remaining boundaryと次sessionの具体的入口を保存してください。
このsessionの実装・tests・docs・結果だけを一commitにし、終了時に必ず同名remoteへpushしてください。
English commit title: `Expose applicable Circular transforms and truthful discovery states`。

```bash
git branch --show-current
git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}'
git push origin HEAD:refs/heads/fix/issue-597-input-session-20260926
git ls-remote --heads origin refs/heads/fix/issue-597-input-session-20260926
git rev-parse HEAD
```

stage/commit前チェックとnon-fast-forward回復は共通作業規約に従います。force pushはしません。
local/remote SHA一致を確認してからbranch/commit/checks/残件を報告してください。
このpromptの次sessionまで連続実行せず、S03の範囲で終了してください。
