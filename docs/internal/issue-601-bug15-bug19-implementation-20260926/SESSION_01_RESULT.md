# S01 — Error disclosure and regex recovery authority candidate

実施日: 2026-09-27（JST）。状態: **AUTHORITY_PENDING**。
二つの承認済みA receiptを既存の静的Product contractへ記録し、検証済み候補を同名remoteへpushした。
この結果は非規範のhandoffであり、Product判断、owner移管、runtimeの受入完了を新たに宣言しない。
S01固有のruntime変更なし。後述の最新dev同期だけを取り込み、S02以降は未開始。PR作成・dev merge・main・release/deployは未実行。
本書の最終commit SHAとremote一致はcommit後のhandoffで報告し、自己参照のためのamendはしない。

## Start and candidate

専用clone: `/tmp/gbdraw-issue601-s01-JnPiSH/repo`。
共有作業ツリー、S00 clone、他sessionのtree/index/未コミット変更/環境は操作していない。
事前に読んだ指針は最新cloneとbyte一致を確認した。cloneのAGENTS.mdとCLAUDE.mdも読了した。

| 項目 | 確認値 |
| --- | --- |
| 開始branch / upstream | `fix/issue-601-bug15-bug19` / `origin/fix/issue-601-bug15-bug19` |
| 開始SHA / 取得remote SHA | `4e3d33391e68ccb7905452641dd0aea45932e7a9`、clean、一致 |
| 開始時origin/dev | `4e3d33391e68ccb7905452641dd0aea45932e7a9` |
| 終了前に取得した最新origin/dev | `34c5104be196bda9d4a034ef4374128edb837982`（PR #612の#598 runtime統合） |
| S00資料commit | `202fe9de554aaa70dc731deb80bf032f26d80061` |
| S00 dev統合commit | `4e3d33391e68ccb7905452641dd0aea45932e7a9`、origin/devのancestor確認成功 |
| S00結果 | [SESSION_00_RESULT.md](./SESSION_00_RESULT.md)、開始treeに存在。先行resultはこれ一件 |
| authority branch / upstream | `product/issue-601-errors-regex-authority-20260926` / `origin/product/issue-601-errors-regex-authority-20260926` |
| authority本体commit | `84db467899d710cd4cb201b62f17ad4005dbd0b9`（指定title/summary） |
| authority最新候補HEAD / 確認remote SHA | `c68d1a00ed561f56caf2a9aead52d1043c35c4cc`、最新devを通常merge、一致、push後clean |
| 候補の変更file | `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md`のみ |
| dev統合 | **AUTHORITY_PENDING**。候補SHAのorigin/dev ancestor確認はexit 1 |
| 実装branchへ戻った直後の基点 | fetch / pull --ff-only後は開始SHAと同一、clean |
| 最終handoffの基点 | 最新dev `34c5104b`。未公開のlocal結果commitだけをこのbaseへ載せ直し、公開remoteのancestryを保持。最終差分は結果一fileだけ |

[Authority候補の差分](https://github.com/satoshikawato/gbdraw/commit/84db467899d710cd4cb201b62f17ad4005dbd0b9)はOIPC一ファイルだけ。
開始時devのrevision **21**とPD-OI-001–045を実際に確認してから、revision **22**、PD-OI-046/047を採番した。
終了前fetchで#598統合を検出し、authorityには通常mergeだけを追加した。最新devのOIPC、map/store、
checker/CIは旧baseとbyte同一。最新baseからtoolsを再抽出し、candidate差分がOIPC一fileだけであること、
18項目・既存全条件・digestとGateを再検証した。runtime/tests/generatedは最新devと同一。
番号を事前予約していない。base allowlistは`satoshikawato`、BD storeは空。
mapの二concern・四hard complete ruleに今回の二concernは未登録であり、map/storeを変更していない。

## Receipt fidelity and existing authority

| PD / concern / scenario | 選択 | 本文・JSON照合 |
| --- | --- | --- |
| `PD-OI-046` / `web.errors.diagnostic-disclosure` / 1 | `A / GUIDANCE_WITH_BOUNDED_DIAGNOSTICS` | 9/9項目完全一致 |
| `PD-OI-047` / `web.rules.rejected-pattern-edit-recovery` / 1 | `A / KEEP_REJECTED_PATTERN_DRAFT` | 9/9項目完全一致 |

正本は[診断公開receipt](./decisions/DECISION_01_ERROR_DISCLOSURE.md)と
[Color field回復receipt](./decisions/DECISION_02_REGEX_EDIT_RECOVERY.md)。
各正本のS00資料commit、開始remote、候補baseのbytesは同一。
候補は本文とJSONを両方そのまま保存し、九項目の全文を文字列単位で照合した。
scenarioRevisionだけは元JSONと同じ整数1。理由、維持条件、退役範囲、残余リスク、owner/dateを翻訳・要約・追加していない。

Receipt SHA-256（UTF-8、末尾改行除外）は正本と候補で一致:

- PD-OI-046: `8617888cb1838521f78640db4653e2dcf88dd427cd466bbbc33d59efacaceae8`
- PD-OI-047: `9cc66bc30f97cc995b9b0002b094ba74cfa00a8e77b163c33074689f0e73b496`

候補OIPC SHA-256: `5451db7da8e4de970f05afd5d19b003b7de77d52eb84c05d71ff85a533d4b1fb`。
既存45件のPD本文、Interpretation/lifecycle、OIPC-C01–C08、Acceptance catalog、
Residual-risk boundaryと既存evidence sectionはすべてbyte同一。変更はmetadataと二record追加だけ。

| 独立要求 | 既存正式authorityのcoverage / 今回の差分 |
| --- | --- |
| 失敗・cancel・stale・superseded時の成功Result/request保護、draftと修正retry | OIPC-C07、PD-OI-016、mapped saved-session continuityを維持。これだけではDetails/Copy/draft回復を全充足しない |
| Alignのlocal draft・underlying error・retry、record-owned direction、atomic History | PD-OI-027/029/031/034、PD-OI-035/039の独立条件を維持。Copy公開範囲の代替authorityではない |
| bounded code/operation/stage/context/副因、keyboard Details、表示中情報限定Copy、clipboard不可時の選択コピー、private/raw非公開、初回Resultとの区別 | 全条件を選択する既存正式recordなし。PD-OI-046に完全receiptを追加 |
| Python rule意味、一Worker・一preparation・prepared reuse・atomic live edit | Web既存契約とPython ownerの保存条件を維持。新しいregex evaluatorへの選択はしていない |
| 拒否pattern field、Not applied、Retry/Revert、close/reopen・一時mode保持、row/revision、対象History/row/document/session/reset解放、成功editだけHistory、非永続、Save/Generate/Export区別 | 全条件を選択する既存正式recordなし。PD-OI-047に既存Color pattern限定の完全receiptを追加 |

要求はANDで保持する。高位option ID、既存正常系tests、current runtimeだけで全coverageとはしない。
PD-OI-037/038/039/044/045、compact canvas、Session排他、科学的出力などの独立条件は変更していない。

## Same-semantic export candidate and writers

最新export branchは`fix/issue-601-export-output-20260926`、SHA
`43836eb77798924bc826d0064d3cc5e219379d91`。
[MASTER_PLAN](https://github.com/satoshikawato/gbdraw/blob/43836eb77798924bc826d0064d3cc5e219379d91/docs/internal/issue-601-export-output-20260926/MASTER_PLAN.md)、
APPROVED_DECISIONS、DECISION_PACK_02、SESSION_01_AUTHORITY、SESSION_03_ERRORS_AND_REGEXを再読した。
計画にはBUG-15/19の担当が残る。session resultはなく、devからの三点diffにruntime/tests/tools/CI変更はない。
**BUG-15/19のowner移管は未確認**であり、計画停止・移管・receipt退役を推測していない。

別candidateは`web.errors.user-facing-diagnostic-disclosure` / scenario 1 /
`A / ACTIONABLE_SUMMARY_WITH_SAFE_DETAILS`。

| 独立要求 | export candidateと今回receiptの比較 |
| --- | --- |
| 既知cause・修正/次action、actual stage、unknown非捏造、有界許可診断、privacy、Align draft/retry、Result/request/direction/History、cancel区別 | 共通して要求。選択IDの差だけを矛盾扱いしない |
| optional Details、keyboard、390 px、local-only、native例外捕捉 | exportはinitial collapseと常時summary、390 px/native/local-onlyも明示。今回receiptと実質矛盾はなく、既存baseの独立privacy/accessibility/native条件も残る |
| Copy diagnosticsの表示情報限定、手動操作、clipboard不可の手動選択、operation/副因、初回と旧Result区別 | 今回receiptに明記。export receiptだけでこれらの全coverageを証明できない |
| 既存Color patternのrejected draft、Retry/Revert、lifecycle、非永続とaccepted rule分離 | 今回の独立receipt。exportのPython rule維持だけでは代替しない |
| consoleのraw cleanup/自由ログ排除、既知修正情報をunknownへ落とさない | 共通して保存する意味。自由rawを診断contextへ移す許可ではない |

取得した全45remote refsのOIPCを検索し、開始時点で三つの診断/回復concernを持つ正式record・remote候補はなかった。
予定authority名とexport側予定名`product/issue-601-decisions-20260926`も作成前にremote不在。
今回の明示依頼に従い、二receiptの候補だけを単独writerとして記録・pushした。
別keyの同義recordを追加していない。key remap、export receiptの退役・supersession、owner移管はしていない。
同義候補の担当間調整は**未完了**。今回のpushは調整完了やdev統合許可を意味しない。
今後export候補を別のactive contractとして追加する収束・公開は停止し、一つの正式authorityで全独立要求を満たす扱いを明示する必要がある。
remote branch不在・result不在は分散writer leaseや担当停止の証拠ではない。
実際の競合candidate/同名writerが現れた場合は再確認し、上書き・force pushで解決しない。

## Shared runtime handoff

| 最新取得branch / SHA | 確認結果 |
| --- | --- |
| `fix/issue-598-alignment-direction-reset-20260926` / `400c746bf2e407f67078a0f691ab1704ad838973` | S02結果を読了。S00時点のS01より進み、committed candidate、direction/Reset、Session、normalizer/UI等のruntimeがある。dev ancestor exit 1 |
| `integrate/issue-598-runtime-dev-20260926` / `2afa57fc661b7835ddd1bbc557933510da2af1ce` | 初回取得時PR #612はOPEN、dev ancestor exit 1。終了前fetchではPR #612が`34c5104b`へ統合済み、integration HEADのdev ancestor exit 0 |
| `fix/issue-602-linear-live-edit-20260926` / `94faf8a98eddf823f5ebfed859d324e258daa689` | S05結果を再読。derived status・compact Editorまで。S06/S07未完了、runtimeのdev ancestor exit 1 |

#598の共有範囲はWorker/client、run-analysis、Align、normalizer、app-setup/index、
record controls、request/Session/Historyとnative planner等へ拡大している。
S02のtraceback cause保持とretry成功時error消去を、後続bounded disclosureへの移行で失わない。
最新devの[DEV_INTEGRATION.md](../issue-598-alignment-direction-reset-20260926/evidence/DEV_INTEGRATION.md)も読了。
これは別担当のaccepted runtime snapshot統合であり、元#598 writer session完了や#602 S06最終受入を宣言しない。
統合evidenceの711 Node・287 affected session・26 browser等を本S01の合格へ流用しない。
#602の共有範囲はindex/app-setup/status/svg-styles/config/active-config/session-request。
S05の既知Legend binding failure、metadata-free Session staging問題を今回直したとはしない。
両branchのlocal合格を全dev staging合格や本S01の受入へ流用しない。
#598の統合基点は`34c5104b`として確定した。#602はなお未統合。後続はこの差分と最新#602結果・
owner移管・作業順を確認し、既存ownerに残る差分だけを実装する。

## Verification: Gate and Review

trusted toolsは確認した最新devの`git archive`からsession専用directoryへ抽出した。
checker、detectors、rules、Product map/BD store、CIをcandidateで変更・弱化していない。
候補を許可するためにcandidate authorityをruntime評価へ用いていない。

| 検証 | 結果 / 証拠 |
| --- | --- |
| 原文/JSON/候補本文の18項目、digest、allowlist、metadata、既存record/条件保持 | PASS。初回`authority-verification.log`、最新baseにも`authority-verification-latest.log`で再検証 |
| OIPC参照先 | local links 22件resolve成功。二receiptのpinned source commit/path実在・内容一致 |
| authority scope・whitespace | PASS。OIPC一ファイルのみ、104 additions / 2 deletions。production/tests/tools/CI/generated差分なし |
| latest-base actual diff CI plan | `policy-documentation` / `DOCUMENTATION_ONLY_PR`。requiredJobs=`web-change-budget`、inheritedEvidence=null。`authority-ci-plan.log` / `authority-ci-plan-latest.log` |
| authority working-tree / actual-head Web gate | **Gate PASS / Review REQUIRED**。四hard rules CONFORMING、candidate authority separation VALID。`authority-gate-precommit.log` / `authority-gate-postcommit.log`、最新base/headの`authority-gate-latest.log` |
| handoff documentation gate | **最新baseもPASS、189 passed / 6107 deselected、83.48s**（`handoff-recipes-latest.log`）。旧base189 passed / 6090 deselected、82.39sは別記録 |
| handoff scope・参照・whitespace・trusted-base Web gate | **PASS / Review CLEAR**。結果一ファイルのみ、local参照・whitespace成功。旧base`handoff-gate-precommit.log`と最新base`handoff-gate-latest-precommit.log`。actual-head gateはcommit後に実行し、最終handoffで報告 |

handoffのstaged diff分類は`documentation`、requiredJobs=`recipes-standard`。authorityとhandoffを同じdiffに混ぜていない。
Review REQUIREDはauthority変更に対する人間のレビュー要求であり、Gate成功やpushをreview承認と扱わない。
S01では新しいruntime、owner/path、互換reader、privilege/dependencyを導入していない。
production owner/pathとOE/PE/CBは変化しない。完全例外packetの条件は発生しない。
production、tests、docs、generatedのdiffを個別にreviewした。
recipe検査が自分のclone内の`gbdraw/bin/linux-x86_64/losat`だけを100755へchmodしたため、
byte変更がないことを確認してtracked 100644へ戻した。最終runtime差分は空。

実行command（専用clone root、logs/cache/outputはsession専用/tmp）:

```bash
git fetch origin
git switch --track -c fix/issue-601-bug15-bug19 origin/fix/issue-601-bug15-bug19
git pull --ff-only
git rev-parse HEAD origin/dev
git rev-parse --abbrev-ref --symbolic-full-name '@{upstream}'
git merge-base --is-ancestor 4e3d33391e68ccb7905452641dd0aea45932e7a9 origin/dev
git switch --no-track -c product/issue-601-errors-regex-authority-20260926 origin/dev
git archive 4e3d33391e68ccb7905452641dd0aea45932e7a9 tools | tar -x -C /tmp/gbdraw-issue601-s01-JnPiSH/trusted-tools
python /tmp/gbdraw-issue601-s01-JnPiSH/verify-authority.py
node /tmp/gbdraw-issue601-s01-JnPiSH/trusted-tools/tools/check-web-change-budget.mjs --base 4e3d33391e68ccb7905452641dd0aea45932e7a9 --head 84db467899d710cd4cb201b62f17ad4005dbd0b9
# 終了前に最新devを検出し、authorityを通常同期して再検証:
git merge --no-edit origin/dev
node /tmp/gbdraw-issue601-s01-JnPiSH/trusted-tools-latest/tools/check-web-change-budget.mjs --base 34c5104be196bda9d4a034ef4374128edb837982 --head c68d1a00ed561f56caf2a9aead52d1043c35c4cc
git push -u origin HEAD:refs/heads/product/issue-601-errors-regex-authority-20260926
git switch fix/issue-601-bug15-bug19
git fetch origin
git pull --ff-only
PYTHONDONTWRITEBYTECODE=1 /tmp/gbdraw-issue601-s01-JnPiSH/venv/bin/python -m pytest tests/ -m 'recipe and not slow' --durations=30 --basetemp /tmp/gbdraw-issue601-s01-JnPiSH/pytest -o cache_dir=/tmp/gbdraw-issue601-s01-JnPiSH/pytest-cache
# 未公開local handoffだけを最新devへ載せ直し、source変更に対して再検証:
git rebase origin/dev
PYTHONDONTWRITEBYTECODE=1 /tmp/gbdraw-issue601-s01-JnPiSH/venv/bin/python -m pytest tests/ -m 'recipe and not slow' --durations=30 --basetemp /tmp/gbdraw-issue601-s01-JnPiSH/pytest-latest -o cache_dir=/tmp/gbdraw-issue601-s01-JnPiSH/pytest-cache-latest
node /tmp/gbdraw-issue601-s01-JnPiSH/trusted-tools-latest/tools/check-web-change-budget.mjs --base 34c5104be196bda9d4a034ef4374128edb837982 --head HEAD
git diff --check
git diff --cached --check
git push origin HEAD:refs/heads/fix/issue-601-bug15-bug19
git ls-remote --heads origin product/issue-601-errors-regex-authority-20260926 fix/issue-601-bug15-bug19 dev
```

actual-head CI分類はtrusted `tools/ci-impact.mjs plan`に、profile=pr / event=pull_request / repository=satoshikawato/gbdraw /
actual base/head/workflow SHA / architectureChange=falseを渡して取得した。PRそのものは作成していない。
環境: Node 26.8.2、Python 3.13.3、pytest 9.0.2、session専用venvは
`/tmp/gbdraw-issue601-s01-JnPiSH/venv`（system packagesをread-only継承、共有editable installなし）。
gbdraw importは専用clone。clean外部ディレクトリのrecipeもtestのPYTHONPATHが当該REPO_ROOTを指定する。
Python 3.13でのlocal recipeをCI Python 3.11 matrixの実行証拠と表現しない。

41 Node・9 Python・2 ChromiumはS00以前のruntime baselineとしてだけ保持する。
開始baseではsource/input tree同一だったが、終了前の#598統合で変わった。必要なrecipe gateを再実行し、
既存runtime baselineを現在dev・今回authority・将来runtimeの合格証明にしていない。
今回のscopeに不要な全suite、browser、wheel/cache-bust生成は繰り返していない。

## Current dev staging and next boundary

開始時dev exact SHA `4e3d33391e68ccb7905452641dd0aea45932e7a9`のremote実状態を一度確認した。
[Tests run 36250628664](https://github.com/satoshikawato/gbdraw/actions/runs/36250628664)はfailure:
CI impact planが`DOCUMENTATION_BASE_EVIDENCE_UNAVAILABLE` / `RUN_NOT_SUCCESSFUL`、
Dev staging / gateが`MISSING_ENVIRONMENT`（CI_IMPACT_PLAN_JSON欠落）、leaf jobsはskip。
同SHAの[Gallery run 36250628605](https://github.com/satoshikawato/gbdraw/actions/runs/36250628605)と
CodeQLはsuccessだが、全dev staging成功の代用にしない。
過去SHA/PR #611のchecksも代用にしていない。workflowの再実行やCI修正は本S01の範囲外。
最新dev `34c5104be196bda9d4a034ef4374128edb837982`も一度読み取り確認した。
[Tests run 36252205658](https://github.com/satoshikawato/gbdraw/actions/runs/36252205658)は確認時点で**in_progress**、
同SHAのGallery run 36252205641とCodeQLはsuccess。全dev staging成功はまだ未確認。
S01に不要なstagingの完了待ちや再dispatchは行わず、実状態の確認を後続に引き継ぐ。
五分未満のCI pollingはしていない。

| S02開始条件 | 状態 |
| --- | --- |
| S00公開・二完全receipt・最新authority inventory | 充足 |
| 二receiptのauthority候補保存・同名remote一致 | 充足。PD-OI-046/047候補をpush済み |
| 正式authorityのdev統合・ancestor確認・実装branchへの取り込み | **未充足 / AUTHORITY_PENDING**。PR作成・dev merge許可なし |
| export側BUG-15/19 owner移管・同義candidateの扱い・一writer確定 | **未充足**。対象scope、引継ぎSHA、担当停止/移管と一つの正式authorityへの扱いを明示する必要あり |
| #598共有runtimeの統合基点 | 充足。PR #612 / dev `34c5104b`とDEV_INTEGRATIONを確認。別writer sessionや#602最終受入は別状態 |
| #602共有runtimeの受入SHAと作業順 | **未充足**。#602 runtimeはdev未統合、owner/作業順と最新handoffが必要 |
| 原監査pattern・画面・公開build | 未確認。現在Color/Label成功や本authority保存をBUG-19 closure証拠としない |

**S02は開始不可。** 次の許可境界は、検証済みauthority候補への人間のレビューとPR作成・dev統合。
同義候補の独立要求を削除せず、一つの正式authorityにどう接続するかとBUG-15/19 owner移管を明示する。
共有runtimeの統合・受入基点が確定してから、必要な最新baseチェックとancestor確認を行う。
今回の候補authority自身でruntimeを自己承認しない。

Authority commit title: `Record approved error disclosure and regex recovery behavior`

Summary: Preserve the two complete product receipts in the existing static contract.

Handoff commit title: `Record authority integration status for issue 601`

Summary: Record the authority candidate and the prerequisite for runtime work.
