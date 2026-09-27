# S01 authority integration and S02 prerequisites

実施日: 2026-09-27（JST）。状態: **AUTHORITY_INTEGRATED / S02_BLOCKED**。
[S01の候補作成結果](./SESSION_01_RESULT.md)は当時の記録として保持する。
今回の利用者の「承認します。devに統合してください。」により、検証済みauthority候補の
PR作成・通常dev mergeを完了し、非規範handoffを記録する。BUG-15/19のowner移管、
export receiptの退役・key remap、#602 runtime統合はこの承認から推測しない。
S02以降は開始していない。

## Exact candidate and integration

| 項目 | 確認値 |
| --- | --- |
| 専用clone | `/tmp/gbdraw-issue601-authority-integration-Xp5su7/repo` |
| 開始authority HEAD / remote | `c68d1a00ed561f56caf2a9aead52d1043c35c4cc`、clean、一致 |
| Authority本体commit | `84db467899d710cd4cb201b62f17ad4005dbd0b9` |
| 開始handoff remote | `b5014ed1d92fb492d878fd38c266cf84cdecc726` |
| 開始origin/dev / trusted base | `9d967f1c72f730b205420d62d383b73127c2a9b1` |
| 最新dev同期済みauthority HEAD | `c1fc4680c4082541f3a18ec586af9455808ecdb5`、通常merge、承認済みOIPC本文は不変 |
| Authority PR | [#615](https://github.com/satoshikawato/gbdraw/pull/615)、base `dev` |
| 正式dev統合commit | `744be7a5943d4a247d027369629058898aa3f33e`、candidate ancestor・OIPC byte一致を確認 |
| 実装branch同期 | `7e3bbf07eb6c308f00dccfcd5661ee50f6811646`、devの通常merge、統合commitのancestor確認成功 |

専用cloneは既存cloneのimmutable objectsを読み取りコピーした後にdissociateし、
`.git/objects/info/alternates`不在を確認した。共有tree/index/環境は操作していない。
AGENTS/CLAUDE/Web CLAUDE/PR skillは読了した共有側原本とSHA-256一致。
dev/mainへ直接pushせず、authority PRを保護された通常経路で統合した。

## Receipt fidelity and scope

OIPC revision **22**、PD-OI-046/047の本文・JSONは二つの原receiptと全18field一致。
正本はS00 commit `202fe9de554aaa70dc731deb80bf032f26d80061`の
[診断公開](./decisions/DECISION_01_ERROR_DISCLOSURE.md)と
[Color pattern回復](./decisions/DECISION_02_REGEX_EDIT_RECOVERY.md)。
理由・維持条件・退役範囲・残余リスク・owner/dateを追加・要約・翻訳していない。

- `PD-OI-046` / `web.errors.diagnostic-disclosure` / scenario 1 /
  `A / GUIDANCE_WITH_BOUNDED_DIAGNOSTICS`:
  receipt SHA-256 `8617888cb1838521f78640db4653e2dcf88dd427cd466bbbc33d59efacaceae8`。
- `PD-OI-047` / `web.rules.rejected-pattern-edit-recovery` / scenario 1 /
  `A / KEEP_REJECTED_PATTERN_DRAFT`:
  receipt SHA-256 `9cc66bc30f97cc995b9b0002b094ba74cfa00a8e77b163c33074689f0e73b496`。
- OIPC SHA-256 `5451db7da8e4de970f05afd5d19b003b7de77d52eb84c05d71ff85a533d4b1fb`。

既存45件のPD本文、Interpretation/lifecycle、cross-surface条件、Acceptance catalog、
Residual-risk/evidence sectionsをbyte比較し、不変を確認。22local linksもresolve成功。
authority追加差分はOIPC一fileのみ、104 additions / 2 deletions、whitespace PASS。
production/tests/checker/CI/Product map/BD store/generated差分は空。
Production、tests、docs、generatedのdiffを別々にreviewした。
production owner/path/compatibilityは不変、OE/PE/CB変化なし。architecture exceptionなし。
候補のauthorityで候補runtimeを自己承認していない。

## Verification and review

| 検証 | 結果 |
| --- | --- |
| Latest-base authority-only local checker | Gate PASS / Review REQUIRED、四hard rules CONFORMING、candidate separation VALID |
| Authority CI plan | policy-documentation / DOCUMENTATION_ONLY_PR、requiredJobs=web-change-budget、inheritedEvidence=null |
| PR #615 exact head必須checks | Web base policy (trusted base) SUCCESS ([run 36282080825](https://github.com/satoshikawato/gbdraw/actions/runs/36282080825))、PR / gate SUCCESS ([Tests run 36282036591](https://github.com/satoshikawato/gbdraw/actions/runs/36282036591)) |
| その他PR checks | CodeQL SUCCESS。対象外leaf jobsはCI plan通りskip。最初のtrusted runはbody編集によるcancel後、後続runがSUCCESS |
| Documentation recipes | 189 passed / 6107 deselected、94.19s、Node 26.8.2 / Python 3.13.3、専用venv・当該clone import確認 |
| 統合・同期 | devにcandidate ancestor、contract SHA-256一致、実装branchに統合commit ancestor |

Recipesはauthority同期済みtreeで実行した。正式dev統合とhandoff同期後の差分は非規範記録のみで、
production/tests/recipe入力は同一のため結果を再利用した。Python 3.11 CI実行と同一視しない。
recipe検査が自分のcloneのlosatだけをchmodしたため、byte不変を確認して100644へ戻した。
Dev stagingはexact SHA `744be7a5943d4a247d027369629058898aa3f33e`で初回確認した。
[Tests run 36282346513](https://github.com/satoshikawato/gbdraw/actions/runs/36282346513)、
[Gallery run 36282346409](https://github.com/satoshikawato/gbdraw/actions/runs/36282346409)、
dynamic run `Push on dev` (36282346369)も確認時点でin_progress。全dev staging合格は未確認。
PR成功や過去SHAの結果をこの合格へ代用せず、再dispatch・CI修正は行わない。
CI監視は五分以上の間隔。

Authority用trusted toolsは開始origin/dev `9d967f1c`、handoff用は統合後dev `744be7a5`の
`git archive`から、それぞれ専用/tmpへ抽出した。
ログ・検証scriptはsession rootにあり、repoへstageしない。
PR文面にはwrite-clear-pull-request skillを適用し、同じtitle/bodyを
`check-pr-language.mjs`で確認してからcreate/editした。
Gate PASSをReview承認と同一視せず、今回の明示承認に基づき通常mergeした。
Protectionを変更・迂回していない。必須statusは
`Web base policy (trusted base)`と`PR / gate`、required review countは0。

主な実行command（専用clone root。receipt検査は当時のauthority候補checkoutで実行）:

```sh
python /tmp/gbdraw-issue601-authority-integration-Xp5su7/verify-authority.py
node /tmp/gbdraw-issue601-authority-integration-Xp5su7/trusted-base/tools/check-web-change-budget.mjs --base 9d967f1c72f730b205420d62d383b73127c2a9b1 --head c1fc4680c4082541f3a18ec586af9455808ecdb5
PYTHONDONTWRITEBYTECODE=1 /tmp/gbdraw-issue601-authority-integration-Xp5su7/venv/bin/python -m pytest tests/ -m 'recipe and not slow' --durations=10 --basetemp /tmp/gbdraw-issue601-authority-integration-Xp5su7/pytest -o cache_dir=/tmp/gbdraw-issue601-authority-integration-Xp5su7/pytest-cache
```

## Remaining owner and shared-runtime boundaries

export branch `fix/issue-601-export-output-20260926`の取得HEADは
`43836eb77798924bc826d0064d3cc5e219379d91`。計画のみでruntime差分・session結果なし。
SESSION_03_ERRORS_AND_REGEXにproducer/adapter/Worker/client/normalizer/Align/UIの担当が残る。
この事実はwriter停止・lease・owner移管の証拠ではない。
同義candidate `web.errors.user-facing-diagnostic-disclosure` / scenario 1 /
`A / ACTIONABLE_SUMMARY_WITH_SAFE_DETAILS`は正式OIPCへ追加していない。

既知cause/next action、actual stage、bounded privacy、Result/request/draft/direction/History、
cancel、keyboard/390px、Python/native意味等は独立して保存する。
export候補だけではmanual Copy、表示情報限定、clipboard不可時の手動選択コピー、
初回と旧Result保持の区別、既存Color field draftのRetry/Revert/lifecycleを全充足しない。
異なるoption IDだけを矛盾とせず、実質的な矛盾は未検出。
今回の正式化はexport receiptの退役やconcern key変更ではない。
今後の同義active contract二重登録を避ける接続方法と担当移管はなお明示判断が必要。

#598はPR #612/#613/#614で最新devへ統合済み。
統合基点は`34c5104b` / `aab5ad43` / `9d967f1c`。
[DEV_INTEGRATION](../issue-598-alignment-direction-reset-20260926/evidence/DEV_INTEGRATION.md)、
[followup](../issue-598-alignment-direction-reset-20260926/evidence/DEV_INTEGRATION_FOLLOWUP.md)、
[S03統合結果](../issue-598-alignment-direction-reset-20260926/evidence/DEV_S03_INTEGRATION.md)を確認した。
成功retryのstale error消去、canonical label/source-bound Reset、unknown/no-candidate、
keyboard/focus、atomic failure/retryを後続で保持する。source結果とintegration検証は区別する。

#602の取得HEADは`0a3390e282ca8a3c6ae263da7185d39085fe01c1`。
[S06結果](https://github.com/satoshikawato/gbdraw/blob/0a3390e282ca8a3c6ae263da7185d39085fe01c1/docs/internal/issue-602-proposal-20260926/results/S06.md)
はS06実装・検証・通常push完了、S07未開始を記録する。
#598 PR #614まで同期しているが、#602のruntime自体はdev未統合。
共有範囲はindex/app-setup/status/svg-styles/config/active-config/session-request等。
このsessionでは#602を直接import・統合しない。受入SHA・owner・共有file作業順は未確定。
S06の42browser成功等を本authority受入や全dev staging合格へ流用しない。
Legend sanitizer retry制約とmetadata-free Session staging問題は未解決として保持する。

## S02 readiness

| 条件 | 状態 |
| --- | --- |
| 二receipt正式authorityのdev統合・ancestor・実装branch取り込み | **充足**、PR #615 / `744be7a5943d4a247d027369629058898aa3f33e` |
| S00・receipt・authority inventory | 充足。取得47remote refsのOIPCに競合recordなし。候補branchだけに二recordを確認してから統合 |
| export側BUG-15/19 owner移管・同義candidateの扱い・一writer | **未充足**、scope/担当停止・移管/引継ぎSHAの明示判断なし |
| #598共有runtimeの統合基点 | 充足、PR #612/#613/#614と各integration結果 |
| #602共有runtime受入SHA・共有file作業順 | **未充足**、S06結果を読了したがruntimeはdev未統合・作業順未確定 |
| 原監査pattern・画面・公開build | 未確認、原報告のclosure証拠にしない |

**S02開始不可。** BUG-15/19移管scope、引継ぎSHA、export担当停止/移管と一writer、
同義candidateの扱い、#602共有runtimeの受入基点と作業順を明示する必要がある。
原監査pattern・画面・公開buildも未確認。authority統合をBUG-19 closureと表現しない。
S01資料の過去状態を書き換えず、この記録を最新handoffとして読む。
main/release/tag/手動deploy/Issue閉鎖は実行していない。
Handoff branchの最新devに対する追加diffは、既存SESSION_01_RESULT.mdと今回の
SESSION_01_INTEGRATION_RESULT.mdの二つの非規範記録だけ。今回の新規編集は後者一file。
旧S01本文を保持し、runtime/tests/tools/CI/generatedは正式devと同一。
最終handoff commit SHA、actual-head Gate/CI分類、remote一致はcommit後の報告に記し、
自己参照のためにamendしない。

Proposed commit title: `Record issue 601 authority integration and runtime prerequisites`

Summary: Record the approved contract merge and the remaining ownership and shared-runtime prerequisites.
