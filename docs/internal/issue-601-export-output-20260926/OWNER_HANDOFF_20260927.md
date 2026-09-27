# BUG-15 / BUG-19 ownership handoff

実施日: 2026-09-27 (JST)。状態: **TRANSFERRED**。
利用者が「BUG-15/19のowner移管と#602共有runtimeの引継ぎは未完了」を選択し、
「完了させてください。」と明示依頼したことに基づく実行担当の移管記録。
Product outcome、receiptの退役・supersession、concern key変更を新たに選択しない。

## Current execution owners

| Scope | 唯一の実行担当 / branch | 順序 |
| --- | --- | --- |
| BUG-15: typed cause、Python render/helper adapter、Worker/client、normalizer、Generate/Align caller、Details/Copy/recovery | `fix/issue-601-bug15-bug19` | 同計画 S02 → S03 → S05 |
| BUG-19: Color/Label Python・Search JS方言案内、既存Color pattern draft、Retry/Revert/lifecycle | `fix/issue-601-bug15-bug19` | 同計画 S03 → S04 → S05 |
| BUG-07: PDF font/manifest/loader、snapshot、配布、PDF固有producer、PDF受入 | `fix/issue-601-export-output-20260926` | この計画 S00 → PDF authority S01 → S02 → PDF受入 S04 |

このbranchのS03によるBUG-15/19実装は移管により停止する。
S01では診断公開の第二active authorityを追加しない。S04はerror/regexの実装や修正を
再所有せず、移管先のpush済み受入を参照してPDFとの接続を確認する。
PDF failure producerはPDF担当に残るが、共通adapter/transport/normalizer/UIは移管先が所有する。
共通境界をPDF担当が別実装・同時編集しない。必要なfinite code/actionは移管先へ引き渡す。

移管先計画: [MASTER_PLAN](../issue-601-bug15-bug19-implementation-20260926/MASTER_PLAN.md)。
正式authorityはdev統合済みOIPC revision 22、PD-OI-046/047、PR #615 /
`744be7a5943d4a247d027369629058898aa3f33e`。

## Receipts and independent requirements

[Decision Pack 02](./DECISION_PACK_02_ERROR_DISCLOSURE.md)と
[APPROVED_DECISIONS](./APPROVED_DECISIONS.md)の元本文・JSONはbyte不変。
`web.errors.user-facing-diagnostic-disclosure` / scenario 1 /
`A / ACTIONABLE_SUMMARY_WITH_SAFE_DETAILS`を退役・key remap・別keyへ書換えしていない。
このbranchの未公開診断authority候補作業を移管先へ引き渡し、別のactive recordを追加しない。
実装は一つの正式authorityで進め、以下を独立した保存・受入条件として保持する。

| 独立条件 | 移管後の扱い |
| --- | --- |
| known causeと修正/next action、actual operation/stage、unknown非捏造、bounded allowlistとprivacy、cleanup主因保持 | 移管先の同じproducer → adapter → Worker/client → normalizerで保持 |
| Detailsは初期collapseの任意表示、summary常時可視、keyboard/390px、native CLI/API、Python rule、local-only | 元承認の条件と既存正式Web契約を維持。高位option IDだけでcoverageとしない |
| Result/request/draft/direction/History、retry、cancel/stale/superseded、transient診断 | 移管先の失敗・回復受入へ引継ぐ |
| manual Copyは表示中のbounded情報のみ、clipboard不可時の手動選択、初回と旧Result区別 | PD-OI-046の追加独立条件を削除しない |
| Color field rejected draft、Not applied、Retry/Revert、accepted rule分離、非永続、row/document/History lifecycle | PD-OI-047を独立して実装・受入。PDF担当には移さない |

この文書は非規範の作業担当記録で、追加のruntime authorityやCI decision storeではない。
新たな製品矛盾が判明した場合は影響するoutcomeだけをProduct Decision Ownerへ提示する。

## Publication and next writer

開始remoteは`43836eb77798924bc826d0064d3cc5e219379d91`、計画のみでBUG-15/19 runtime差分なし。
最新dev `c922fc38aac78da9be83342c09ac0164ecef6ff6`を通常mergeし、正式authorityを取り込んだ。
今回の変更は本書とMASTER_PLAN/S01/S03/S04の担当指示のみ。元receipt、runtime、tests、
checker/CI/authorityは編集しない。次の実行者は最新remoteの本書を先に読む。
分散processの終了や恒久leaseを推測した記録ではない。担当指示は今回の明示依頼により変更し、
各branchのpush前にremote更新を再確認する。同名branchへの同時writerは禁止。

#602受入SHAは`9a14db3af8fafd6f8374cfe2a8197a8927cf43f7` (S07完了)。
共有runtimeの取り込み・file順序は移管先のSESSION_01_HANDOFF_RESULT.mdに記録する。
PDF branchは共通fileの編集前にその結果と最新devを取得し、同じownerへ合わせる。

## Verification

最新dev由来checker、base `c922fc38aac78da9be83342c09ac0164ecef6ff6`で
移管文書のworking-tree Gate PASS / Review CLEAR。五つの文書だけを新規編集し、
参照・whitespace、元承認三fileのbyte保持、production/tests/tools/CIのdev同一性を確認。
Documentation recipe gateは先行の189 passed / 6107 deselected (94.19s)を再利用する。
検証tree `9d967f1c`からの追加差分はinternal文書と独立Node architecture testsだけで、
Python/recipe code・fixtures・公開recipe入力・環境・受入条件は不変。
Log: `/tmp/gbdraw-issue601-authority-integration-Xp5su7/handoff-recipes.log`。
Actual-head gate/CI分類、commit/remote SHAはcommit後のhandoffで確認し、自己参照amendしない。
