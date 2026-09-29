# 実行状況

計画: [総合実装計画](00_MASTER_PLAN.md)。
作業ブランチ: fix/gui-feedback-remediation-20260928。
再利用 worktree: /mnt/c/Users/genom/GitHub/gbdraw/.worktrees/gui-feedback-remediation-20260928。
初期 dev: 57cef3ba47f4b7790a9145f2ce6988422a00710e。

2026-09-28: 計画書とセッション別プロンプトを作成。runtime・規約・checkerは未変更。
基準のソース調査と限定resolver確認は計画に記載した。S00の再現可能なbaseline採取は未実施。
Product契約の同時更新例外、計算重複防止契約も未実装であり、以下の順に実施する。

計画書作成時の検証: 10 Markdownの相対リンク・書式・branch指定・要求IDを確認済み。git diff --cached --checkは成功。node tools/check-web-change-budget.mjs --base origin/devはGate PASS / Review CLEAR。独立レビューで計算契約の範囲、旧Session、規約移行、BGC手順を確認した。runtime差分がないため製品テストはこの計画書commitでは実行していない。

2026-09-29: S00 を実施。runtime・規約・checker・期待出力は未変更。
dev と main を `git archive` で `$S00_BASELINE_DIR=/home/kawato/gbdraw-baselines/gui-remediation-20260928/` へ抽出し、各 snapshot 自身の source から browser wheel を準備した（tracked source 不変を manifest で確認）。
計測 harness は results/s00/ にあり、各 pass は直列に実行した（coverage dev/main、timing dev/main）。raw data は `$S00_BASELINE_DIR/perf/`、機能証拠は `$S00_BASELINE_DIR/evidence/`。
主要な結果と未確定事項は results/S00.md 第 1 節と第 11 節。commit title: "Record S00 baseline evidence and decision scope"。

2026-09-29: S01 を実施（docs-only）。Static Product Contract co-change 経路（base checker が実装した時点で有効）と CW-01〜06 を規範化した。検証は results/S01.md 第 5 節。commit title: "Allow reviewed Product Contract co-changes and define computation ownership"。

2026-09-29: S02 を実施（checker と fixture test のみ）。co-change 条件を companion 拒否と runtime＋guard 拒否の 2 か所へ同一条件で適用し、Review 理由を必須にした。governance test 221 passed、変異 4 件すべて検出。commit title: "Admit reviewed Product Contract co-changes in the Web checker"。

2026-09-29: PR #641（a9333e79）を S02 の後に merge（70bca139）。S03 を実施し、Pending status と表示専用 intent bookkeeping を除去した。Generate 1 回の label 表構築を 1 回にし、前処理中の Cancel と stale を判定するようにした。CW test を追加した。Node 1,115 passed。移行した Playwright は 52 passed、2 failed（#641 由来の既存失敗）。commit title: "Remove derived Pending status and guard duplicate table builds"。

2026-09-29: S04 を実施。S00 第 11.4 節の判断 1〜3 を実装した。Plot title・title font・定義 font を既存 mode-profile manager の per-mode 値にし、`linearRecordLayout`・`linearTypographyLinked` の欠落を fresh 既定にした。Circular/Linear request の projection が他 mode の値を発明しないようにし、Gallery publication は未使用 mode を fresh 既定で書く。Node 1,115 passed。browser は 124 passed、5 failed（追加 test の selector 1 件は修正済み、残り 4 件は PR #641 head でも失敗する既存失敗）。commit title: "Keep title and fonts per mode and restore fresh defaults for omitted layout"。

2026-09-29: S05 を実施。Similarity alignment の直接 RBH を、committed request の orthogroup resource から Worker adapter が取得するようにした（UI catalog の edge は常に空だった）。PR #641 の recipe-only 確定（`projectGeneratedProteinRecipe`）は、新規 LOSATP Generate 後の Align と record rotation を browser で失敗させていたため撤回した（Owner-delegated、第 2 節）。Review と plan inspector は record 名・accession と gene・protein ID を表示する。Node 1,116 passed、Python 561 passed、browser は既存の #641 由来の失敗を除き passed。commit title: "Resolve alignment RBH from committed orthogroup evidence and show biological names"。

2026-09-29: S06 を実施。Preview の Editor を上端から開き、検索（最大 39.5rem）と toolbar を Editor 幅を除いた残り幅に置いた。Align と Review を popup・drawer の同じ行へ移した。Align と Lock の常時説明を help-tip（hover・focus・tap、accessible description）へ、比較の「Current:」状態と要約を削除して `aria-pressed` にした。全 browser 実行で S04 の mode 復元が Circular 定義を live 更新する退行を見つけ、別 commit で修正した。commit title: "Open the Editor from the Preview top and move always-on help into help-tips"。

2026-09-29: S07 を実施。検索非表示時に Preview の canvas が 200 px に縮む grid 配置を固定行で直した（dev からの既存不具合で、docs・Gallery capture の中心合わせ失敗の原因）。#641 の bounded helper reply を test の Worker tracker が失敗として数えていた問題を直した。Gallery を owner tool で再生成し、S06 の影響を受けた tutorial media 7 枚と docs capture 11 scenario を再生成・目視した。最終候補で S00 harness、統合 journey、CW 変異 2 件、最終チェックを実施。比較切替の settlement 予算が 1 行不合格（間欠的な render frame、S03 と同頻度）で、Owner 判断を残す。commit title: "Keep the Preview canvas in the flexible row when the search is hidden"、"Count one Worker settlement per bounded helper reply in browser tests"、"Refresh Gallery artifacts and recapture comparison and Align tutorial media"、"Regenerate documentation captures for the pressed comparison state"、"Align tests with the pressed comparison state and the refreshed Vnig Session"、"Record S07 integration acceptance and handoff"。

2026-09-29: S07 の追加作業（Owner 指示）。比較切替の settlement 予算の未達を調べ、原因が大きな Result の後の V8 major GC であることをharness の診断（b39b083f）で確かめた。catalog の nested 値の共有（308c5402）で live heap を 37% 減らしたが、固定定義の 1 行は未達のままで、Owner が受容した。PR #641 の失敗を #641 の branch で修正して push し（a9333e79..d5347f69）、CI の pass の後に Owner が merge した（b1744142）。`origin/dev` を本 branch へ merge した（d1bd76f6）。結果は results/S07.md 第 12 節。

| Session | 状態 | 証拠・次の条件 |
| --- | --- | --- |
| S00 | 完了（証拠採取のみ） | commit d7fbe23b（証拠）、9b671081（Owner 判断）。[results/S00.md](results/S00.md)。開始 a351c01d、dev 57cef3ba・main 4556e04e の snapshot/wheel、性能 baseline、G03/G06/G07 再現、R01/R02 範囲。旧 Session 方針と receipt 文言は第 11.4 節で取得済み。commit SHA は次セッションで追記 |
| S01 | commit 済み・dev 未反映 | commit 19c9c939。[results/S01.md](results/S01.md)。docs-only。base checker で Gate PASS・Review REQUIRED。dev 取り込みは push/merge 承認待ち |
| S02 | commit 済み・dev 未反映 | commit 9e7a056f。[results/S02.md](results/S02.md)。checker＋fixture のみ。S01 未 merge のため S01 commit を base に局所検証。PR は S01 の dev 取り込み後 |
| S03 | commit 済み・dev 未反映 | commit 3bd092af。[results/S03.md](results/S03.md)。PR #641 を merge（70bca139）後に実装。Pending と表示専用 intent を除去し、PD-OI-037/049 を co-change で revision 2 にした。CW-01〜04 の自動検証と変異 3 件を検出。性能は全操作が提案予算内 |
| S04 | commit 済み・dev 未反映 | commit 2e0fdfba（follow-up 636d09f2）。[results/S04.md](results/S04.md)。S03（3bd092af）の上。判断 1〜3 を既存 mode-profiles・config・projection・publication の owner で実装。v41/v44 Vnig の受入と Save→fresh Load・拒否 Load を確認。#641 由来の既存 browser 失敗 4 件を特定（S07 へ）。Gallery 再生成は S07 |
| S05 | commit 済み・dev 未反映 | commit f04c599a。[results/S05.md](results/S05.md)。S04 の上。直接 edge を committed orthogroup resource から Worker adapter へ渡し、CLI と束縛 helper を共有。#641 の recipe-only 確定を撤回して新規 LOSATP 後の Align と record rotation を回復。Review・inspector の名前を共有 helper へ。実 Gallery と新規 LOSATP で livA→racM、parA→racL を確認 |
| S06 | commit 済み・dev 未反映 | commit f4a59295。[results/S06.md](results/S06.md)。S05 と S04 follow-up（636d09f2）の上。Contract revision 28（PD-OI-024 rev 3、PD-OI-054 rev 2）を co-change。Editor を Preview 上端から開き検索は残り幅（最大 39.5rem）。Align/Review を同じ行へ。2 つの常時説明を help-tip へ、Current 状態と要約を削除して `aria-pressed`。docs・capture の screenshot 再生成は S07 |
| S07 | commit 済み・dev 未反映（settlement 予算 1 行の未達は Owner が受容） | [results/S07.md](results/S07.md)。commit d728fc3f・690c610f・8051e545・ac55a72a・e8c0d7e1 と結果記録。統合 journey、S00 harness（1 行を除き予算内）、CW 変異 2 件検出、Gallery・docs の再生成。PR は S01 → S02 → #641 → S03〜S07 の順（第 9 節） |

各セッション終了時に、結果文書へのリンク、実際のSHA、検証結果、次の条件をこの表へ反映する。
完了したcommitのSHAは次のセッションで追記してよい。自分自身のcommit SHAを文書へ埋め込むためのamendを繰り返さない。
