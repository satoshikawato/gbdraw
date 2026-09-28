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

| Session | 状態 | 証拠・次の条件 |
| --- | --- | --- |
| S00 | 完了（証拠採取のみ） | [results/S00.md](results/S00.md)。開始 a351c01d、dev 57cef3ba・main 4556e04e の snapshot/wheel、性能 baseline、G03/G06/G07 再現、R01/R02 範囲。S04 の値変更は第 11.3 節の判断 1〜3 待ち。commit SHA は次セッションで追記 |
| S01 | 未着手 | S00後。規約のみの取り込み |
| S02 | 未着手 | S01がdevへ反映後。checkerのみの取り込み |
| S03 | 未着手 | S02反映後。応答・CW検証 |
| S04 | 未着手 | S00の旧Session方針とS02反映後 |
| S05 | 未着手 | canonical edge/provenanceの実装とbrowser確認 |
| S06 | 未着手 | S02反映後。Product契約とUIを同時改修 |
| S07 | 未着手 | S03〜S06の統合検証 |

各セッション終了時に、結果文書へのリンク、実際のSHA、検証結果、次の条件をこの表へ反映する。
完了したcommitのSHAは次のセッションで追記してよい。自分自身のcommit SHAを文書へ埋め込むためのamendを繰り返さない。
