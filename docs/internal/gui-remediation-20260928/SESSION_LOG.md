# 実行状況

計画: [総合実装計画](00_MASTER_PLAN.md)。
作業ブランチ: fix/gui-feedback-remediation-20260928。
再利用 worktree: /mnt/c/Users/genom/GitHub/gbdraw/.worktrees/gui-feedback-remediation-20260928。
初期 dev: 57cef3ba47f4b7790a9145f2ce6988422a00710e。

2026-09-28: 計画書とセッション別プロンプトを作成。runtime・規約・checkerは未変更。
基準のソース調査と限定resolver確認は計画に記載した。S00の再現可能なbaseline採取は未実施。
Product契約の同時更新例外、計算重複防止契約も未実装であり、以下の順に実施する。

計画書作成時の検証: 10 Markdownの相対リンク・書式・branch指定・要求IDを確認済み。git diff --cached --checkは成功。node tools/check-web-change-budget.mjs --base origin/devはGate PASS / Review CLEAR。独立レビューで計算契約の範囲、旧Session、規約移行、BGC手順を確認した。runtime差分がないため製品テストはこの計画書commitでは実行していない。

| Session | 状態 | 証拠・次の条件 |
| --- | --- | --- |
| S00 | 未着手 | baseline・旧Session・方向の再現、判断範囲 |
| S01 | 未着手 | S00後。規約のみの取り込み |
| S02 | 未着手 | S01がdevへ反映後。checkerのみの取り込み |
| S03 | 未着手 | S02反映後。応答・CW検証 |
| S04 | 未着手 | S00の旧Session方針とS02反映後 |
| S05 | 未着手 | canonical edge/provenanceの実装とbrowser確認 |
| S06 | 未着手 | S02反映後。Product契約とUIを同時改修 |
| S07 | 未着手 | S03〜S06の統合検証 |

各セッション終了時に、結果文書へのリンク、実際のSHA、検証結果、次の条件をこの表へ反映する。
完了したcommitのSHAは次のセッションで追記してよい。自分自身のcommit SHAを文書へ埋め込むためのamendを繰り返さない。
