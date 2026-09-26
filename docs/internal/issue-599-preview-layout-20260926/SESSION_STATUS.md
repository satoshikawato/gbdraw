# Issue #599 — セッション進捗

対象ブランチ: `fix/issue-599-preview-layout-20260926`。
計画基準: `origin/dev@d457b7189b137185a8dec800819a312c30b969fa`。

2026-09-26の初期状態。製品結果の3案 A は承認済み。runtime、恒久 authority の追加、修正後の受入検証は未着手である。
担当者は開始前に前セッションの commit/push と単一 writer を確認し、終了時にこの表と担当結果文書を更新する。

| セッション | 状態 | 結果文書 | 次の開始条件 |
| --- | --- | --- | --- |
| 計画作成 | 文書作成済み。公開 SHA は Git 履歴・終了報告で確認 | MASTER_PLAN / 個別 prompts / 承認 Pack | 計画 commit が同名 remote branch に存在 |
| S00 / authority | 未着手 | 終了時に `results/S00.md` を作成 | 全必要 authority の dev merge を確認するまで S01 runtime は開始不可 |
| S01 / 配置継承 | 未着手 | 終了時に `results/S01.md` を作成 | S01 commit/push、C01–C08/R01 の担当範囲の証拠 |
| S02 / 操作説明 | 未着手 | 終了時に `results/S02.md` を作成 | S02 commit/push、A01–A02 の証拠 |
| S03 / chrome | 未着手 | 終了時に `results/S03.md` を作成 | S03 commit/push、P01–P03 の証拠 |
| S04 / 統合 | 未着手 | 終了時に `results/S04.md` を作成 | 全受入と適用 gate、文書・diff review。未達の merge後 staging は別記 |

結果文書は事実だけを記載する。入力 SHA、検査対象 SHA/環境、authority base、コマンドと結果、artifact、制限、未完了条件、次の担当者の開始条件を残す。
結果文書を含む自身の commit SHA を事前に埋め込む必要はない。公開 SHA は `git log` と終了報告で追跡する。

[総合計画](MASTER_PLAN.md) / [承認済み製品判断](00_APPROVED_PRODUCT_DECISIONS.md)。

## 計画作成時の検証

文書11ファイル、関心別Pack3件、セッション別prompt5件を確認した。承認本文は元の選択Aと全文一致し、JSONの全9フィールドと一致する。local参照51件、bash code block15件のsyntaxを確認した。基準観測JSONは未変更、観測runnerは説明文のみ更新し、実行処理のASTは同一である。runtimeは変更していないため修正後の受入試験は未実施。
