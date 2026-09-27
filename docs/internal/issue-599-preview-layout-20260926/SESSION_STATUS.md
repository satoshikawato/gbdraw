# Issue #599 — セッション進捗

対象ブランチ: `fix/issue-599-preview-layout-20260926`。
計画基準: `origin/dev@d457b7189b137185a8dec800819a312c30b969fa`。

2026-09-28更新。製品結果の3案 A を OIPC revision 25 の PD-OI-052–054 に記録・検証し、[PR #634](https://github.com/satoshikawato/gbdraw/pull/634) の必須 checks 成功後に dev `007388567222638b707fbb16fe82dbeba61551c9` へ merge 済み。**S00 完了、次は S01**。runtime と修正後の受入検証は未着手である。
担当者は開始前に前セッションの commit/push と単一 writer を確認し、終了時にこの表と担当結果文書を更新する。

| セッション | 状態 | 結果文書 | 次の開始条件 |
| --- | --- | --- | --- |
| 計画作成 | 文書作成済み。公開 SHA は Git 履歴・終了報告で確認 | MASTER_PLAN / 個別 prompts / 承認 Pack | 計画 commit が同名 remote branch に存在 |
| S00 / authority | 完了。OIPC単独PR #634はrequired checks成功・dev merge済み。同期後のdocs-only Gate PASS | [results/S00.md](results/S00.md) | retained checkout で最新devの完全な3 outcomeと関連authorityを確認してS01へ |
| S01 / 配置継承 | 未着手 | 終了時に `results/S01.md` を作成 | S01 commit/push、C01–C08/R01 の担当範囲の証拠 |
| S02 / 操作説明 | 未着手 | 終了時に `results/S02.md` を作成 | S02 commit/push、A01–A02 の証拠 |
| S03 / chrome | 未着手 | 終了時に `results/S03.md` を作成 | S03 commit/push、P01–P03 の証拠 |
| S04 / 統合 | 未着手 | 終了時に `results/S04.md` を作成 | 全受入と適用 gate、文書・diff review。未達の merge後 staging は別記 |

結果文書は事実だけを記載する。入力 SHA、検査対象 SHA/環境、authority base、コマンドと結果、artifact、制限、未完了条件、次の担当者の開始条件を残す。
結果文書を含む自身の commit SHA を事前に埋め込む必要はない。公開 SHA は `git log` と終了報告で追跡する。

[総合計画](MASTER_PLAN.md) / [承認済み製品判断](00_APPROVED_PRODUCT_DECISIONS.md)。

## 計画作成時の検証

文書11ファイル、関心別Pack3件、セッション別prompt5件を確認した。承認本文は元の選択Aと全文一致し、JSONの全9フィールドと一致する。local参照51件、bash code block15件のsyntaxを確認した。基準観測JSONは未変更、観測runnerは説明文のみ更新し、実行処理のASTは同一である。runtimeは変更していないため修正後の受入試験は未実施。

## S00 handoff and checkout reuse

S00 の authority base は `origin/dev@ecb96a065d808addff0fb088f9eb6b6ba3c92a6d`、入力 HEAD は計画 commit `bcc4e0aa5ffcf4bdf9a12952caa8bf5de2cfee08`。最新 dev は通常 merge で取得し、その上で3ファイルの authority/docs 差分を作成した。既存51 record と acceptance catalog、runtime、checker、workflow、tests、generated assets は dev と比較して不変。詳細・検査コマンド・各 preservation contribution・初回 Gate 失敗とその解消 は [S00](results/S00.md) に記録する。

依頼者の追記により、**次回以降は `/mnt/c/users/genom/github/gbdraw-issue599-s00` を再利用する**。再開時に旧 `/tmp` checkout が失われていたため公開済み branch をこの永続パスへ復元した。毎回の clone は行わない。総合計画の共通開始手順と S01–S04 prompt を更新した。各担当者は前セッション終了、単一writer、clean tree、同名branch/upstream、fetch / ff-only pullを確認する。共有作業ツリーは変更しない。

初回の混在差分の Gate FAIL は免除せず、同じ S00 内で OIPC 単独候補を準備した。追加の明示承認により、候補 `32e2af3c347ea4c526d69ff0217ae07db76d7720` を authority branch へ公開し、PR #634 の required checks 成功後に dev merge した。dev `007388567222638b707fbb16fe82dbeba61551c9` は実装ブランチへ通常merge `bcc19d9613990cfc9409af3a1e0620888ee1f194` で取り込み済み。後続 S01 は最新devの3つの完全なoutcomeと関連authorityを確認する。S00 を再実行する必要はない。
