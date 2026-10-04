# S02 INSTRUCTION PROMPT: Product契約とruntimeの同時変更をCIで許可する

S01の規約を実装してください。既存checkerの限定的な分岐とfixtureを改修し、Product契約・runtime・通常tests/docsを一緒にレビューできるようにします。製品runtimeとProduct契約内容は変更しないセッションです。

## 前提と場所

**fix/gui-feedback-remediation-20260928**、/mnt/c/Users/genom/GitHub/gbdraw/.worktrees/gui-feedback-remediation-20260928を再利用。
[総合計画](../00_MASTER_PLAN.md)、[SESSION_LOG](../SESSION_LOG.md)、results/S01.mdとリポジトリ指示を読む。
originをfetchし、S01がorigin/devに入った実際のSHAを確認して同じbranchへ取り込む。候補規約がworktreeにあるだけでは前提成立としない。S01の移行条項に従い、S02有効後は旧Contractの先行merge手続きだけを置換する。製品結果・保持条件を上書きするものではない。

## 実装所有者

tools/check-web-change-budget.mjsと既存のchecker fixture tests、特にtests/web/architecture-contracts.test.mjs。
現行コードのproductContractAuthorityPath、candidateAuthorityCompanionPaths、changedGuards/changedAuthoritiesを確認する。

1. static Product Contractの全companion拒否とruntime+guard拒否の双方に、同じ狭い許可条件を適用する。二つの場所に異なる例外規則を実装しない。
2. 変更したgovernance/guardが正確なstatic Product Contractだけであり、companionがruntimeと関連する非guard tests/docsである場合に受付可能とする。
3. contract変更はReview REQUIREDにする。既存human review/receiptを利用し、Markdown意味の新parserは追加しない。
4. その他のchecker/detector/workflow/map/BD/rule/allowlist変更の隔離、trusted-base実行、mapped hard coverageは維持する。
5. reportが「すべてのcandidate authorityはinertのみ」「必ずcontract-only」の旧説明を出し続けないよう、静的な手動レビュー例外と機械authorityの違いを表示する。

## 必須fixture

| 入力差分 | 合否 |
| --- | --- |
| static Contract + runtime + 通常tests/docs | Gateは他の条件を満たせば通過、Review REQUIRED |
| static Contractのみ | 既存経路を維持 |
| runtimeのみ | 既存経路を維持 |
| static Contract + runtime + checker/detector/workflow | 拒否 |
| static Contract + runtime + map/BD/rules/allowlist | 拒否 |
| 類似した別ファイル名・pathの削除/移動で例外を拡張 | 拒否または従来分類。exact pathを迂回しない |
| mapped checkpointの候補テストだけを差し替え | 既存のhard coverage不足を検出 |
| static Contract変更のreview理由を落とす変異 | fixture失敗 |

「完全なhuman choiceがない場合の製品承認」は人間のReview責務。現在parserが読まない内容を自動検証済みと報告しない。機械Gateと人間Reviewの役割を結果に明示する。

## 検証とhandoff

既存architecture/product-impact fixtureを対象にnode --testで実行し、synthetic repository差分でbase checker経路を検証する。新しいchecker自身を使ったunit検証と、実際のPRを判定するtrusted-baseの適用時点を区別する。

results/S02.mdとSESSION_LOGへ許可/拒否matrix、コマンド、base/headを記録。担当差分をcommitし、英語title/summaryを残す。規範文書・static Contract・runtimeをchecker-only差分へ混ぜない。

承認された経路でdevへ取り込んだ後、そのSHAを同じ作業branchへ取り込む。S03以降はこの改正済みbaseで製品契約・runtime・通常testを同時更新する。追加clone/branchは不要。
