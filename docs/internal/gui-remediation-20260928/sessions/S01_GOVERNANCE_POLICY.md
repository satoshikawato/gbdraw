# S01 INSTRUCTION PROMPT: 同時更新の規約と計算責務の契約

製品契約とruntimeを同じPRで更新できる限定例外、および重複計算を防ぐ規範を、文書だけで定義してください。現行checkerはchecker実装とauthority文書の同時変更を拒否するため、このセッションではchecker/runtimeを変更しません。

## 作業場所と入力

**fix/gui-feedback-remediation-20260928** と /mnt/c/Users/genom/GitHub/gbdraw/.worktrees/gui-feedback-remediation-20260928 を継続使用する。clone/worktree/branchを追加しない。
[総合計画](../00_MASTER_PLAN.md)、[SESSION_LOG](../SESSION_LOG.md)、results/S00.md、AGENTS/CLAUDE/Web CLAUDEと規範文書を読む。

## 文書の責任範囲

- AGENTS.mdの一律candidate-authority禁止の表現を、基準ルールが認める限定例外と整合させる。
- docs/internal/PRODUCT_IMPACT_RATCHET.md、WEB_CHANGE_POLICY.md。
- .github/pull_request_template.mdの該当手順。
- docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.mdの計算責務契約。
- gbdraw/web/CLAUDE.mdには具体ownerへの短い参照を置く。

静的Product Contract本体はこのセッションでは変更しない。現行checkerが同ファイルと全companionを拒否するため、同ファイルのlifecycle文はS03の最初の許可された同時更新で改訂する。古い「ファイルが存在しない」というbootstrap記述は事実に合わせ、今回の移行条件を明確にする。

## 定義する同時更新ルール

1. 対象は正確にdocs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.mdと、その結果を実現するruntime、関連通常tests/docs。
2. static Contractはguard inventoryに残り、人間によるReview REQUIREDを維持する。明示されたProduct選択を既存receiptで記録する。欠けたrationale/retirement/riskを推測して埋めない。
3. checker/detector/workflow、map/BD/architecture rules/allowlistの変更は例外に含めない。基準側のcheckerと機械authorityを使う。
4. candidate Markdownを機械評価のauthorityへ自動昇格させない。今回の例外は同時変更の受付と手動レビュー経路であり、全Product意味を自動検証できるとの約束ではない。
5. mapped hard contractの候補変更だけを安全証拠にする禁止は維持する。通常runtime-only/contract-onlyの経路も保持する。
6. この規約改正を取り込んだ後にS02を取り込む。製品契約と実装を分ける手順はその後の通常案件には要求しない。S02有効後は新base規則が旧Contract lifecycle/個別receiptの先行merge手続きだけを置き換えると明記する。製品結果・保持条件・互換性は変更せず、Contractの旧手続き文はS03で整合させる。

新しい署名システム、承認サービス、JSON registry、汎用policy interpreterを作る提案に拡大しない。

## 定義する計算責務

総合計画CW-01〜06を規範へ落とす。意味、owner、対象、合否基準を記載する。同一操作・同一意味入力の中間値の重複構築を禁止し、status/selectionだけの全feature走査を禁止する。

自動検証の初期範囲はlabel override表とselector metadata/index、指定UI操作に限定する。独立したbefore/after観測は値が同じ場合も正当とし、新しいGenerate、retry、必須validation/sanitizer/renderも維持する。cache導入を義務にせず、consumer削除とowner内の操作ローカル再利用を優先する。未計測を0として通す検証を禁止し、positive controlと変異による感度確認を要求する。

## 検証と終了

規範同士・PR template・AGENTSに矛盾がないかレビューする。現行base checkerでdocs-only差分を確認し、候補checkerで自己許可しない。実装の具体式はS02へ渡す。

results/S01.mdとSESSION_LOGを更新し、文書だけをcommitする。英語title/summary、適用パス、基準SHAを残す。
dev取り込みは明示された公開/merge承認に従う。未承認なら具体的差分と検証済みPR案まで完成させる。S02はS01が基準devへ入るまで開始しない。待機中は独立調査を進めてよいがruntimeをこの差分へ混ぜない。
