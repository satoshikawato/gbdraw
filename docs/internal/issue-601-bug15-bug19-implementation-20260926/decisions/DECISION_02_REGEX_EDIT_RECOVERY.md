# 承認済みProduct Decision — 拒否された既存Color ruleの入力回復

Status: Product Decision Ownerが2026-09-26に選択Aを承認済み。
Product Decision Owner: satoshikawato。
対象: Issue #601 の BUG-15 / BUG-19。B案は採択していない。

## 選択された製品動作

拒否された既存Color ruleのpattern textを未確定field draftとして保持し、原因・未適用状態・Retry/Revertを提供する。

本文の各条件が実装の境界である。別の関心への承認を推測したり、残余リスクを追加したりしない。
このファイルは承認内容のレビュー・引継ぎ記録であり、新しいCI decision storeや正式Product authorityではない。
正式な未登録concernのauthorityは、既存OPTION_INTEGRITY_PRODUCT_CONTRACT.mdへの記録がdevに統合された後に有効となる。
[総合計画書](../MASTER_PLAN.md)のS01→S02の順序を守る。

## 承認された本文

以下の本文は選択Aの記入済みreceiptをそのまま保存したもの。
receipt SHA-256（UTF-8、末尾改行を除く）: 9cc66bc30f97cc995b9b0002b094ba74cfa00a8e77b163c33074689f0e73b496。

~~~text
PRODUCT_DECISION
Concern: web.rules.rejected-pattern-edit-recovery
Scenario revision: 1
Choice: A / KEEP_REJECTED_PATTERN_DRAFT
Rationale: 有効なルールと図を保護しながら、入力ミスや一時的な検証失敗から同じ文字列を修正・再試行できるようにする。
Must preserve: Color/Label の Python regex semantics、既存一 Worker・一 preparation・atomic live commit、優先順位、prepared reuse、valid target と Generate の同値、canonical rule/Result/History、stale/cancel 隔離を維持する。対象は既存 Color rule の pattern field。拒否された text を field に保持して原因と Not applied、Save/Generate は last accepted rule、Export は現在 Result を使うことを示す。keyboard の編集/Retry/Revert、正しい syntax/runtime 分類、同 document の drawer close/reopen・一時 mode 切替での draft 保持、row/revision の現在性を保つ。対象 rule の Undo/Redo 置換、row 削除、document/session 成功置換、reset で draft を解放する。成功 edit だけ History へ記録し、未確定 text は Session/diagnostics/console に自動保存・公開しない。TSV、新規 rule、preset、Search の意味は変えない。
May retire: 対象 field の failure 後に未確定 pattern text を無条件で accepted 値へ戻す表示。不正 rule の拒否、last accepted rule の保護、正常な live edit は退役しない。
Accepted residual risk: 表示 text と accepted rule が一時的に異なり、Save/Generate は accepted rule、Export は現在 Result を使う。Not applied と対象説明、Revert を提供する。Session/document や対象 rule の History 置換後に未確定 draft は保持しない。
Owner: satoshikawato
Decision date: 2026-09-26
~~~

## 機械表現

本文の九項目だけを、同じ意味のJSONへ機械的に変換した。
scenarioRevisionだけを整数にし、rationale・mustPreserve・mayRetire・acceptedResidualRiskは原文のまま保存した。
これはレビュー用の表示であり、独立したJSON registryではない。

~~~json
{
  "concern": "web.rules.rejected-pattern-edit-recovery",
  "scenarioRevision": 1,
  "choice": "A / KEEP_REJECTED_PATTERN_DRAFT",
  "rationale": "有効なルールと図を保護しながら、入力ミスや一時的な検証失敗から同じ文字列を修正・再試行できるようにする。",
  "mustPreserve": "Color/Label の Python regex semantics、既存一 Worker・一 preparation・atomic live commit、優先順位、prepared reuse、valid target と Generate の同値、canonical rule/Result/History、stale/cancel 隔離を維持する。対象は既存 Color rule の pattern field。拒否された text を field に保持して原因と Not applied、Save/Generate は last accepted rule、Export は現在 Result を使うことを示す。keyboard の編集/Retry/Revert、正しい syntax/runtime 分類、同 document の drawer close/reopen・一時 mode 切替での draft 保持、row/revision の現在性を保つ。対象 rule の Undo/Redo 置換、row 削除、document/session 成功置換、reset で draft を解放する。成功 edit だけ History へ記録し、未確定 text は Session/diagnostics/console に自動保存・公開しない。TSV、新規 rule、preset、Search の意味は変えない。",
  "mayRetire": "対象 field の failure 後に未確定 pattern text を無条件で accepted 値へ戻す表示。不正 rule の拒否、last accepted rule の保護、正常な live edit は退役しない。",
  "acceptedResidualRisk": "表示 text と accepted rule が一時的に異なり、Save/Generate は accepted rule、Export は現在 Result を使う。Not applied と対象説明、Revert を提供する。Session/document や対象 rule の History 置換後に未確定 draft は保持しない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
~~~
