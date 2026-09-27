# Issue #602 実装セッション用 INSTRUCTION PROMPTS

各ファイルは新しい実装セッションへ全文を渡して使う。実行順はS00→S01→S02→S03→S04→S05→S06→S07。各本文は対象、署名済み結果、前提、責務、検証、handoffを説明する。

[総合計画](../00_MASTER_PLAN.md)は実装契約の所有文書、[署名記録](../06_SIGNED_PRODUCT_DECISIONS.md)はProduct選択の正確な全文と機械表現。01-B、02-A、03-A、04-A、05-Aが採択済みで、選択の再承認は不要。依存runtimeの開始にはS00の正式authorityがorigin/devへmergeされている必要がある。

| Session | 指示書 | 成果 |
| --- | --- | --- |
| S00 | [Authority](SESSION_00_AUTHORITY.md) | 正式Contractへのexact serialization。runtimeを含まない |
| S01 | [Auto disclosure](SESSION_01_AUTO_DISCLOSURE.md) | Auto初期値を維持し、非表示理由とLabels導線を表示 |
| S02 | [Definition default](SESSION_02_DEFINITION_DEFAULT.md) | Web fresh/reset Lock ON、保存値と他surfaceを保持 |
| S03 | [Canonical Status projection](SESSION_03_CANONICAL_STATUS_PROJECTION.md) | request共用projectionとResult/draft比較基準 |
| S04 | [Application feedback](SESSION_04_APPLICATION_FEEDBACK.md) | Pending/live/review適用説明とGenerate/Save/Export意味 |
| S05 | [Compact Editor](SESSION_05_COMPACT_EDITOR.md) | 同一Editorと図の上下配置 |
| S06 | [Compact alignment review](SESSION_06_COMPACT_ALIGNMENT_REVIEW.md) | local draft/retryを保ち、reviewとcanvasを同時操作 |
| S07 | [Acceptance and docs](SESSION_07_ACCEPTANCE_AND_DOCS.md) | 統合検証、既存public文書更新、最終diff review |

実装ブランチは、計画commitを含む既存の `fix/issue-602-linear-live-edit-20260926`。S01で正式authorityがmerge済みのorigin/devをこのブランチへmergeし、S02–S07も同じブランチを継続する。計画commitを保持し、ブランチを作り直さない。upstreamとpush先は同名の `origin/fix/issue-602-linear-live-edit-20260926` のみとする。S00だけは別のauthority-only worktreeを使用し、計画文書を候補へ持ち込まない。S00の成果記録は実装ブランチの計画checkoutに保存する。実装指示の作成は後続のpush/PR/merge権限を含まない。

一つのセッションが終わっただけで全体完了とは扱わない。検証証拠にはcommand、対象HEAD、入力、環境、受入条件を対応付け、次の変更で無効になった範囲だけを更新する。
