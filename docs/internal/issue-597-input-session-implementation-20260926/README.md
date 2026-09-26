# Issue #597 — Input, source management, and Session implementation

このディレクトリは gbdraw Web の BUG-02 の record transform controls と BUG-20 の
大規模 Session の応答性・整合性を実装するための引継ぎ資料である。
実装担当者は会話履歴を読む必要がない。
BUG-01 は2026-09-26のユーザー指示で対象外とし、専用の計画・承認記録・S02 prompt を削除した。

1. [総合計画書](./MASTER_PLAN.md): 問題、承認済み outcome、アーキテクチャ、全受け入れ条件、セッション順序。
2. [セッション共通作業規約](./SESSION_WORKFLOW.md): 専用ブランチの取得、独立 checkout、検証、コミット・プッシュ、競合回復。
3. [セッション別 INSTRUCTION PROMPTS](./sessions/README.md): S00 → S01 → S03〜S08 を順番に一つずつ実行する。
4. [変形欄の発見性](./decisions/02_RECORD_DISCOVERY.md)、[Session operation](./decisions/03_SESSION_OPERATIONS.md): 完全な人の receipt と inert machine representation。
5. [基準動作の証拠](./evidence/README.md)、[実施結果の記録先](./results/README.md)。

実装・結果保存・push の対象は `fix/issue-597-input-session-20260926`。
各セッションは remote の最新同名ブランチを取得して使う。
既存の shared checkout や別タスクの branch/worktree を切り替えない。
Product outcome は承認済み。runtime 開始には必要な authority の dev への統合と性能方式の検証が必要である。
この計画を保存する初回コミットには実行時コードや有効 authority の変更は含まれない。
