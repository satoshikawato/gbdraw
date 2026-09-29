# Session INSTRUCTION PROMPTS

各ファイルを一つずつ新しい実装セッションへ渡す。前提結果を読んでから開始し、S00 → S01 → S03〜S08を順番に実行する。S02はBUG-01除外に伴い削除した。
全セッションは `fix/issue-597-input-session-20260926` の最新remoteを独立checkoutへ取得する。
担当結果を一commitにし、終了時に同名branchへpushしてremoteSHAを確認する。
共通作業規約と総合計画は各promptの必読資料。会話履歴や未保存の口頭合意は必要ない。

| Session | Prompt |
| --- | --- |
| S00 | [Authority intake and approved decision serialization](./S00_INSTRUCTION_PROMPT.md) |
| S01 | [Baseline and import transport evidence](./S01_INSTRUCTION_PROMPT.md) |
| S03 | [Discovery state and transform disclosure](./S03_INSTRUCTION_PROMPT.md) |
| S04 | [Session operation consistency](./S04_INSTRUCTION_PROMPT.md) |
| S05 | [Measured Session import Worker](./S05_INSTRUCTION_PROMPT.md) |
| S06 | [Time-bounded projection and restoration](./S06_INSTRUCTION_PROMPT.md) |
| S07 | [Reproducible discovery and Session documentation](./S07_INSTRUCTION_PROMPT.md) |
| S08 | [Integrated acceptance and final head handoff](./S08_INSTRUCTION_PROMPT.md) |

S00のauthority integrationは別のmaintainer delivery。base未反映でもS01の独立measurementを完了できる。
必要authority未反映のruntimeは停止し、準備/evidenceの成果をcommit/pushして残る境界を報告する。
