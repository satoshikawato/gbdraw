# Issue #599 — 承認済み製品判断の一覧

Product Decision Owner `satoshikawato` は2026-09-26に下記3件それぞれの推奨 A を明示承認した。
判断本文、完全な選択 outcome、preservation/retirement/risk、非実行 JSON のレビュー表現は、関心ごとに分けた Pack に保存する。
この一覧は新しい一括 authority ではない。各 concern の責任と受入条件を個別に追跡する。

| Concern | 承認された選択 | Decision Pack |
| --- | --- | --- |
| `web.composition-decoration-continuity` | `A / CARRY-MATCHED-DECORATION-DELTAS` | [装飾配置の継承](DECISION_01_COMPOSITION_CONTINUITY.md) |
| `web.layout-edit-affordance` | `A / EXPLICIT-MODE-WITH-DISCOVERABLE-TARGETS` | [Layout edit の発見](DECISION_02_LAYOUT_AFFORDANCE.md) |
| `web.preview-search-placement` | `A / DOCKED-SEARCH-AND-CONTROLS` | [検索・toolbar の専用行](DECISION_03_PREVIEW_CHROME.md) |

再署名待ちの提案ではない。残る作業は [S00](SESSION_00_AUTHORITY_INSTRUCTION_PROMPT.md) による恒久契約への正確な記録と、dev merge の確認である。
このディレクトリの JSON code block は承認回答の表示用で、runtime や checker は参照しない。
機械フィールドは PRODUCT_DECISION 本文の1対1の転記であり、新しい実行契約ではない。

[総合計画](MASTER_PLAN.md) と [進捗](SESSION_STATUS.md) を参照する。
