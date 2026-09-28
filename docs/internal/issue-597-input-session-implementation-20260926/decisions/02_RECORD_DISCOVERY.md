# Approved Product Decision — Record transform discoverability

Status: human choice approved; base authority integration pending.

- 承認者: `satoshikawato`。決定日: `2026-09-26`。
- 選択: `A / REVEAL_APPLICABLE_SINGLE_RECORD_CONTROLS`、concern `diagram-generation.circular-transform-discoverability`、scenario revision `1`。
- 人の承認記録: プロジェクト所有者が2026-09-26に、以下の本文を持つ推奨 outcome の採用を明示的に承認した。承認の原文は「推奨案で承認します」。この記録は下記の具体的な選択と全文に結び付く。
- 承認対象資料の SHA-256: `57debd7625f007c99849b8de5f992508b369b7a2277c01ed1ff51e7c95354d93`。
- 既存 receipt の rationale、preservation、retirement、risk、owner、date を変更していない。手書き署名は捏造しない。
- この Markdown は承認内容の実行計画用記録であり、checker が読む第二の decision registry ではない。
- runtime 開始前に S00 の経路で既存の static Product Contract または適用可能な durable store に統合する。
- Product 選択の再承認は不要。Architecture Ratchet の exact-head exception、性能適合、privileged operator 許可、PR merge は別の条件である。

## Complete human receipt

```text
PRODUCT_DECISION
Concern: diagram-generation.circular-transform-discoverability
Scenario revision: 1
Choice: A / REVEAL_APPLICABLE_SINGLE_RECORD_CONTROLS
Rationale: source の探索状況と一件用 crop の適用条件を区別して示し、編集可能になった一件用 controls は selection の直後に見えるようにする。通常 upload の自動探索と保存済みプレビューの軽い閲覧を両立する。
Must preserve: valid native upload の自動 record discovery と Generate 前の rotation controls、既存 parser/helper 境界、exact source-bound identity、explicit single/grid/batch、fresh shared-canvas default、saved explicit choices、一件 crop と topology/start/reverse の適用条件、手動 close/expand、元の focus、keyboard/390 px、preview-only Load の Python Worker 0、active draft と saved artifact の分離、失敗時の旧 Result、Retry/Replace/Remove/Inspect/Generate の継続。
May retire: 適用可能になった一件用 section が常に collapsed で始まる挙動、valid fresh upload に manual Load が必須であるかのような prompt、実行中でない deferred discovery を loading と表す UI。全 record subset editing や複数 source support の選択は含まない。
Accepted residual risk: applicable になった時に一件用 section が展開されて pane 高さが変わる。操作元の focus と scroll anchor を維持し、無関係な更新で再展開しない。科学的意味の変更、先頭 record の自動選択、grouping の自動切替、preview-only Load による Python 初期化は受容しない。
Owner: satoshikawato
Decision date: 2026-09-26
```

## Serialized representation for review

以下は人の receipt を機械項目へ忠実に変換した inert representation。
認識済み authority file に直接追加したものではなく、S00 で target namespace と schema に合わせて取り込む。

```json
{
  "concern": "diagram-generation.circular-transform-discoverability",
  "scenarioRevision": 1,
  "choice": "A / REVEAL_APPLICABLE_SINGLE_RECORD_CONTROLS",
  "rationale": "source の探索状況と一件用 crop の適用条件を区別して示し、編集可能になった一件用 controls は selection の直後に見えるようにする。通常 upload の自動探索と保存済みプレビューの軽い閲覧を両立する。",
  "mustPreserve": "valid native upload の自動 record discovery と Generate 前の rotation controls、既存 parser/helper 境界、exact source-bound identity、explicit single/grid/batch、fresh shared-canvas default、saved explicit choices、一件 crop と topology/start/reverse の適用条件、手動 close/expand、元の focus、keyboard/390 px、preview-only Load の Python Worker 0、active draft と saved artifact の分離、失敗時の旧 Result、Retry/Replace/Remove/Inspect/Generate の継続。",
  "mayRetire": "適用可能になった一件用 section が常に collapsed で始まる挙動、valid fresh upload に manual Load が必須であるかのような prompt、実行中でない deferred discovery を loading と表す UI。全 record subset editing や複数 source support の選択は含まない。",
  "acceptedResidualRisk": "applicable になった時に一件用 section が展開されて pane 高さが変わる。操作元の focus と scroll anchor を維持し、無関係な更新で再展開しない。科学的意味の変更、先頭 record の自動選択、grouping の自動切替、preview-only Load による Python 初期化は受容しない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

## Implementation boundary

全体仕様・受け入れ条件は [MASTER_PLAN.md](../MASTER_PLAN.md)、作業規約は
[SESSION_WORKFLOW.md](../SESSION_WORKFLOW.md)。この decision は別 concern の選択を代行しない。
