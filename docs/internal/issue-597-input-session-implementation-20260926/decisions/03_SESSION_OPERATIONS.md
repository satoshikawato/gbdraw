# Approved Product Decision — Session operation consistency

Status: human choice approved; base authority integration pending.

- 承認者: `satoshikawato`。決定日: `2026-09-26`。
- 選択: `A / EXCLUSIVE_SEMANTIC_SESSION_OPERATION`、concern `web.session-operation-consistency`、scenario revision `1`。
- 人の承認記録: プロジェクト所有者が2026-09-26に、以下の本文を持つ推奨 outcome の採用を明示的に承認した。承認の原文は「推奨案で承認します」。この記録は下記の具体的な選択と全文に結び付く。
- 承認対象資料の SHA-256: `7cef6306aefc5ff935d949efef7f76964eae2d24df8281c6ca371df5c2e87633`。
- 既存 receipt の rationale、preservation、retirement、risk、owner、date を変更していない。手書き署名は捏造しない。
- この Markdown は承認内容の実行計画用記録であり、checker が読む第二の decision registry ではない。
- runtime 開始前に S00 の経路で既存の static Product Contract または適用可能な durable store に統合する。
- Product 選択の再承認は不要。Architecture Ratchet の exact-head exception、性能適合、privileged operator 許可、PR merge は別の条件である。

## Complete human receipt

```text
PRODUCT_DECISION
Concern: web.session-operation-consistency
Scenario revision: 1
Choice: A / EXCLUSIVE_SEMANTIC_SESSION_OPERATION
Rationale: Save/Load は一つの整合した document に対する操作として完了させ、異なる時点の Result、draft、source、cache を混合しない。処理中は閲覧を維持し、semantic edits を終了後に再開する明確な workflow を優先する。
Must preserve: 主スレッドの応答と閲覧・scroll・pan/zoom・検索、visible pending/busy reasons、同時 Save の join と一度の download、title/size/repeat-download 取消、committed Result と editable draft の分離、supported Sessions と settings-only、atomic Load、failed/canceled/stale/teardown recovery、旧 request/resources/Result/History、source bytesと全 cache/evidence/provenance、JSON/gzip と CLI/Python replay、privacyとsize/sanitization constraints、現行 performance gates。
May retire: Save/Load pending 中の source/editor/History/Reset/Generate 等の semantic mutation と、mutation entry point によって偶然編集可能または silent no-op になる振る舞い。Generate/automatic reflow 中の Save/Load 開始も busy reason 付きで停止し、完了後の再試行を提供する。通常編集・閲覧・成功後の操作は廃止しない。
Accepted residual risk: 長い Save/Load の間、document 編集は一時停止する。閲覧、status、bounded completion、error/retry を維持する。無期限 lock、main-thread freeze、データの省略、checkpoint混合、追加memoryの未計測、既存gateの弱化は受容しない。
Owner: satoshikawato
Decision date: 2026-09-26
```

## Serialized representation for review

以下は人の receipt を機械項目へ忠実に変換した inert representation。
認識済み authority file に直接追加したものではなく、S00 で target namespace と schema に合わせて取り込む。

```json
{
  "concern": "web.session-operation-consistency",
  "scenarioRevision": 1,
  "choice": "A / EXCLUSIVE_SEMANTIC_SESSION_OPERATION",
  "rationale": "Save/Load は一つの整合した document に対する操作として完了させ、異なる時点の Result、draft、source、cache を混合しない。処理中は閲覧を維持し、semantic edits を終了後に再開する明確な workflow を優先する。",
  "mustPreserve": "主スレッドの応答と閲覧・scroll・pan/zoom・検索、visible pending/busy reasons、同時 Save の join と一度の download、title/size/repeat-download 取消、committed Result と editable draft の分離、supported Sessions と settings-only、atomic Load、failed/canceled/stale/teardown recovery、旧 request/resources/Result/History、source bytesと全 cache/evidence/provenance、JSON/gzip と CLI/Python replay、privacyとsize/sanitization constraints、現行 performance gates。",
  "mayRetire": "Save/Load pending 中の source/editor/History/Reset/Generate 等の semantic mutation と、mutation entry point によって偶然編集可能または silent no-op になる振る舞い。Generate/automatic reflow 中の Save/Load 開始も busy reason 付きで停止し、完了後の再試行を提供する。通常編集・閲覧・成功後の操作は廃止しない。",
  "acceptedResidualRisk": "長い Save/Load の間、document 編集は一時停止する。閲覧、status、bounded completion、error/retry を維持する。無期限 lock、main-thread freeze、データの省略、checkpoint混合、追加memoryの未計測、既存gateの弱化は受容しない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

## Implementation boundary

全体仕様・受け入れ条件は [MASTER_PLAN.md](../MASTER_PLAN.md)、作業規約は
[SESSION_WORKFLOW.md](../SESSION_WORKFLOW.md)。この decision は別 concern の選択を代行しない。
