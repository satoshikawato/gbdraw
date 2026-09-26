# Approved Product Decision — Circular source collection

Status: human choice approved; base authority integration pending.

- 承認者: `satoshikawato`。決定日: `2026-09-26`。
- 選択: `A / INDEPENDENT_ORDERED_SOURCES`、concern `diagram-generation.circular-source-collection`、scenario revision `1`。
- 人の承認記録: プロジェクト所有者が2026-09-26に、以下の本文を持つ推奨 outcome の採用を明示的に承認した。承認の原文は「推奨案で承認します」。この記録は下記の具体的な選択と全文に結び付く。
- 承認対象資料の SHA-256: `44ec756a4f0be6475e21508fe853db6fa3760bae4bf938d4a592338c20903851`。
- 既存 receipt の rationale、preservation、retirement、risk、owner、date を変更していない。手書き署名は捏造しない。
- この Markdown は承認内容の実行計画用記録であり、checker が読む第二の decision registry ではない。
- runtime 開始前に S00 の経路で既存の static Product Contract または適用可能な durable store に統合する。
- Product 選択の再承認は不要。Architecture Ratchet の exact-head exception、性能適合、privileged operator 許可、PR merge は別の条件である。

## Complete human receipt

```text
PRODUCT_DECISION
Concern: diagram-generation.circular-source-collection
Scenario revision: 1
Choice: A / INDEPENDENT_ORDERED_SOURCES
Rationale: Circular の複数 GenBank は元のファイルを管理単位とし、一ファイルの更新で別ファイルの設定を失わない操作を提供する。source と record の区別を入力、生成、保存、復元の全工程で保つ。
Must preserve: 元の source 内 record order と全 records、同名・重複 accession・同一 bytes の別 instance、既存単一 GenBank と単一 GFF3/FASTA、single/grid/batch と保存済み明示選択、single-record crop、source-bound transforms、depth sparse alignment、比較 intent、対象外 source の状態、active draft と committed Result の分離、Undo/Redo、supported released Sessions と CLI/Python replay。
May retire: 新 writer の scalar c_gb inventory と複数原ファイルを一つの合成入力に見せる current presentation。released scalar/composite readers は廃止しない。branch-only Session 44/bindings 2 artifacts は current writer に再生成し、その組だけの新 compatibility reader は作らない。
Accepted residual risk: bindings namespace に source collection と released input の bounded normalization が必要になる保守負担を受容する。大きい record lists と交換対象の dependency reconciliation は focused gates で管理する。科学的出力の誤対応、データ喪失、source bytes 改変、別 source の intent 喪失、Architecture Ratchet failure は受容しない。
Owner: satoshikawato
Decision date: 2026-09-26
```

## Serialized representation for review

以下は人の receipt を機械項目へ忠実に変換した inert representation。
認識済み authority file に直接追加したものではなく、S00 で target namespace と schema に合わせて取り込む。

```json
{
  "concern": "diagram-generation.circular-source-collection",
  "scenarioRevision": 1,
  "choice": "A / INDEPENDENT_ORDERED_SOURCES",
  "rationale": "Circular の複数 GenBank は元のファイルを管理単位とし、一ファイルの更新で別ファイルの設定を失わない操作を提供する。source と record の区別を入力、生成、保存、復元の全工程で保つ。",
  "mustPreserve": "元の source 内 record order と全 records、同名・重複 accession・同一 bytes の別 instance、既存単一 GenBank と単一 GFF3/FASTA、single/grid/batch と保存済み明示選択、single-record crop、source-bound transforms、depth sparse alignment、比較 intent、対象外 source の状態、active draft と committed Result の分離、Undo/Redo、supported released Sessions と CLI/Python replay。",
  "mayRetire": "新 writer の scalar c_gb inventory と複数原ファイルを一つの合成入力に見せる current presentation。released scalar/composite readers は廃止しない。branch-only Session 44/bindings 2 artifacts は current writer に再生成し、その組だけの新 compatibility reader は作らない。",
  "acceptedResidualRisk": "bindings namespace に source collection と released input の bounded normalization が必要になる保守負担を受容する。大きい record lists と交換対象の dependency reconciliation は focused gates で管理する。科学的出力の誤対応、データ喪失、source bytes 改変、別 source の intent 喪失、Architecture Ratchet failure は受容しない。",
  "owner": "satoshikawato",
  "decisionDate": "2026-09-26"
}
```

## Implementation boundary

全体仕様・受け入れ条件は [MASTER_PLAN.md](../MASTER_PLAN.md)、作業規約は
[SESSION_WORKFLOW.md](../SESSION_WORKFLOW.md)。この decision は別 concern の選択を代行しない。
