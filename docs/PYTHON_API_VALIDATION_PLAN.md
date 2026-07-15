# Python API Validation and Typed Request Plan

- 作成日: 2026-07-15
- 対象: `DiagramOptions` 監査後の validation backlog と session bridge の再開条件
- 対象バージョン: `0.14.0b0` 以降
- 状態: Phase 1 と Phase 2/E0 完了。canonical request/schema 設計へ移行

## 0. 目的

既存の Python API 改修と follow-up は完了した。次の優先事項は、公開 high-level
builder が利用できない option を黙って無視しないことと、単一 Circular の depth 入力を
欠落なく検証することである。その契約を固定した後、session bridge の前提となる CLI
非依存 typed request model を設計する。

### 0.1 Phase 1 実施結果

2026-07-15 に high-level builder validation を実装した。公開 symbol、field、default、
low-level assembler signature は変更していない。

- wrong-mode field: Circular builder 2 種 × Linear-only 21 field と、Linear builder ×
  Circular-only 10 field の 52 case を追加した。
- single Circular depth: ambiguous/lossy input 9 case を追加した。
- Python API targeted suite: `507 passed`。
- full non-slow suite: `1343 passed, 16 skipped, 6 deselected`。
- `ruff check gbdraw/`: 成功。
- reference SVG comparison: 13 case 成功、生成 13 case は opt-in がないため skip。
- `tests/reference_outputs/`: 差分なし。

## 1. Phase 1: high-level builder の fail-fast validation（完了）

### 1.1 mode-specific option

- `build_circular_diagram` と `build_circular_multi_diagram` は、Linear-only field が
  non-default の場合に `ValidationError` を送出する。
- `build_linear_diagram` は、Circular-only field が non-default の場合に
  `ValidationError` を送出する。
- default 値は共有 `DiagramOptions` の互換性のため許可する。
- 例外には builder mode と該当 field 名を含め、修正方法を判断できるようにする。
- low-level `assemble_*` の signature と既定値は変更しない。

対象 field は [`DIAGRAM_OPTIONS_AUDIT.md`](DIAGRAM_OPTIONS_AUDIT.md) の分類を source of
truth とする。Linear-only は `depth_track_heights` と comparison/protein/collinearity 21 field、
Circular-only は conservation と definition control 10 field である。

### 1.2 単一 Circular depth input

`build_circular_diagram` は次を assembly 前に拒否する。

1. singular (`depth_table` / `depth_file`) と plural (`depth_tables` / `depth_files`) の併用。
2. `depth_tables` または `depth_files` の要素数が 1 以外。
3. table form と file form の併用。

これにより、plural sequence の index 0 だけを使って残りを黙って捨てる挙動を廃止し、
既存の `assemble_circular_diagram_from_records` の単一 record validation と揃える。

### 1.3 互換性

これは beta API の behavior correction として release note に記載する。正しい mode に
指定した option、全 default、単一要素の plural depth input の結果は変更しない。
`DiagramOptions` の field、型、default、`gbdraw.api.__all__` は変更しない。

### 1.4 完了条件

- wrong-mode field 31 件の non-default case が field 名を含む `ValidationError` になる。
- 全 default の `DiagramOptions` は 3 builder で従来どおり利用できる。
- 単一 Circular depth の conflict、0 件、2 件以上を assembly 前に検出する。
- 正しい mode の non-default forwarding test が維持される。
- public contract、既定 reference SVG、CLI behavior に差分がない。

## 2. Phase 2: CLI 非依存 typed request model の設計ゲート

Phase 1 完了後、session version 27～30 の GUI/CLI sidecar を対象に compatibility matrix を
作る。最初の成果物は公開 symbol ではなく、次の設計判断を含む ADR 更新とする。

1. GenBank、GFF3+FASTA、records table と selector/region/reverse を表す record input spec。
2. Circular と Linear を分離した render request。
3. output request と temporary materialization lifetime。
4. `config` / `files` から復元できず `cliInvocation.args` だけにある情報の一覧。
5. legacy session を lossless に変換できない場合の明示的な非対応境界。

public session symbol は、request が CLI option 名、argument index、parser module に依存せず、
current GUI session と canonical CLI sidecar の round-trip を証明できた後だけ追加する。

### 2.1 Phase E0 実施結果

2026-07-15 に version 27～30 の
[`Python Session Compatibility Matrix`](PYTHON_SESSION_COMPATIBILITY_MATRIX.md) を作成し、
[`session ADR`](ADR_PYTHON_SESSION_API.md) を更新した。

- version 27～30 の受理は envelope compatibility であり、lossless typed conversion の
  保証ではない。
- CLI sidecar は `cliInvocation.args` と argument index を持つ binding が完全な場合だけ、
  現行 internal path で lossless replay できる。
- GUI-only fallback は table/edit state、custom slot、output policy、cache/result の一部を
  復元しないため、typed request への lossless conversion はできない。
- version 30 には version bump 後の shape 追加があり、canonical payload を同じ version に
  追加しない。
- request model owner を `gbdraw.api.requests`、canonical payload の予定 field を
  `renderRequest`、次 schema version を 31 と決定した。
- legacy version 27～30 replay は internal のまま維持し、public session bridge の gate は
  閉じたままとする。

次の作業は [`Python Session Canonical Request Plan`](PYTHON_SESSION_CANONICAL_REQUEST_PLAN.md)
に引き継いだ。Phase R1 の request model と Phase R2 の normalization/rendering は完了し、
Phase R3 の canonical payload codec が次の実装である。

## 3. 検証

```bash
python -m pytest tests/test_api_library_usage.py tests/test_public_contract.py -q
python -m pytest tests/test_output_comparison.py::TestOutputComparison -q
ruff check gbdraw/
git diff --exit-code -- tests/reference_outputs/
```

Phase 1 の targeted gate 後、変更範囲に応じて full non-slow suite を実行する。

[Python API](./PYTHON_API.md) | [API improvement plan](./PYTHON_API_IMPROVEMENT_PLAN.md) | [Follow-up plan](./PYTHON_API_FOLLOWUP_PLAN.md)
