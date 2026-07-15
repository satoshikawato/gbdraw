# Python API Improvement Plan

- 作成日: 2026-07-14
- 更新日: 2026-07-15
- 対象バージョン: `0.14.0b0`
- 監査スナップショット: branch `cli_tsv_inputs`, HEAD `8f952d4` の作業ツリー
- 実装スナップショット: branch `cli_tsv_inputs`, HEAD `67cb66d` の作業ツリー
- 対象: `gbdraw/api/`, Python API から利用する入出力・設定・interactive SVG・session 経路、関連テストと文書
- 目的: 静的描画エンジンとして既に高い機能到達度を持つ `gbdraw.api` を、CLI/Web の主要ワークフローも安全に利用できる公開 API にする

## 0. 実装結果

本計画に基づく API 改修は完了した。Phase 4 の session bridge は、CLI 引数を
library contract に漏らさないため設計ゲートで停止し、非公開を維持する判断を
[ADR](ADR_PYTHON_SESSION_API.md) に記録した。これは完了の定義 7 で認めた結果である。

| 領域 | 状態 | 実装結果 |
|---|---|---|
| P1 anchor mode | 完了 | `all`、`one_to_one`、`rbh` と alias の正規化結果を collinearity owner まで転送 |
| P1 typed config | 完了 | `LabelsFilteringConfig.raw` 内の 3 種のラベル DataFrame を override 後も保持 |
| P1 export | 完了 | `save_figure_to` を strict にし、実在 path だけを返し、変換失敗を型付き例外に統一 |
| P2 interactive SVG | 完了 | `InteractiveSvgContext`、`enrich_svg`、`build_interactive_svg_context` を公開し、bytes API へ context を追加 |
| P2 table input | 完了 | records/conservation/circular-track の型・reader と 3 種の label reader、DataFrame/file options を公開 |
| P2 multi-record Circular | 完了 | `CircularMultiRecordOptions` と `build_circular_multi_diagram` を追加 |
| P3 session | 設計ゲート完了 | version 31 canonical request の両 mode round-trip が成立するまで session orchestration を内部に維持 |
| P3 documentation | 完了 | capability matrix、主要 recipe、制約を追加し、全 Python code block を pytest から実行 |

実装後の公開面は `gbdraw.api.__all__` が 87 symbol、`DiagramOptions` が 70 field、
`CircularMultiRecordOptions` が 5 field である。契約 snapshot には意図した追加だけを反映した。

検証結果:

- 実装時 follow-up targeted suite: `446 passed`
- 実装時 follow-up full non-slow suite: `1282 passed, 16 skipped, 6 deselected`
- export の最終修正後の関連 suite: `41 passed`
- repository 全体の Ruff: 既存 85 件を解消し、`ruff check gbdraw/` が成功
- reference SVG: 13 件の SHA-256 は不変。生成 test は opt-in、比較失敗時の actual SVG は
  temporary directory へ保存
- release readiness: release note と `DiagramOptions` 70-field audit を追加

2026-07-15 再検証:

- 計画指定 targeted suite（Python API 文書内の全 code block を含む）: `507 passed`
- full non-slow suite: `1383 passed, 16 skipped, 6 deselected`
- `ruff check gbdraw/`: 成功
- `tests/reference_outputs/`: tracked diff なし。comparison test も成功
- public contract: 87 symbol、`DiagramOptions` 70 field、
  `CircularMultiRecordOptions` 5 field を再確認

今回得られた保守手順は
[`maintain-python-api`](skills/maintain-python-api/SKILL.md) SKILL に整理した。

## 1. 結論と目標状態

監査時の Python API は、Circular/Linear の静的描画についてはほぼ全機能へ到達できた。CLI も最終的には次の公開アセンブラを呼んでいた。

- `assemble_circular_diagram_from_record`
- `assemble_circular_diagram_from_records`
- `assemble_linear_diagram_from_records`

不足していたのは描画能力そのものより、次の公開契約とワークフローである。

1. 公開されているが指定値が反映されない、または結果が不正確な API の修正。
2. 型付き config、ラベル用 DataFrame、interactive metadata を安全に合成する経路。
3. records/conservation/track table と session を、内部 module へ依存せず利用する入口。
4. 51～63 個の引数を持つ低位アセンブラへ直接降りなくても高度機能を使える構成。
5. 実装済みの高度機能を発見・再現できるドキュメント。

今回の実装で 1、2、4、5 と、3 のうち table 入力を解消した。session は不安定な
CLI contract を公開しない方が安全と判断し、ADR の再開条件を満たすまで内部に維持する。

公開 namespace だけを利用する現在の到達点は次のとおりである。

```text
入力を読む／メモリ上の SeqRecord を受け取る
    → 表・設定・比較データを検証する
    → Circular/Linear 図を組み立てる
    → static または rich interactive SVG を生成する
    → 実際に生成できた形式と metadata を結果として受け取る
```

session の保存・再実行は、CLI 非依存 typed request model の設計後にこの流れへ追加する。

## 2. スコープ

### 2.1 対象

- `gbdraw.api` の公開 symbol、型、既定値、例外契約。
- Circular/Linear の static diagram 機能。
- GenBank、GFF3+FASTA、record selector、reverse complement、region。
- feature color/visibility/shape、label filtering/override。
- GC content、GC skew、depth、conservation、comparison、protein comparison。
- circular/linear track slot。
- SVG、interactive SVG、PNG、PDF、EPS、PS。
- records table、conservation table、circular track table。
- GUI/CLI session の読込、検証、再実行に必要な Python 側機能の設計ゲート。
- API 契約、回帰、ドキュメント実行テスト。

### 2.2 非目標

- Web UI の DOM 操作を Python へ移植しない。
- drag & drop、undo/redo、popup、ブラウザ内の即時 preview を Python API の責務にしない。
- CLI parser を Python API として公開しない。
- 既存の全設定を新しい抽象化へ一括移行しない。
- 公開 API のためだけに renderer や configurator を複製しない。
- 互換性のない全面再設計を 1 PR で行わない。

## 3. 監査開始時の基準値

### 3.1 監査開始時の公開面

- `gbdraw.api.__all__`: 70 symbol。
- `DiagramOptions`: 64 field。
- `assemble_circular_diagram_from_record`: private 相当引数を除き 51 引数。
- `assemble_circular_diagram_from_records`: 58 引数。
- `assemble_linear_diagram_from_records`: 63 引数。
- `build_circular_diagram` と `build_linear_diagram` は `DiagramOptions` を受け取る簡易入口。

### 3.2 監査時に確認した機能

| 領域 | 現状 | 改修方針 |
|---|---|---|
| Circular static rendering | ほぼフル | 非回帰を維持 |
| Linear static rendering | ほぼフル | 非回帰を維持 |
| file/in-memory input | 利用可能 | table 入力を公開面へ追加 |
| labels/style/config | 到達可能だが low-level | DataFrame 用の明示的な入口を追加 |
| binary export | 利用可能 | 成否と戻り値の契約を修正 |
| rich interactive SVG | 内部 API が必要 | context と builder を公開 |
| session | 内部 module のみ | CLI 非依存の公開境界を設計 |
| Web post-editing | JavaScript 専用 | 非目標として維持 |

### 3.3 監査時のテスト結果

代表的な API、track、depth、conservation、comparison、interactive、selector、CLI table のテストを実行し、292 件が成功した。

```text
292 passed
```

この値を改修開始時の targeted baseline とする。全体テストの基準値は各 PR の開始時に別途記録する。

## 4. 実施原則

### 4.1 正確さを機能追加より先にする

- 公開引数が無視される状態を残したまま symbol を増やさない。
- 戻り値は「要求したもの」ではなく「実際に生成・検証できたもの」を表す。
- fallback、skip、部分成功は warning と構造化された結果で区別する。
- silent fallback が必要な CLI と、失敗検出が必要な library API の契約を分ける。

### 4.2 薄い公開層を維持する

- `gbdraw.api` は既存 owner の関数・型を再利用する。
- CLI private helper をそのまま import して公開しない。
- parser、normalization、rendering、I/O の依存方向を逆転させない。
- 1 機能の公開に registry、plugin system、汎用 request DSL を導入しない。

### 4.3 互換性

- 既存 import path、既定値、static SVG の意味を維持する。
- additive な field/symbol 追加を優先する。
- 挙動修正で既存利用者に影響する場合は、変更理由と移行方法を release note に明記する。
- `assemble_*` は互換層として維持し、先に削除・非公開化しない。
- 既定の `collinearity_anchor_mode="rbh"` は維持する。

### 4.4 コード量

- 既存 helper の re-export または小さな adapter を優先する。
- 新しい option class は、重複引数または raw config 操作を実際に減らす場合だけ追加する。
- 同一ロジックを CLI と API に複製しない。
- `docs/CODEBASE_REDUCTION_PLAN.md` の YAGNI/KISS 方針と矛盾する変更は採用しない。

## 5. 優先度付き所見（監査時）

この節の「現状」は監査スナップショット時点を表す。実装後の状態は「0. 実装結果」を参照する。

### P1: 公開引数 `collinearity_anchor_mode` が反映されない

現状:

- `DiagramOptions.collinearity_anchor_mode` は `all`、`one_to_one`、`rbh` 相当を受け取る型になっている。
- `assemble_linear_diagram_from_records` は値を検証した後、常に `rbh` を使う。
- `build_linear_diagram` も常に `rbh` を下位へ渡す。

目標:

- 指定された正規化済み mode を collinearity 構築へ渡す。
- 既定値は `rbh` のまま維持する。
- CLI が意図的に `rbh` 固定であるなら、CLI だけが明示的に `rbh` を渡す。

完了条件:

- `all`、`one_to_one`、`rbh` の各値が `build_orthogroup_collinearity_blocks(..., edge_mode=...)` まで届く。
- alias を含む正規化テストが通る。
- 「`all` を渡しても `rbh` になる」現行 characterization test を正しい期待値へ更新する。
- 既定値を省略した結果は改修前と同じになる。

### P1: 型付き config と override の合成でラベル用 DataFrame が失われる

現状:

- `apply_config_overrides(GbdrawConfig, ...)` は `dataclasses.asdict()` を使用する。
- `LabelsFilteringConfig.raw` に保持した `whitelist_df`、`qualifier_priority_df`、`label_override_df` が一段深い `raw` の下へ移動する。
- downstream は `cfg.labels.filtering.as_dict()` の top level を参照するため、これらの DataFrame を見つけられない。

推奨修正:

1. `GbdrawConfig` 全体の汎用 serializer は追加しない。
2. `apply_config_overrides` で `asdict()` 後の `labels.filtering` を、`deepcopy(config.labels.filtering.as_dict())` で置換する。
3. DataFrame の identity ではなく内容と downstream の抽出結果で検証する。

完了条件:

- typed config のみ、dict config のみ、typed config + override の結果が意味的に一致する。
- whitelist、qualifier priority、label override の各 DataFrame が保持される。
- override 対象外の filtering key が失われない。
- Circular/Linear の両方で label 選択結果を確認する。

### P1: `save_figure_to` が未生成ファイルの path を返す

現状:

- CairoSVG がない場合も、要求された PNG/PDF/EPS/PS の path を返す。
- 変換時に例外が発生しても log のみで、未生成 path を含む list を返す。
- library caller は戻り値だけでは成功を判断できない。

推奨契約:

- 明示的に要求された形式を生成できない場合、`ValidationError` または export 専用の `GbdrawError` を送出する。
- 成功時の戻り値には、存在を確認した path だけを含める。
- CLI の「CairoSVG がない形式を warning して skip」という既存挙動は `save_figure` 側に残してよい。
- Pyodide で browser conversion に委ねる経路は、ローカル file path を成功扱いしない。

完了条件:

- SVG-only、interactive SVG、各 binary format、複数形式、overwrite のテストがある。
- CairoSVG unavailable、converter exception、部分変換失敗の結果が明示的である。
- 戻り値の全 path が関数終了時に存在する。
- `render_to_bytes` と例外方針が整合する。

### P2: rich interactive SVG の公開契約が不完全

現状:

- `save_figure_to` は `InteractiveSvgContext` を型注釈で要求できる。
- `InteractiveSvgContext` と `enrich_svg` は `gbdraw.api` から export されていない。
- Circular/Linear CLI の context builder は private で、利用者が feature/match metadata を再現しにくい。
- `render_to_bytes(canvas, "interactive_svg")` は context を受け取れない。

目標 API:

```python
context = build_interactive_svg_context(
    records,
    selected_features_set=features,
    feature_table=feature_table,
    color_table=color_table,
    default_colors=default_colors,
    orthogroups=orthogroups,
)

data = render_to_bytes(canvas, "interactive_svg", interactive_context=context)
paths = save_figure_to(
    canvas,
    "interactive_svg",
    interactive_context=context,
)
```

設計方針:

- `InteractiveSvgContext` と `enrich_svg` を `gbdraw.api` から re-export する。
- context 構築処理の owner を CLI module から CLI 非依存 module へ移す。
- Circular/Linear の差は必要な metadata 入力で表現し、巨大な mode switch を作らない。
- `render_to_bytes` へ optional keyword-only `interactive_context` を追加する。

完了条件:

- feature popup と pairwise match popup に必要な metadata が公開 API だけで生成できる。
- static SVG の byte 列は変更しない。
- context なしの interactive SVG は従来どおり生成できる。
- CLI interactive SVG と API interactive SVG の主要 metadata schema が一致する。

### P2: CLI table 入力が公開 API から利用しにくい

対象:

- records table
- conservation table
- circular track table
- label whitelist
- qualifier priority
- label override table

改修方針:

1. `RecordsTable`、`ConservationTable`、`CircularTrackTable` と各 reader を `gbdraw.api.io` から薄く公開する。
2. label table reader も同じ namespace へ追加する。
3. `DiagramOptions` にラベル用 DataFrame の明示 field を追加する。
4. 既存の raw `config` 経路は互換性のため残す。
5. file path と DataFrame の両方が指定された場合は、他の API と同様に `ValidationError` にする。

候補 field:

```python
label_whitelist_table: DataFrame | None = None
label_whitelist_file: str | None = None
qualifier_priority_table: DataFrame | None = None
qualifier_priority_file: str | None = None
label_override_table: DataFrame | None = None
label_override_file: str | None = None
```

完了条件:

- CLI table と API reader が同じ validation を使う。
- path 解決、row-scoped selector、region、reverse flag、multi-record position の結果が一致する。
- 表を使うために `gbdraw.io.*` や `gbdraw.config.*` を import する必要がない。

### P2: multi-record Circular の高位 API がない

現状:

- `assemble_circular_diagram_from_records` では全機能を使える。
- `build_circular_diagram` は単一 `SeqRecord` のみを受け取る。
- `DiagramOptions` には multi-record size mode、gap、position がない。

改修方針:

- 既存関数の入力型を曖昧に拡張せず、`build_circular_multi_diagram(records, *, options, layout)` を追加する。
- multi-record 固有値だけを持つ小さな `CircularMultiRecordOptions` を追加する。
- `assemble_circular_diagram_from_records` は互換層・高度利用者向けとして維持する。

候補 field:

```python
multi_record_size_mode: Literal["linear", "auto", "equal", "sqrt"] = "auto"
multi_record_min_radius_ratio: float = 0.55
multi_record_column_gap_ratio: float = 0.10
multi_record_row_gap_ratio: float = 0.05
multi_record_positions: Sequence[str] | None = None
```

完了条件:

- 現行 CLI の multi-record canvas と同じ配置を公開 API で再現できる。
- 単一 record の `build_circular_diagram` は変更しない。
- 共有 legend、plot title、record definition、track geometry metadata が維持される。

### P3: session は内部機能で、library workflow が未定義

現状:

- `load_session`、`validate_session`、`session_to_cli_args` 等は `gbdraw.session_io` にある。
- `session_to_cli_args` は argparse/CLI invocation を前提とし、公開 Python API としては抽象度が低い。
- LOSAT cache、embedded file、GUI state を含むため、単純な re-export だけでは不十分である。

実施前の設計ゲート:

次のどちらを公開契約とするか、1ページ以内の ADR で決定する。

1. session を materialize し、公開 options/input model へ変換する。
2. `render_session(...)` が内部 orchestration を呼び、構造化結果を返す。

推奨は 1 とする。CLI 引数へ変換する経路を公開 API の中心にすると、CLI option 名が library contract へ漏れるためである。

段階的な実施:

1. `SessionDocument` 相当の検証済み payload と embedded file materialization を公開する。
2. Circular/Linear の input/options へ変換する pure function を用意する。
3. LOSAT cache は optional metadata として保持し、cache がなくても描画可能にする。
4. session 保存は、入力・options・生成結果を受け取る builder として公開する。

完了条件:

- GUI session と CLI sidecar session の両方を公開 API だけで検証・再実行できる。
- 一時ファイルの lifetime と cleanup が context manager または明示的 result object で管理される。
- session version 不一致と壊れた embedded file を型のある例外で識別できる。
- CLI module を `gbdraw.api` から import しない。

### P3: ドキュメントが実装能力を表していない

現状:

- `docs/PYTHON_API.md` は Circular の最小例が中心。
- Linear は短い説明のみ。
- depth、conservation、protein comparison、multi-record、track slot、interactive SVG、table/session の実行例がない。

追加する実行可能な例:

1. Circular single record。
2. Circular multi-record canvas。
3. Linear nucleotide BLAST comparison。
4. Protein pairwise/orthogroup/collinear のいずれか 1 つと、precomputed comparison。
5. depth、conservation、custom track slot。
6. label/color/visibility table。
7. rich interactive SVG。
8. session load/render/save。
9. in-memory bytes と binary export の失敗処理。

完了条件:

- 各 code block を pytest から実行する。
- optional dependency/external binary が必要な例は、fixture または precomputed data で CI 実行可能にする。
- `gbdraw.api` 外の import が必要な例は、理由を明記した低位拡張例に限定する。

## 6. 実施フェーズ

### Phase 0: 契約と失敗モードの固定

状態: 完了。

目的:

- 修正前の正常系と、修正すべき誤動作を characterization test で区別する。

実施内容:

1. `gbdraw.api.__all__`、signature、dataclass field/default の snapshot を維持する。
2. CLI → assembler の代表的な引数転送テストを整理する。
3. 次の failing regression test を先に追加する。
   - anchor mode が指定値まで届く。
   - typed config + override で label DataFrame が残る。
   - export failure が未生成 path を返さない。
4. interactive context あり/なしの metadata schema を固定する。
5. full test、lint、reference SVG の改修前基準値を記録する。

完了条件:

- 正常系の baseline test が通る。
- 3 件の既知不具合を表す test が修正前に意図どおり失敗する。

### Phase 1: correctness 修正

状態: 完了。

順序:

1. `collinearity_anchor_mode` を反映する。
2. typed config の filtering raw data を保持する。
3. `save_figure_to` の失敗・戻り値契約を修正する。

各項目を独立 PR にする。SVG geometry を変更する修正ではないため、reference SVG に差分が出た場合は原因を調査し、機械的に更新しない。

完了条件:

- P1 の完了条件をすべて満たす。
- public contract snapshot の変更は意図した signature/exception contract に限定される。

### Phase 2: interactive と table 入力の公開

状態: 完了。

順序:

1. `InteractiveSvgContext` と `enrich_svg` の re-export。
2. CLI 非依存の interactive context builder。
3. `render_to_bytes(..., interactive_context=...)`。
4. CLI table model/reader の re-export。
5. label DataFrame/file の明示 options。

完了条件:

- rich interactive SVG と table-driven input の利用例が `gbdraw.api` import のみで動く。
- CLI と API が同じ parser/builder を共有する。

### Phase 3: 高位 API の不足を補う

状態: 高位 API 追加は完了。multi-record 固有 option だけを分離し、既存
`DiagramOptions` の全面分割は行わなかった。field 利用状況の監査は後続作業とする。

順序:

1. `CircularMultiRecordOptions`。
2. `build_circular_multi_diagram`。
3. 現行 `DiagramOptions` の field を用途別に一覧化し、重複・dead field を監査する。
4. option class 分割は、呼び出しコードと実装コードの合計が明確に簡潔になる場合のみ実施する。

この phase では既存 `DiagramOptions` を廃止しない。

### Phase 4: session bridge

状態: 設計ゲート完了、実装停止。判断と再開条件は
[`ADR_PYTHON_SESSION_API.md`](ADR_PYTHON_SESSION_API.md) に記録した。

順序:

1. ADR で公開モデルを決定する。
2. session validation/materialization の公開。
3. session → public input/options 変換。
4. render と session 保存。
5. LOSAT cache の optional round-trip。

session 変換が CLI parser への依存を避けられない場合は、無理に公開せず phase を停止し、依存境界の整理を別計画にする。

### Phase 5: ドキュメントと安定化

状態: API 文書、実行 test、公開契約 snapshot は完了。repository 全体の既存 lint debt と
release note は後続作業とする。

実施内容:

1. `docs/PYTHON_API.md` を capability matrix と実行可能 recipe で拡張する。
2. CLI option と Python API の対応表を追加する。
3. exception、optional dependency、external binary、version pinning を説明する。
4. public API contract test を CI の必須 test にする。
5. beta 期間終了前に deprecation/compatibility 方針を確定する。

## 7. PR 分割案

| PR | 内容 | 規模 | 主なリスク |
|---|---|---:|---|
| 1 | anchor mode の転送修正と回帰テスト | S | collinearity 出力差分 |
| 2 | typed config filtering の lossless override | S | DataFrame deepcopy と比較 |
| 3 | export の strict failure contract | S–M | 既存 caller の期待値変更 |
| 4 | interactive 型の re-export と bytes API | S | public contract 追加 |
| 5 | interactive context builder の共通化 | M | CLI/API metadata 差分 |
| 6 | table reader と label table options の公開 | M | path/selector validation |
| 7 | Circular multi-record 高位 API | M | option forwarding 漏れ |
| 8 | session 公開モデル ADR と materialization | M | lifetime、security、互換性 |
| 9 | session render/save bridge | L | CLI/Web/Python の依存境界 |
| 10 | Python API docs と実行テスト | M | optional dependency の CI 差 |

PR 1～3 を release-blocking correctness、PR 4～7 を API completeness、PR 8～10 を workflow completeness とする。

実装結果として PR 1～7 相当と PR 10 相当は完了した。PR 8 は ADR の設計ゲートまで、
PR 9 は ADR の再開条件を満たさないため実施していない。この表は今後も論理的な review
単位として残し、実際の commit 数を規定しない。

## 8. テスト計画

### 8.1 targeted test

最低限、各 PR で関連する次の test を実行する。

```bash
python -m pytest \
  tests/test_api_library_usage.py \
  tests/test_public_contract.py \
  tests/test_comparisons.py \
  tests/test_feature_shapes.py \
  tests/test_gc_content_percent.py \
  tests/test_circular_conservation.py \
  tests/test_circular_track_slots.py \
  tests/test_linear_track_slots.py \
  tests/test_depth_track.py \
  tests/test_interactive_svg_cli_format.py \
  tests/test_linear_selectors.py \
  tests/test_cli_tables.py \
  -q -m "not slow"
```

### 8.2 full gate

```bash
python -m pytest tests/ -v -m "not slow"
ruff check gbdraw/
```

出力または配置に影響する変更では、該当する reference SVG test と browser test も実行する。

`TestGenerateReferences` は `--update-reference-outputs` を指定した場合だけ書き込む。
通常の full gate は tracked reference SVG を読み取るだけである。

### 8.3 追加する契約テスト

- 全 public symbol の import。
- public function signature と dataclass default。
- CLI と API の同一入力に対する主要 SVG group/metadata の一致。
- file input と DataFrame input の同値性。
- missing dependency、invalid table、invalid selector、converter failure の例外型。
- session version と embedded file validation。
- documentation code block の実行。

## 9. 互換性とリリース方針

### 9.1 additive change

次は通常の minor/beta 更新で追加できる。

- symbol の re-export。
- optional keyword-only 引数。
- option dataclass の optional field。
- 新しい high-level builder。
- 新しい例外 subclass。ただし既存 `GbdrawError` で捕捉可能にする。

### 9.2 behavior correction

次は release note へ明記する。

- `collinearity_anchor_mode` が初めて指定どおり反映される。
- `save_figure_to` が missing converter を例外にする。
- 戻り値から未生成 path が除外される。

### 9.3 deprecation

- 低位 `assemble_*` の削除は本計画に含めない。
- dead option が見つかった場合も、即時削除せず warning → documented replacement → major release の順にする。
- 内部 module import の互換性は保証しないが、公開代替を提供してから文書から除く。

## 10. リスクと対策

| リスク | 対策 |
|---|---|
| API を増やして CLI と二重実装になる | owner を共有し、API は adapter/re-export に限定する |
| option class が増えてさらに複雑になる | multi-record と label tables の不足分だけ追加し、全面分割を急がない |
| export の strict 化が既存 pipeline を壊す | `save_figure` の CLI 挙動を維持し、`save_figure_to` の変更を release note に明記する |
| anchor mode 修正で図が変わる | default `rbh` を維持し、非 default 指定時だけ変更する |
| session API が CLI 引数へ密結合する | ADR と dependency gate を設け、条件を満たせなければ phase を停止する |
| DataFrame を config に保持することで deepcopy cost が増える | config 構築時のみ処理し、render loop 内で複製しない |
| interactive metadata schema が Web と乖離する | CLI/API/Web が共通 builder または schema test を使う |

## 11. 完了の定義

本計画は、次をすべて満たしたとき完了とする。

1. 公開済み option が黙って無視されない。
2. typed config + override で label/filter DataFrame が保持される。
3. export API の戻り値が実際の生成結果と一致する。
4. 公開 API だけで rich interactive SVG を生成できる。
5. records/conservation/track/label table を公開 API から読める。
6. Circular multi-record canvas を high-level builder から作れる。
7. session を CLI parser に依存せず検証・再実行できる、または非公開とする合理的な ADR がある。
8. static diagram の既存既定出力に意図しない差分がない。
9. public contract、targeted test、full non-slow test、lint が通る。
10. 高度機能の主要 recipe が文書から実行され、CI で検証される。

2026-07-15 最終判定:

- 1～10 は documented mode の範囲で完了した。
- public contract、targeted test、full non-slow test、repository 全体 lint は成功した。
- wrong-mode の non-default field と単一 Circular depth input の fail-fast validation は
  [`PYTHON_API_VALIDATION_PLAN.md`](PYTHON_API_VALIDATION_PLAN.md) で完了した。

## 12. 初回実装単位（完了）

最初の実装単位は次の 3 PR とする。

1. `collinearity_anchor_mode` の転送修正。
2. `apply_config_overrides` の filtering raw data 保持。
3. `save_figure_to` の strict failure と実在 path 契約。

この 3 件は、API surface を大きく増やさず、現行利用者が結果を信頼できる状態にする。完了後に interactive と table の公開へ進む。

3 件とも完了し、続く interactive、table、multi-record の公開まで実装した。

## 13. 次の作業

詳細な実施計画は
[`PYTHON_API_FOLLOWUP_PLAN.md`](PYTHON_API_FOLLOWUP_PLAN.md) に記録した。

優先順は次のとおりとする。

1. **完了:** reference SVG 生成を明示 opt-in にし、通常の full suite を comparison-only にする。
   生成 test が回帰 fixture を先に更新して差分を隠さないことを確認する。
2. **完了:** `0.14.0b0` の release note に anchor mode と strict export の behavior correction、
   新しい public symbol、session 非公開判断を記載する。
3. **完了:** repository 全体 Ruff の既存 85 件を、API 改修とは独立した cleanup として baseline 化し、
   段階的に解消する。
4. **完了:** `DiagramOptions` 70 field の用途・利用箇所を監査する。分割は API と実装の合計コード量が
   明確に減る場合だけ提案する。
5. **E0 完了、公開ゲート継続:** version 27～30 の compatibility matrix で lossless typed
   conversion が不可能と判定した。session-independent request model と version 31 canonical
   payload を先に設計し、それまでは公開 symbol を追加しない。

監査後の validation backlog と E0 は
[`PYTHON_API_VALIDATION_PLAN.md`](PYTHON_API_VALIDATION_PLAN.md) で完了した。次の request/schema
実装は [`PYTHON_SESSION_CANONICAL_REQUEST_PLAN.md`](PYTHON_SESSION_CANONICAL_REQUEST_PLAN.md)
に記録した。
