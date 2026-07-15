# Python Session Canonical Request Plan

- 作成日: 2026-07-15
- 対象: Workstream E1～E6、session version 31 の canonical typed request
- 前提: [`Python Session Compatibility Matrix`](PYTHON_SESSION_COMPATIBILITY_MATRIX.md)
- 状態: Phase R1～R2 完了。Phase R3 canonical payload codec が次の実装

## 0. 目的と境界

version 27～30 の CLI argument replay を公開 API にせず、Circular/Linear の入力、diagram
option、comparison、output policy を CLI 非依存の型で表現する。その型から render できることを
先に証明し、同じ型を version 31 の `renderRequest` に保存する。

この計画では次を守る。

- legacy `session_to_cli_args()` は version 27～30 の internal compatibility path として残す。
- request model は session JSON、CLI option 名、argument index を知らない。
- schema codec と embedded resource materialization は request model から分離する。
- `gbdraw.api.__all__` と public contract fixture は R1～R5 では変更しない。
- version 31 writer は Python/Web の双方が読める段階まで有効化しない。

## 1. target architecture

```text
in-memory/file caller ──> typed request ──> request renderer ──> render result
                              ^
                              |
v31 renderRequest ──> resource materialization + schema codec

v27–30 session ──> existing internal CLI replay（public conversion なし）
```

責務は次の三層に分ける。

| Layer | Owner | CLI/session dependency |
|---|---|---|
| Runtime request/value objects | `gbdraw.api.requests` | なし |
| Canonical JSON codec | new internal session request codec | session schema のみ。CLI parser なし |
| Legacy replay | `gbdraw.session_io` | 現行 CLI args/bindings を維持 |

## 2. Phase R1: session-independent request model

状態: 完了。`gbdraw.api.requests` に incubation model を追加したが、root
`gbdraw.api.__all__` からは export していない。

実施結果:

- `GenBankInputSource`、`GffFastaInputSource`、`InMemoryRecordSource`。
- `RecordInput`、`RecordPresentation`。
- `CircularDiagramRequest`、`LinearDiagramRequest`、`RenderOutputRequest`。
- materialized path、selector/region owner、double reverse、grid placement、format/output path、
  mode-specific `DiagramOptions` を fail-fast validation する。
- mode-specific option validation は `gbdraw.api.options` の単一 owner を builder/request で
  再利用し、field list を複製していない。
- request module の CLI/session owner import prohibition を含む 30 test が成功した。
- API/request/public-contract targeted suite は 245 test が成功した。
- `ruff check gbdraw/` は成功した。

最初の change では JSON/session integration を行わず、型と validation だけを追加する。

### 2.1 Record input

候補型:

- `GenBankRecordInput`: materialized `Path`、optional `RecordSelector`、optional `RegionSpec`。
- `GffFastaRecordInput`: GFF3 `Path` と FASTA `Path`、同じ selector/region metadata。
- `RecordPresentation`: label、subtitle、reverse complement、optional grid placement。
- `RecordInput`: 上記 source と presentation を結合した immutable value object。

source kind は tagged union で表し、空 path、GFF/FASTA の片側欠落、負の grid position、
selector/region の矛盾を `ValidationError` にする。records table はこの型の sequence へ
normalize し、table row 自体を runtime request contract にしない。

### 2.2 Mode-specific request

- `CircularDiagramRequest`
  - 1 件以上の `RecordInput`。
  - `DiagramOptions`。
  - multi-record の場合だけ `CircularMultiRecordOptions`。
- `LinearDiagramRequest`
  - 1 件以上の `RecordInput`。
  - `DiagramOptions`。
  - comparison input は既存の typed `DiagramOptions` field を使用する。
- `RenderOutputRequest`
  - output prefix、format token、overwrite policy、interactive metadata policy。

Circular と Linear は別 dataclass にする。mode-specific option validation は既存 high-level
builder と同じ owner を再利用し、同じ field list を複製しない。

### 2.3 Comparison input

R1 では comparison 専用 dataclass を追加しない。`DiagramOptions` は nucleotide BLAST、
precomputed protein comparison、orthogroup、collinearity の typed field と validation owner を
すでに持つため、同じ値を request に複製すると二つの precedence が生じる。

canonical JSON codec では少なくとも次を別 kind/resource として区別する。

- nucleotide BLAST result。
- precomputed protein comparison。
- generated protein/orthogroup/collinearity request。
- collinearity block/result。

runtime `LinearDiagramRequest` はこれらを `DiagramOptions` から受け取り、codec が schema 上の
`comparisons` との相互変換を所有する。専用 runtime model は、R3 で code/validation の重複を
減らせる証拠がある場合だけ追加する。LOSAT process execution と cache serialization は後
phase に分離する。

### 2.4 R1 gate

- `gbdraw.api.requests` は `gbdraw.circular`、`gbdraw.linear`、CLI parser、`session_io` を
  import しない。
- runtime model に CLI option 名、argument index、session slot がない。
- GenBank、GFF3+FASTA、selector、region、reverse、placement の validation test がある。
- Circular/Linear wrong-mode request が assembly 前に失敗する。
- `gbdraw.api.__all__` と `tests/fixtures/public_contract.json` は不変。

## 3. Phase R2: request normalization and rendering

状態: 完了。`gbdraw.api.request_render` に internal adapter を追加し、root public export は
行っていない。

実施結果:

- file/in-memory source を copy/load し、selector、reverse、region、label/subtitle の順で
  1 `RecordInput` → 1 `SeqRecord` に正規化する。
- GFF3+FASTA は既存 candidate-feature resolver を再利用し、selected feature と
  color/visibility table 由来の feature type を保持する。
- Circular single/multi と Linear を既存 high-level builder へ渡す。
- canonical grid row/column を既存 `multi_record_positions` の row/order contract に変換する。
- output/interactive policy を strict `save_figure_to()` と metadata builder に渡し、実在する
  output path だけを `RequestRenderResult` に返す。
- request/render/API/session/public-contract targeted suite は 281 test が成功した。
- broad Python API targeted suite は 547 test が成功した。
- real Circular request → SVG smoke test と `ruff check gbdraw/` が成功した。

typed request から既存 public I/O/build/render owner を呼ぶ adapter を作る。

1. record source を `load_gbks()` / `load_gff_fasta()` で読む。
2. selector、region、reverse を既存 helper で適用する。
3. record label/subtitle/placement を mode-specific builder input に normalize する。
4. `build_circular_diagram` / `build_circular_multi_diagram` / `build_linear_diagram` を呼ぶ。
5. `RenderOutputRequest` を `save_figure_to` / interactive SVG owner に渡す。
6. 実在 path、warnings、mode、metadata を持つ typed render result を返す。

file read と output write の side effect を adapter 境界に集める。request validation と
normalization は pure function として個別 test できるようにする。

## 4. Phase R3: canonical payload schema 1

version 31 の top-level `renderRequest` は次を持つ。

```json
{
  "schema": 1,
  "mode": "linear",
  "records": [],
  "diagramOptions": {},
  "layout": {},
  "comparisons": [],
  "output": {}
}
```

embedded file は CLI argument index や `files.*` slot path ではなく、stable resource ID で
参照する。session document には resource ID → embedded entry の canonical mapping を追加する。
runtime request は materialization 後の `Path` だけを受け取り、resource ID を保持しない。

schema rule:

- unknown schema/mode/kind は typed conversion error。
- unknown optional field は forward compatibility policy を明示して扱う。
- required field の欠落を default 推測しない。
- `DiagramOptions` の enum/number/default は public dataclass contract と一致させる。
- DataFrame input は canonical TSV/resource に変換するか、保存不能として明示 error にする。
- output destination directory は session から信頼せず、replay caller が指定する。
- `results`、`editorState`、LOSAT cache は request 本体と分離する。

codec は `request -> JSON-compatible payload -> request` の round-trip test を両 mode で持つ。

## 5. Phase R4: validated document and materialization

候補 lifecycle:

```python
document = load_session_document(path)
with materialize_session(document) as materialized:
    request = session_to_request(materialized)
```

- materialization context が temporary directory と全 resource path を所有する。
- basename sanitization、directory containment、base64、declared size、depth codec validation を
  既存 owner から再利用する。
- partial failure、duplicate resource ID/name、missing reference、cleanup failure を test する。
- request path は context exit 後に利用不能であることを docs/type docstring で明記する。
- version 27～30 に canonical request がなければ public conversion は明示 error とする。

## 6. Phase R5: session version 31 and Web migration

version bump は次を一つの change として行う。

1. Python `CURRENT_SESSION_VERSION` と Web `SESSION_VERSION` を 31 にする。
2. Python CLI sidecar writer が canonical `renderRequest` と resources を生成する。
3. Web save が current editor state から同じ canonical payload を生成する。
4. Web load は canonical payload を editor projection へ反映する。
5. `renderRequest` がある場合はそれを authority とし、`cliInvocation` や GUI config との
   silent merge を行わない。
6. version 27～30 の current internal replay/read compatibility を維持する。
7. Python/Web schema fixture と browser round-trip test を同じ PR で更新する。

version 27～30 を自動で canonical request に書き換える migration は、各 capability が
lossless であると証明できた path に限定する。それ以外は legacy document のまま開き、
保存前に unsupported state を明示する。

## 7. Phase R6: public session bridge

R1～R5 完了後にのみ public export を行う。

- document load/inspect API。
- materialization context manager。
- mode-specific `session_to_request()`。
- request render API と structured result。
- canonical session save API。
- typed format/version/resource/conversion/render errors。

公開時に `gbdraw.api.__all__`、contract fixture、Python API guide、release note、session ADR を
同時に更新する。legacy CLI argument replay helper は export しない。

## 8. verification gates

各 phase:

```bash
python -m pytest \
  tests/test_api_requests.py \
  tests/test_api_request_render.py \
  tests/test_session_io.py \
  tests/test_api_library_usage.py \
  tests/test_public_contract.py \
  -q
ruff check gbdraw/
git diff --exit-code -- tests/reference_outputs/
```

R2 以降は Circular/Linear の representative SVG comparison、R5 は Web session tests と
Playwright、R6 は full non-slow suite を追加する。reference output の更新は intentional
geometry change がない限り行わない。

## 9. stop conditions

次のいずれかが発生した phase では public export に進まない。

- runtime request が CLI option/argument index/session slot を必要とする。
- GUI current state から canonical request を lossless に作れない。
- version 31 Python/Web writer の payload が一致しない。
- temporary resource の cleanup owner が一意でない。
- save/load/render round-trip が mode または input kind で欠落する。

[Session ADR](ADR_PYTHON_SESSION_API.md) | [Compatibility matrix](PYTHON_SESSION_COMPATIBILITY_MATRIX.md) | [Follow-up plan](PYTHON_API_FOLLOWUP_PLAN.md)
