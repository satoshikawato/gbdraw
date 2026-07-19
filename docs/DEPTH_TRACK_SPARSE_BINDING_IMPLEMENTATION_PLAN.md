# Sparse Depth Track Binding Implementation Plan

- 作成日: 2026-07-19
- 対象バージョン: `0.14.0b0` 以降
- 設計スナップショット: branch `multi_records_per_line`, HEAD `5eff948` の作業ツリー
- 状態: Phase 0–6実装済み
- 対象: Linear/Circularで共有するDepth入力の正規化、論理系列binding、Linear Custom Track Slots、Web/session
- 関連計画: `docs/LINEAR_TRACK_OCCUPANCY_LAYOUT_IMPLEMENTATION_PLAN.md`

## 0. 実装状況（2026-07-19）

| Phase | 状態 | 現在の範囲 |
|---|---|---|
| 0. Characterization | 実装済み | sparse normalization、logical aggregation、Linear/Circular binding、CLI placeholder、Web/session、browserのfocused regression testを追加した。 |
| 1. Logical identity | 実装済み | `DepthTrackSpec`と`DepthTrackData`へ`track_index`を伝搬し、count、height、shared axis、representative selectionをlogical index基準へ変更した。 |
| 2. Binding and aggregation | 実装済み | Linear default/customとCircular multi-recordでrecord-local lookupを使用する。欠損cellはpaintせず、非整数・負数・globalに存在しないslot indexをerrorにする。custom凡例もenabled slotを正本とする。 |
| 3. Web state | 実装済み | Linear/Circular Depthをglobal series column x record rowとして扱い、cell clearとseries removalを分離した。argument生成はrecord数分のplaceholderを保持し、全recordで空のseriesをerrorにする。 |
| 4. Session normalization | 実装済み | settings/canonical projectionで`null` cellを保持し、saved row、metadata、slot参照から求めたlogical widthまでrowをpaddingする。Depth shape自体は維持し、後続のLinear slot schema v2移行でsession versionを33へ更新した。complete sessionをpreflight後にだけlive stateへcommitする。 |
| 5. Occupancy planner | 実装済み | Default/Customを共通のper-record plannerへ接続し、欠損Depthをpaintなし・reserveありとして扱う。recordごとのfeature/numeric occupancyとcomparison exclusion edgeを同じplanから解決する。 |
| 6. Cleanup/docs/deployment | 実装済み | positional/legacy offset ownerを削除し、Tutorial、CLI Reference、FAQ、release notes、Web helpを更新した。browser wheelはpackaging/browser検証用に再生成し、cache-bust更新だけをdeploy operationとして分離した。 |

この更新時点でsparse binding、session state、Linear occupancy plannerへの接続まで完了した。
dense inputのlogical bindingとCircular radial layoutのgeometryは変更対象外のままとする。

## 1. 結論

この修正では、Depthの`track_index`をrecord内listの位置として扱わない。
`track_index`はdiagram全体で共有する論理Depth系列の識別子とし、各recordにその系列のデータが
存在するかどうかを別の状態として扱う。

```text
global logical depth tracks
          |
          +-- track_index=0  Sample A
          +-- track_index=1  Sample B
          |
          v
record x logical-track sparse binding
          |
          +-- record 0: {0: A0}
          +-- record 1: {1: B1}
          +-- record 2: {}
          |
          v
per-record data preparation
          |
          v
default/custom track rendering from the same binding contract
```

recordに指定系列のDepthデータがない場合は正常な欠損として扱う。

- Depth groupをpaintしない。
- 別系列や別recordのデータを代用しない。
- coverage 0として描かない。
- diagram全体のslot構造と予約領域は維持する。
- globalに存在しない`track_index`だけをvalidation errorとする。

初回修正では公開CLI/API/session形式を変更しない。既存の列順を論理identityとして内部まで明示的に
伝搬させる。永続UUIDや新しい公開`track_id`は導入しない。

Depth identity修正とLinear全体の縦配置変更は別phaseにする。先にdense出力を変えずにbindingを直し、
その後で既存のLinear occupancy計画へ接続する。

## 2. 現状の問題

### 2.1 入力契約は疎なmatrixを許す

Depthの公開入力はrecord-major matrixであり、cellは`DataFrame | None`または`str | None`である。
CLIとWebも、recordに対応するDepth TSVがない場合を空値として保持する。

```text
                 logical track 0   logical track 1
record 0               A                None
record 1              None               B
```

この入力は、record 0だけにA、record 1だけにBを描く意味を持つ。

### 2.2 正規化が論理列を圧縮する

`gbdraw/analysis/depth_tracks.py::_tracks_from_record_major_values()`は、各論理列を走査した後、
`None`をskipして存在する`DepthTrackSpec`だけをrecord rowへ`append()`する。

その結果、上のmatrixは次のlistへ変わる。

```text
record 0: [A]
record 1: [B]
```

この時点でAが論理track 0、Bが論理track 1だった情報をlist位置から復元できない。

### 2.3 後段が圧縮後のlist位置をidentityとして使う

現在は次の処理が`enumerate(row)`または`row[track_index]`に依存する。

- `depth_track_count()`
- `depth_track_heights()`
- `depth_track_data_count()`
- per-track shared-axis最大値
- legend representativeの選択
- Linear default Depth配置
- Linear custom slotからDepth dataへのbinding
- Circularで共有するDepth count/data helper

そのため、欠損rowが空なら例外になり、前方に穴があれば別系列を誤って描く可能性がある。
単純に例外を`continue`へ変えても系列取り違えは残る。

### 2.4 DefaultとCustomで欠損処理が異なる

Linear default経路はrecord内に存在するDepth dataだけを列挙するため、空rowでは何も描かない。
Custom Track Slots経路はdiagram共通のslotを全recordへ適用し、record-local listを
`shared_depth_tracks[track_index]`で参照する。

したがって、部分的なDepth入力はdefaultでは成功し、Custom Track Slotsを有効にすると失敗する。
DepthをFeatureより上へ移動する操作はCustom経路を有効にするため、上下順が原因のように見えるが、
同じ欠損入力はDepthがAbove/Belowのどちらでも失敗する。

### 2.5 Webのseries操作とrecord-local file操作が混在する

Webでは`adv.depth_tracks`とDepth slotがdiagram-globalである一方、`linearSeqs[*].depth`はrecord-localである。
現在のadd/remove処理は一つのrecordの配列をspliceしながら、global metadataとglobal slot indexも更新する。
別recordの同じ列が同期して変更されない場合、保存された配列位置とmetadataの対応がずれる。

## 3. 設計上の不変条件

### 3.1 Logical identity invariant

1. 各`DepthTrackSpec`と`DepthTrackData`は明示的な非負整数`track_index`を持つ。`bool`、小数、非数文字列、負数は拒否し、暗黙の`int()`丸めを行わない。
2. record row内のlist位置をlogical identityとして使用しない。
3. 同一record row内で同じ`track_index`を重複させない。
4. 同一logical trackのlabel、color、height、tick設定は全recordで共通とする。
5. 縦方向slot orderを変更しても`track_index`とDepth sourceのbindingは変わらない。

### 3.2 Missing-data invariant

1. per-record cellの欠損は正常値である。
2. 欠損recordでは該当Depth group、axis、tick、pathを生成しない。
3. 欠損を別recordのsourceで暗黙補完しない。
4. 欠損をdepth 0のseriesへ変換しない。
5. global slotのreserve bandは維持し、record間のstack alignmentを崩さない。

### 3.3 Validation invariant

1. 全recordで空のlogical track列は入力正規化時にerrorとする。
2. slotがglobal logical track範囲外を参照した場合は描画前にerrorとする。
3. globalには存在するが特定recordだけ欠損している場合はerrorにしない。
4. validation messageはslot ID、logical track index、必要な場合はrecord ID/indexを含める。

### 3.4 Compatibility invariant

1. dense Depth入力のSVG geometry、DOM ID、color、legend、axisを変更しない。
2. legacy `depth_table`、`depth_file`、`depth_tables`、`depth_files`を維持する。
3. `--depth_track`のtrack-major CLI構造を維持する。
4. `DiagramOptions.depth_track_tables/files`のrecord-major matrixを維持する。
5. identity修正だけではsession versionとcanonical request schemaを上げない。
6. CircularのDepth identity処理は同時に修正するが、Circular layoutは変更しない。

## 4. 推奨内部モデル

### 4.1 DepthTrackSpecとDepthTrackData

既存dataclassへ`track_index`を追加する。

```python
@dataclass(frozen=True)
class DepthTrackSpec:
    id: str
    label: str
    table: DataFrame
    track_index: int = 0
    fill_color: str | None = None
    height: float | None = None
    large_tick_interval: float | None = None
    small_tick_interval: float | None = None
    tick_font_size: float | None = None


@dataclass(frozen=True)
class DepthTrackData:
    id: str
    label: str
    df: DataFrame
    config: DepthConfigurator
    track_index: int = 0
    height: float | None = None
```

default値は既存の内部生成やtest fixtureとの互換を保つために置く。正規化されたmulti-track入力では必ず
明示値を設定する。

### 4.2 Record row

record rowは存在するdataだけを持つ疎なlistのままでよい。ただしlookup前にlogical index mapへ変換する。

```python
def index_depth_track_row(
    row: Sequence[DepthTrackData],
) -> dict[int, DepthTrackData]:
    ...
```

またはSpec/Dataそれぞれに型付きhelperを用意する。汎用frameworkや新しい継承階層は作らない。

### 4.3 Diagram-global helpers

次のhelperをlogical index基準にする。

```text
depth_track_count(rows)
  = max(track.track_index for all present tracks) + 1

depth_track_heights(rows)
  = first configured height for each logical index

representative_depth_tracks(rows)
  = first present data for each logical index, ordered by index
```

正規化前のraw matrixで全recordが空の列を検出するため、global track列が完全に欠落した状態は
normalized modelへ入れない。

### 4.4 永続track IDを初回に導入しない理由

現行のCLI、Python API、Web state、sessionはすべて系列順を公開契約としている。今回必要なのは、
その列番号をrecord-local compactionから独立させることである。

UUIDまたは任意文字列IDを追加すると、次も同時に変更する必要がある。

- CLI slot syntax
- Python public options
- Web upload state
- settings/session migration
- canonical request schema
- gallery session

列の挿入・削除をglobal operationとして実装すれば、整数logical indexで今回の欠陥を根治できる。
将来、系列順とは独立した永続identityが必要になった時点で別ADRとして検討する。

## 5. 実装phase

### Phase 0: Characterization and failing tests

描画コードを変更する前に、疎なDepth bindingをtestで固定する。

1. `[[A], []]`: 後続recordがtrack 0欠損。
2. `[[], [A]]`: 先頭recordがtrack 0欠損。
3. `[[A, None], [None, B]]`: diagonal sparse matrix。
4. `[[A, B], [A, None], [None, None]]`: ragged rowとempty record。
5. `[[A, B]]`を複数recordへshared expansion。
6. dense 2 records x 2 tracks。
7. 全recordで空のlogical column。
8. slotがglobal count以上のindexを参照。

各sparse caseを次のlayoutで検証する。

- Depth、Features、Axis: Depthが同じsideのFeatureより外側。
- Depth、Axis、Features: DepthがAbove、FeatureがBelow。
- Features on Axis、Depth: DepthがBelow。
- slot orderを反転してもsource bindingが変わらない。

このphaseでは`tests/reference_outputs/`を更新しない。

### Phase 1: Logical identity propagation

`gbdraw/analysis/depth_tracks.py`を中心に変更する。

1. `_tracks_from_record_major_values()`で、`None`をskipする前の列番号を`DepthTrackSpec.track_index`へ保存する。
2. `build_depth_track_dataframes()`で`track_index`を`DepthTrackData`へ伝搬する。
3. `depth_track_count()`と`depth_track_data_count()`を明示index基準にする。
4. `depth_track_heights()`を明示index基準にする。
5. shared-axis最大値をlogical indexごとに集約する。
6. row内index重複と負値をvalidationする。
7. Circularを含む全`DepthTrackData`生成箇所へ明示indexを伝搬する。

Phase 1完了時点ではrendererの配置式を変更しない。

### Phase 2: Binding and aggregation cutover

Linear assemblyでrecordごとの`depth_by_index`を一度だけ作る。

Custom Track Slotsでは次の規則を使う。

```python
depth_track = depth_by_index.get(track_index)
if depth_track is None:
    continue
```

この`continue`は、Phase 1でidentityが保存された後にだけ導入する。現在の圧縮listに対する
range checkを単純に削除してはならない。

Default経路でも`enumerate(shared_depth_tracks)`をlogical identityとして使わず、
`depth_track.track_index`をDepth offset選択へ渡す。

さらに次を切り替える。

- APIのslot index validation: diagram-global logical countを使用する。
- legend: 全rowからlogical trackごとの代表dataを集める。custom stackではenabled Depth slotだけをslot順に出し、`legend_label`を優先する。
- height: logical indexごとに解決する。
- shared axis: logical indexごとに同じmaxを共有する。
- geometry metadata: slot indexとdata logical indexを混同しない。

Phase 2完了時点で、今回の例外とsparse diagonalの系列取り違えを解消する。

### Phase 3: Web state and operation semantics

WebではLinear/CircularともDepth stateを「global series columns x record rows」として扱う。

#### Add series

- `adv.depth_tracks`へmetadata列を追加する。
- 全`linearSeqs[*].depth`へ同じindexの`null` cellを追加する。
- managed Depth slotを一つ追加する。

#### Clear one record source

- 対象cellを`null`にする。
- rowをspliceしない。
- metadataとslot indexを変更しない。

#### Remove series

- 全record rowから同じ列をspliceする。
- `adv.depth_tracks`から同じ列を削除する。
- 既存Depth slotのindexをdiagram-globalに一度だけreindexする。
- 削除対象を参照するmanual slotはdisableし、理由をUIへ表示する。

#### Count and validation

- `uploadedDepthFileCount(row)`の最大値をlogical track countとして使わない。
- active logical indexは全record rowの同一列をunionして求める。
- slot selectorはglobalに存在するseriesだけを候補にする。
- partial coverageはwarningまたは`2/3 records`表示とし、generationを禁止しない。
- globalにsourceが一つもないseriesだけをinvalidとする。

#### Argument generation

- `--depth_track`ごとにrecord数と同じ数のpath/empty placeholderを保持する。
- globalに空の列を黙って圧縮しない。
- metadata argumentの順序をlogical track columnと一致させる。

combined multi-reference TSVを複数recordへ使う場合の`Apply to all records`操作は有用だが、
core bug fixの必須条件にはしない。

### Phase 4: Session normalization

serialized shapeは現行と同じため、identity修正だけではsession versionを上げない。

load時にlogical widthを次の最大値から求める。

```text
max(
    maximum saved record-row length,
    saved depth_tracks metadata length,
    maximum referenced Depth slot track_index + 1,
)
```

各record rowをこの幅まで`null`でpaddingする。保存済みの列順は変更しない。

Canonical requestはすでにDepth file/table matrix内の`null`を保持できるため、schemaを変更しない。
settings JSON、canonical session、gallery sessionについて、save/load後に次が同一であることを検証する。

- recordごとの欠損cell
- Depth metadata order
- slot `track_index`
- slot order
- axis index

session全体をpreflight validation/migrationした後にだけlive stateへcommitする。canonical data、slot geometry、
referenceが不正な場合はimportを失敗させ、現在のstateを変更しない。

過去にrecord-local spliceで既に失われたidentityは完全には復元できない。migrationは保存された位置を
そのまま論理列と解釈し、推測によるsource移動は行わない。

### Phase 5: Linear occupancy planner integration

Phase 0から4を、`docs/LINEAR_TRACK_OCCUPANCY_LAYOUT_IMPLEMENTATION_PLAN.md`の前提contractとする。

既存計画のPhase 2でDefault/Customを同じper-record plannerへ送る際、Depth slotへ次を渡す。

```text
logical track identity
data_available
paint footprint
reserve footprint
```

欠損cellでは`data_available=False`、paintなし、reserveありとする。これによりDefault/Customの
欠損semanticsと縦配置ownerを一致させる。

Phase 5ではSVG Y座標が変わる可能性があるため、Phase 1から4とは別のreview単位にする。
identity修正とgeometry差分を同じcommitへ混ぜない。

### Phase 6: Cleanup, documentation, and deployment

parity確認後に次を行う。

- record-local list位置をidentityとして使う互換helperを削除する。
- Custom専用Depth range checkを削除する。
- Web helpで`track_index`とvertical orderの違いを明記する。
- Depth tutorialとLinear layout tutorialへpartial coverageの意味を追記する。
- CLI Referenceの`--depth_track` placeholder例を検証する。
- packaging/browser検証用のbrowser wheelを再生成する。
- deploy時だけcache-bust tokenを更新する。

## 6. ファイル別の変更計画

### Core Depth model

- `gbdraw/analysis/depth_tracks.py`
  - `DepthTrackSpec.track_index`と`DepthTrackData.track_index`を追加する。
  - normalization、count、height、shared-axis、representative helperをlogical index基準へ変更する。
- `gbdraw/api/diagram.py`
  - global Depth countとslot validationをlogical index基準へ変更する。
  - Circular/Linearのprecomputed Depth dataへindexを伝搬する。
- `gbdraw/interface.py`
  - public `DepthTrackOptions`から作るmatrix列順がlogical indexになることをtestで固定する。
- `gbdraw/cli_utils/common.py`
  - empty placeholderを保持する現行contractを維持し、end-to-end testを追加する。

### Linear assembly and layout

- `gbdraw/diagrams/linear/assemble.py`
  - record-local `depth_by_index`を作る。
  - default/custom両経路をlogical lookupへ切り替える。
  - legend representativeを全rowから集約し、custom時はenabled slotを凡例の正本とする。
- `gbdraw/diagrams/linear/track_slots.py`
  - layout offsetはlogical indexで保持する。
  - Phase 5で`data_available`とreserve policyをper-record plannerへ渡す。
- `gbdraw/diagrams/linear/builders.py`
  - Depth placementへlogical indexを渡す。
- `gbdraw/diagrams/linear/positioning.py`
  - Phase 5までは互換配置を維持し、その後にresolved planへcutoverする。

### Circular consumers

- `gbdraw/diagrams/circular/assemble.py`
- `gbdraw/api/diagram.py`

共有Depth modelの明示indexへ追従する。Circular radial layoutの配置規則は変更しない。

### Web and session

- `gbdraw/web/js/app/depth-track-state.js`
  - clear-cellとremove-seriesを別操作にする。
  - logical column単位のreindex helperを追加する。
- `gbdraw/web/js/app/app-setup.js`
  - Linear add/removeを全record共通のseries operationへ変更する。
- `gbdraw/web/js/app/linear-track-slots.js`
  - active logical indicesとslot selectorを同期する。
- `gbdraw/web/js/app/run-analysis.js`
  - sparse matrix placeholderとmetadata順を保持する。
- `gbdraw/web/js/services/config.js`
  - load時にrecord rowsをlogical widthへpaddingする。
- `gbdraw/web/js/services/session-request.js`
  - `null` matrixのround-tripをtestで固定する。
- `gbdraw/session_io.py`
- `gbdraw/session_request_codec.py`
  - schema変更なしでsparse matrixが失われないことを検証する。

### Tests

- `tests/test_depth_track.py`
- `tests/test_linear_track_slots.py`
- `tests/test_session_request_codec.py`
- `tests/test_session_io.py`
- `tests/test_web_packaging.py`
- `tests/web/depth-track-state.test.mjs`
- `tests/web/depth-track-session.playwright.spec.js`

Phase 5のper-record vertical planner testは既存計画どおり`tests/test_linear_vertical_layout.py`へ置く。

## 7. Test plan

### 7.1 Pure normalization tests

最低限、次をassertする。

| Raw matrix | Normalized record 0 | Normalized record 1 | Global count |
|---|---|---|---:|
| `[[A], [None]]` | track 0 | empty | 1 |
| `[[None], [A]]` | empty | track 0 | 1 |
| `[[A, None], [None, B]]` | track 0 | track 1 | 2 |
| `[[A, B], [A, None]]` | tracks 0, 1 | track 0 | 2 |
| `[[None], [None]]` | error | error | N/A |

diagonal caseでは、record 1のBがtrack 0として扱われないことを必ず直接assertする。

### 7.2 Aggregation tests

- heightがlogical indexごとに選択される。
- colorとlabelが別系列へ移らない。
- shared-axis最大値が同じlogical indexのrecord間だけで共有される。
- first recordにtrackがなくてもlegend entryが作られる。
- 欠損recordの値をshared maxへ含めない。

### 7.3 Python integration and SVG tests

- partial track 0 + Custom Depth above same-side Features。
- partial track 0 + Custom Depth below Features。
- diagonal two-trackを異なるcolor/heightで描画する。
- 対象record/trackのSVG groupだけが存在する。
- 欠損recordにDepth axis/pathが存在しない。
- global slot reserveによりrecord stackが崩れない。
- slot order変更後も同じDepth dataが同じrecordへ描かれる。
- dense inputの既存group ID、axis、legend、geometryが不変である。

### 7.4 CLI end-to-end tests

mock forwardingだけでなく、一時Depth TSVと小さいGenBank/GFF3+FASTAを使って最後までrenderする。

```text
--depth_track record0.tsv ''
--linear_track_slot depth:depth@track_index=0
--linear_track_slot features:features
--linear_track_axis_index 2
```

さらにtwo-track diagonal placeholderを検証する。

### 7.5 Web unit and browser tests

- 一つのrecordのfileをclearしても他recordのcolumn indexが変わらない。
- series削除は全record、metadata、slotを同時にreindexする。
- sparse argsがpath/empty placeholderを保持する。
- 2 recordsの片方だけDepth upload、Custom Slots有効、DepthをFeatureより上へ移動、Generate成功。
- error panelが表示されない。
- DepthがあるrecordだけSVG groupを持つ。
- diagonal two-trackのcolor/labelが交差しない。
- Save、Load、Generate後もmatrix、slot order、axis indexが同じである。

### 7.6 Circular regression tests

- single-record multi-depth。
- circular multi-record partial Depth。
- custom circular slotによるtrack selection。
- shared-axisのlogical index分離。
- radial geometryと既存reference outputが不変である。

### 7.7 Reference output policy

Phase 1から4ではdense出力のreference SVG差分を0とする。`tests/reference_outputs/`を更新しない。

Phase 5で意図したgeometry変更が発生した場合だけ、既存occupancy計画に従って差分をreviewする。
差分はY transform、canvas height、comparison pathへ限定し、X座標、bp mapping、feature identity、colorを
変更しない。

## 8. Rollout order

次のreview単位に分ける。

1. Characterization tests。
2. Python logical identityとaggregation修正。
3. Linear/Circular renderer binding修正。
4. Web global-series operationとsession normalization。
5. Browser verificationと利用者向け文書更新。
6. Linear occupancy plannerへの統合。
7. Legacy positional ownerの削除。

長期feature flagは追加しない。Phase 2完了時点でsparse bindingを新しい唯一のcontractとする。
Phase 5のgeometry cutoverはshadow measurementとfocused testsを通してから行う。

## 9. 主なriskとmitigation

| Risk | Mitigation |
|---|---|
| 例外だけ止まり、別系列を描く | `track_index`をSpec/Dataへ保存し、diagonal sparse testで直接検出する |
| shared-axis、height、legendだけ旧list位置を使い続ける | aggregation helperをlogical index基準へまとめ、系列ごとのcolor/maxをtestする |
| Circularへ共有model変更が波及する | identityだけ同時修正し、radial geometryは非目標としてreference差分0を要求する |
| Webのrecord-local削除が他recordをずらす | clear-cellとremove-seriesを分離し、series削除をglobal atomic operationにする |
| 古いsessionの短いrowを誤解釈する | 保存順を維持し、最大logical widthまで`null` paddingする |
| 欠損recordだけstackが縮む | paintなし、reserveありをplanner contractにする |
| identity修正とlayout差分の原因が混ざる | Phase 1から4ではSVG parity、Phase 5を別review単位にする |
| record数xtrack数の性能低下 | recordごとにindex mapを一度作り、`O(records * tracks)`を維持する |

## 10. Non-goals

- 全recordへのDepth TSV提供を必須にしない。
- 欠損Depthを0 coverageとして描かない。
- Depth file内容からrecordへの割当を曖昧に自動推定しない。
- 初回修正でUUIDまたは任意文字列のpublic track IDを追加しない。
- 初回修正でsession/canonical request schemaを上げない。
- Circular layoutをLinear occupancy plannerへ統合しない。
- 一つのDepth bugのために汎用track plugin frameworkを導入しない。
- identity修正と全Linear rendererのgeometry cutoverを一つのcommitへまとめない。

## 11. Acceptance criteria

実装は次をすべて満たしたときに完了とする。

1. screenshot相当のpartial Depth + Depth above Featuresが成功する。
2. 同じ入力でDepth below Featuresも成功する。
3. slot order変更はY位置だけを変え、Depth source bindingを変えない。
4. diagonal sparse matrixで各recordが正しいlogical trackだけを描く。
5. per-record欠損はskipされ、globalに存在しないindexだけが早期errorになる。
6. 欠損recordにDepth group、axis、tick、pathが生成されない。
7. height、color、label、legend、shared-axisがlogical track間で混線しない。
8. Default/Customが同じDepth binding contractを使う。
9. 欠損slotのreserve policyがDefault/Customで一致する。
10. dense/shared/legacy入力の既存SVGがidentity phaseでは不変である。
11. Circularの既存Depth geometryが不変である。
12. settingsとcanonical sessionがsparse matrixを無損失round-tripする。
13. Webのclear-cellとremove-seriesが別操作として動作する。
14. focused Python、JS、browser、session、reference comparison testが通る。
15. Depth bindingとaggregationが`O(records * tracks)`を維持する。
16. Linear occupancy計画の完了後、Default/Customが同じper-record plannerを使う。

## 12. 完了後の文書更新

実装完了後に次を更新する。

- `docs/TUTORIALS/6_Depth_Quantitative_Tracks.md`
- `docs/TUTORIALS/7_Linear_Layout.md`
- `docs/CLI_Reference.md`
- `docs/FAQ.md`
- release notes
- Web help text

この計画書は実装中のsource of truthとし、phase完了時に状態と未完了項目を更新する。
