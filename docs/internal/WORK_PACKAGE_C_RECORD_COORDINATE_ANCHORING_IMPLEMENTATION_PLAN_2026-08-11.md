# Work package C: レコード単位の環状座標アンカー実装計画

- Date: 2026-08-11
- Status: ready for implementation
- Baseline: `docs_renovation` / `a2c69de112ae`
- Source: `docs/internal/gbdraw_v0.14.0_codex_roadmap.md` の Work package C
- Scope: Circular/Linear renderer、typed request、package-root Python API、CLI、Web、session、interactive SVG、documentation

この文書は、roadmap の「Circular display orientation」を、レコード単位の環状座標アンカーとして実装するための完全な計画である。元の Work package C が Circular の12時位置だけを対象としている箇所、`--circular_top_coordinate`、`top_coordinate`、既定値1を前提とする箇所は、本計画が置き換える。

文書作成時点では production code、test、fixture、Gallery artifact を変更していない。

## 1. 結論

一つの mode-neutral な `RecordDisplayTransform` を作り、Circular と Linear の全座標依存 consumer に渡す。

- Circular mode は指定した source coordinate を12時位置へ置く。
- Linear mode は同じ source coordinate を左端へ置き、レコード末端から先頭へ折り返して表示する。
- 入力 `SeqRecord`、配列、feature location、annotation table、depth table、comparison table は回転・再採番しない。
- tick、popup、download、interactive metadata は source coordinate を表示する。
- topology と表示起点はファイル単位でなく、解決後の個々の biological record が所有する。
- default は追加シフトなしであり、既存 Circular/Linear reference SVG を変えない。

`SeqRecord` を物理回転する実装、SVG group 全体を回転する実装、Circular/Linear に別々の座標式を追加する実装は採用しない。

## 2. 目標と対象範囲

### 2.1 目標

1. 一つの GenBank/GBFF 内に circular と linear の record が混在しても、各 record の topology と表示起点を個別に設定できる。
2. 完全な circular record を Circular と Linear のどちらでも同じ source coordinate にアンカーできる。
3. feature、label、tick、GC/skew、depth、annotation、comparison の全 geometry が同じ transform に従う。
4. source-coordinate identity、feature identity、match identity、session replay を維持する。
5. Web、CLI、typed API、package-root Python API、saved session で同じ意味を表現できる。
6. 起点変更は presentation-only change とし、LOSAT や source解析を再実行しない。

### 2.2 対象

- GenBank topology の検出と per-record override
- 1-based source coordinate による display start
- Circular single、grid、batch
- Linear single-record row、multi-record row、比較付きレイアウト
- reverse complement との合成
- ordinary、multipart、source-origin-spanning、display-cut-spanning feature
- horizontal/radial/embedded/external label と leader line
- ruler、tick、record definition の coordinate表示
- dinucleotide content/skew、depth
- Circular comparison ring、Linear nucleotide/protein comparison、Similarity group、Collinear block
- coordinate/feature annotation と automatic feature underlay
- interactive SVG、feature/match popup、sequence action、stable semantic hook
- canonical request schema、session writer/reader、Web editable draft、undo/redo
- CLI records table と単一 record convenience flags
- typed API と package-root API
- user documentation、release notes、必要な screenshot/capture contract

### 2.3 対象外

- source GenBank/GFF/FASTA の書換えや export
- biological sequence の再採番
- SnapGene `.dna` file の読書き
- circular region crop、origin をまたぐ crop、部分配列を circular record として扱う機能
- gapped alignment の塩基単位再構成、CIGAR/BTOP/aligned-sequence input の追加
- topology 推定アルゴリズム。GenBank metadata 以外は unknown とする
- Circular/Linear mode の自動切替
- plasmid restriction analysis、ORF prediction、assembly correction
- Work package E の自動 preview 再描画そのもの

## 3. 実装前に固定する製品判断

### D1. これは表示座標 transform であり、sequence editing ではない

`display start` は最初に表示する source base を選ぶ。source record の sequence、feature location、qualifier、record ID、comparison row、depth position は変更しない。

SnapGene の「起点を設定/回転」に近い操作感を提供するが、gbdraw は編集済み sequence file を作らない。source coordinate を1へ再採番せず、popup と tick は元の座標を示す。

### D2. diagram mode、record topology、display start を分離する

| 概念 | 値 | 責務 |
| --- | --- | --- |
| diagram mode | `circular` / `linear` | canvas と最終 geometry の形 |
| detected topology | `circular` / `linear` / `unknown` | source metadata の観測値 |
| topology override | `None` / `True` / `False` | user が record 単位で指定した意図 |
| effective circular | boolean | modular wrap と非default display start の許可 |
| display start | `None` または 1-based source base | 表示位置0に置く明示アンカー |

Circular diagram に含まれること自体は record を circular にしない。`effective circular=False` の record も従来どおり Circular diagram へ描けるが、表示起点の変更はできない。

### D3. default は `None = 追加シフトなし`

`start_coordinate` の public/default 値は `None` とする。これは現在の selection、crop、reverse-complement 処理後の先頭をそのまま使う identity path である。

明示値 `K` は source base `K` を最初に表示することを意味する。forward の完全 record では `K=1` が identity だが、reverse-complement 済み record では identity の先頭 source base は通常 `L` である。したがって `1` を無条件の default にしない。

Web の空欄は未設定を表す。placeholder または補助表示で `Current start: 1`、reverse complement 時は `Current start: L` のように現在値を示す。Reset は数値1を代入せず、override を `None` へ戻す。

### D4. topology は検出値と override を分ける

- GenBank/GBFF の `record.annotations["topology"] == "circular"` は detected circular。
- `linear` は detected linear。
- absent、空文字、その他の値、GFF3/FASTA、FASTA は unknown。
- user override があれば detected 値より優先する。
- override がなく unknown の record は、起点変更の許可判定では non-circular とする。
- Circular mode だから unknown を circular と推定する処理は入れない。

Web は単純な checkbox を表示してよいが、内部 state は `detectedTopology` と nullable `topologyOverride` を別々に保持する。override 後には `Reset to detected topology` を提供する。

### D5. setting owner は record である

`RecordDisplayOptions` は `RecordPresentation` の field ではなく、`RecordInput` の peer とする。topology と modular transform は label/placement ではないためである。

`RecordCardinality.ALL` の規則は次とする。

- default `is_circular=None, start_coordinate=None` は各展開 record で個別に topology を検出する。
- 明示した同じ override/start を全 record に適用する uniform fan-out は許可し、全 record で検証する。
- record ごとに異なる値を指定する場合は、selector 付き `EXACTLY_ONE` の `RecordInput` を複数作る。
- Web と CLI records table は、一つの source resource を共有する複数の exact-one record に materialize する。

record ID は重複し得るため、Web draft key は record ID だけを使わない。source UID と selector/local record index の組を使う。

### D6. v0.14.0 では crop と modular display start を併用しない

部分配列は complete circular record ではない。`RecordInput.region` または `RecordCollectionOptions.regions` が有効な record に明示 `start_coordinate` がある場合は `ValidationError` とする。

source record の circular metadata は消さなくてよいが、UI は crop 中の display-start control を disabled にし、理由を示す。将来 circular slice を導入する場合は別の coordinate contract を作る。

reverse complement は対応する。処理順は selection、既存 reverse complement、region validation、modular display transform とする。明示アンカーは最後まで source coordinate である。

### D7. base coordinate と boundary coordinate を別 API にする

一つの `position` 関数へ異なる座標規約を混ぜない。記号は次のとおりとする。

- `L`: complete source record length
- `a`: 1-based source base の明示 display start
- `s`: display orientation。forward は `+1`、reverse complement は `-1`
- `b`: 1-based source base、`1 <= b <= L`
- `u`: 0-based source boundary。source base `b` の boundary は `b-1` と `b`

明示アンカーがあるとき、source base から 0-based display base index への写像は次である。

```text
display_base_index(b) = mod(s * (b - a), L)
```

display offset 0 に対応する source boundary は orientation で異なる。

```text
anchor_boundary = a - 1   when s = +1
anchor_boundary = a       when s = -1

display_boundary_offset(u) = mod(s * (u - anchor_boundary), L)
```

scalar boundary の modulo 結果は `[0, L)` とする。fragment projector が返す geometry endpoint だけは `L` を使用できる。0/L boundary、full-length span、zero-length span を endpoint の一致だけで判定しない。

明示アンカーが `None` の場合は既存 record-local geometry をそのまま返す fast path とし、浮動小数点計算や fragment allocation を増やさない。

### D8. 入力座標ごとに adapter を一度だけ持つ

| Input | 入力規約 | transform へ渡す規約 | 表示 metadata |
| --- | --- | --- | --- |
| Biopython feature | 0-based half-open boundary | oriented interval/parts | source location parts |
| annotation coordinate | 1-based inclusive | half-open source interval | table の source values |
| BLAST q/s endpoints | 1-based inclusive、方向付き | oriented half-open span + original endpoints | original outfmt values |
| depth | 1-based base position | source base index/sample | original position |
| GC/skew window | record-local sample/window | source/display series point | source/window values |
| Circular/Linear tick | source label value | source base/boundary geometry | source coordinate label |
| feature shortcut | source location parts | 1-based source base | selected feature identity |

各 reader/drawer が独自に `-1`、`+1`、`360 * position / length` を追加しない。adapter の unit test で 1 bp feature、base 1、base L、boundary 0/L、reverse complement を固定する。

### D9. interval と series は seam で分割する

`RecordDisplayTransform` は scalar mapping だけでなく、最低限次を返す。

- `project_interval()` → 1個以上の `DisplayFragment`
- `project_parts()` → biological part と artificial seam を保持した fragment 列
- `project_series()` → display順へ並べた連続 series segment
- `alignment_cut_breakpoints()` → comparison の共通 parameter 上の cut

`DisplayFragment` は source/local span、display start/end、orientation、biological start/end の包含、artificial start/end、part index を持つ。

- feature の人工 seam には arrowhead、bracket cap、line cap を付けない。
- biological terminal を含む fragment だけが方向性 terminal を描く。
- Linear series は seam で path を閉じ、端点値を補間して左右端へ複製する。全幅を横切る diagonal を作らない。
- Circular path は mapped first point から開始し、source順を理由に不要な chord を作らない。

### D10. feature shortcut は 5′ end と covered-base midpoint

UI label は曖昧な `feature start` ではなく次とする。

- `Use selected feature 5′ end`
- `Use selected feature midpoint`

5′ end は biological order の最初の source base とする。midpoint は min/max envelope の平均ではなく、ordered location parts の covered bases をたどる。長さ `N` の feature では 5′端から `floor((N - 1) / 2)` base 進んだ塩基を選ぶ。これにより偶数長の丸めが strand に対して決定的になる。

unknown/mixed strand、empty location、record length 不明、selected record 不一致の場合は shortcut を disabled にし、推測しない。shortcut は draft を変更するだけで、committed preview は次の成功した Generate まで変えない。

### D11. cut-crossing HSP は既存 endpoint geometry の線形補間として分割する

現在の comparison input は q/s endpoints、alignment length、gap count を持つが、CIGAR や aligned sequence を持たない。したがって gapped HSP の cut に対応する相手側の厳密な塩基位置は復元できない。

v0.14.0 では、現在の endpoint-based ribbon と同じ幾何学的仮定を仕様化する。HSP を `t in [0,1]` の線形 parameter とし、query/subject の cut が発生する全 `t` を合併し、同じ `t` 区間で双方を同時に分割する。query と subject を独立分割して直積を作らない。

この分割は表示上の近似であり、塩基単位 alignment の主張ではない。将来 CIGAR/BTOP を受け入れる場合だけ exact mapping へ置き換える。

一つの source HSP から複数 SVG path が生じても、共通 `data-gbdraw-match-id` と catalog entry を使う。各 path は unique DOM ID と `data-gbdraw-match-fragment` を持つ。interactive collector は同一 match ID の metadata 一致を検証して全 fragment を一つの popup/sequence action へ束ねる。fragment ごとに別 match identity を作らない。

### D12. topology-aware annotation validation は record 解決後に行う

`wraps_origin=true` を diagram mode だけで拒否しない。parse は構文を検証し、record 解決後に effective circular と crop の有無を使って意味検証する。

- complete effective-circular record は Linear mode でも source-origin-spanning annotation を許可する。
- effective non-circular record または cropped record は明示エラーにする。
- ordinary source interval が display seam をまたぐ場合の内部 fragment 化は、user が `wraps_origin=true` を指定したこととは別概念である。

### D13. display start は analysis/cache identity に入れない

display transform は sequence content、selected protein、LOSAT arguments、source comparison rowsを変えない。LOSATN、TLOSATX、LOSATP、derived grouping/cache identity に display start を含めない。

Generate は既存 comparison evidence を再利用し、renderer だけを実行する。失敗、cancel、stale result は last successful preview、feature catalog、exportを置き換えない。

## 4. 現状の根拠

実装開始時の Phase 0 で行番号と baseline を取り直す。次の owner と制約は計画時点の `a2c69de112ae` で確認済みである。

| 現状 | Owner | 影響 |
| --- | --- | --- |
| per-record 設定は label、subtitle、reverse complement、grid placement だけ | `gbdraw/api/requests.py` の `RecordPresentation` / `RecordInput` | topology/display start の owner がない |
| 全 mode が同じ record resolver を通る | `gbdraw/api/record_planning.py::resolve_record_inputs()` | topology/anchor validation を一箇所へ置ける |
| planner materialization と Circular batch item が `RecordInput` を明示コピーする | `gbdraw/api/request_render.py` | 新 field を追加してもコピー点を更新しないと session/batch で消える |
| source/display mapping は `gbdraw_coord_base + gbdraw_coord_step * i` の affine model | `gbdraw/core/record_metadata.py`、`gbdraw/io/record_select.py`、`gbdraw/io/regions.py` | modular shift は現 modelだけでは表せない |
| Circular planner は topology linear を警告するだけ | `gbdraw/api/request_render.py::_warn_circular_topologies()` | diagram mode と biological topology が分離されていない |
| GenBank loader は record ごとの topology annotation を保持する | `gbdraw/io/genome.py`、`tests/test_genome_loading.py` | mixed-topology source は loader変更なしで扱える |
| GFF3/FASTA は reliable topology を持たない | genome/Web discovery path | unknown と explicit override が必要 |
| Linear Web の一つの `linearSeqs[]` entry は source card で、Automatic 時は複数 canonical record へ展開する | `gbdraw/web/js/state.js`、`app/record-discovery.js`、`services/config.js` | checkbox を card に直接置けない |
| Circular Web discovery は selector、ID、length を保持するが topology を捨てる | `gbdraw/web/js/app/record-discovery.js`、`app/python-helpers.js` | discovery result の共通拡張が必要 |
| records table は一行一 displayed record | `gbdraw/io/cli_tables.py` | mixed file の正式な per-record CLI surface に使える |
| package-root `draw_circular()` / `draw_linear()` は `SeqRecord` から default `RecordInput` を合成する | `gbdraw/interface.py` | beginner API に明示 per-record adapter が必要 |
| canonical request schema 5 は record/presentation fields を exact decode する | `gbdraw/session_request_codec.py` | additive fieldでも schema bump が必要 |
| session 40 と schema 5 は `main` first-parent に存在する | `main` の `gbdraw/session_io.py` と `session_request_codec.py` | 40/5 は positive compatibility input として維持する |
| decoder に `schema == CANONICAL_REQUEST_SCHEMA`、Web/session reader に `version == CURRENT` がある | Python/Web session modules | 定数だけ上げると 5/40 を誤って legacy 扱いする |
| Circular geometry に direct `360 * position / length` が多数ある | `gbdraw/svg/circular_features.py`、`circular_ticks.py`、`labels/circular.py` など | SVG group rotation では全 consumerを直せない |
| Linear ruler/definition は affine base/step を直接読む | `render/groups/linear/seq_record.py`、`definition.py`、`diagrams/linear/assemble.py` | source label と display x を分離する必要がある |
| Linear annotation reader/resolver は `wraps_origin` を mode で拒否する | `gbdraw/annotations/io.py`、`resolve.py` | record topology 後の validation へ移す必要がある |
| HSP は endpoints から一つの ribbon/path を作る | `render/groups/linear/pairwise_match.py` | cut-crossing は共通 t-space projection が必要 |
| interactive match collector は match ID の重複を拒否する | `gbdraw/web_support/feature_catalog.py`、`render/interactive_svg.py` | fragment group supportを comparison と同時に入れる |
| feature identity は source location parts を既に保持する | `gbdraw/web_support/feature_metadata.py`、`features/ids.py` | display geometry と biological identity を分けられる |

## 5. Public contract

### 5.1 Typed request

`gbdraw/api/requests.py` に frozen dataclass を追加する。

```python
@dataclass(frozen=True)
class RecordDisplayOptions:
    is_circular: bool | None = None
    start_coordinate: int | None = None
```

意味は次のとおり。

- `is_circular=None`: source metadata から検出する。
- `True` / `False`: user override。
- `start_coordinate=None`: 追加 shift なし。
- positive integer: 1-based source base を最初に表示する。

validation:

- `is_circular` は strict bool または `None`。整数 0/1 を bool として受けない。
- `start_coordinate` は bool でない正の整数または `None`。
- `is_circular=False` と明示 start の組は request construction 時に拒否する。
- upper bound、auto topology、region/collection transformとの競合は record loading 後に planner が検証する。
- `RecordInput.display` は既存 positional constructor を壊さないよう、現在の全 field の末尾へ追加する。

`RecordDisplayOptions` は `gbdraw.api` から export する。`RecordPresentation` へ compatibility alias を置かない。

### 5.2 Resolved plan contract

planner は record ごとに immutable な resolved context を作る。

```python
@dataclass(frozen=True)
class ResolvedRecordDisplay:
    detected_topology: Literal["circular", "linear", "unknown"]
    is_circular: bool
    start_coordinate: int | None
    current_start_coordinate: int
    orientation_step: Literal[-1, 1]
```

`ResolvedRecordProvenance` と request plan の record-aligned tuple に保持する。render plan は同じ順序の `RecordDisplayTransform` tuple を所有する。

この resolved context と transform を `SeqRecord.annotations` へ書かない。source metadata と transient display state を混ぜず、batch projection、session encoding、interactive context が plan から同じ値を受け取る。

### 5.3 Canonical request schema 6

schema 6 の各 record は `display` object を持つ。

```json
{
  "recordKey": "record-1",
  "source": {"kind": "genbank", "resourceId": "..."},
  "selector": null,
  "region": null,
  "presentation": {
    "label": null,
    "subtitle": null,
    "reverseComplement": false,
    "gridRow": null,
    "gridColumn": null
  },
  "display": {
    "isCircular": null,
    "startCoordinate": null
  }
}
```

wire contract:

- `isCircular` は boolean または `null`。`null` は auto detection。
- `startCoordinate` は positive integer または `null`。
- schema 6 では object と2 fieldを必須とし、unknown field を拒否する。
- schema 5 decode は `is_circular=None, start_coordinate=None` を補う。
- schemas 1/2 の既存migrationを維持し、3/4は引き続き拒否する。
- supported schemas は `{1, 2, 5, 6}`。future 7 は拒否する。

`schema == CANONICAL_REQUEST_SCHEMA` を current-format 判定に使わない。少なくとも次の historical era predicate を導入する。

```python
GROUPING_SCHEMA = 5
RECORD_DISPLAY_SCHEMA = 6
```

schema 5/6 は同じ top-level grouping/output shape を使い、record display field の有無だけを version threshold で分ける。

`_materialized_record_inputs()`、`_is_materialized_exact_one_request()`、Circular batch `item_plans()`、Web record projection の全 explicit-copy point を更新する。

### 5.4 Session version 41

session 41 は canonical schema 6 と Web の per-record editable draft を保存する。

accepted sessions は `27–33, 39, 40, 41` とする。34–38 は development-only のまま拒否する。session 40 は current-authority era の positive input であり、41への bump を理由に legacy field migrationを再実行しない。

historical threshold を導入する。

```text
CURRENT_AUTHORITY_SESSION_MIN_VERSION = 40
RECORD_DISPLAY_SESSION_VERSION = 41
```

次を明示的に監査する。

- Python `== CURRENT_SESSION_VERSION` / `< CURRENT_SESSION_VERSION`
- Web `=== SESSION_VERSION` / `!== SESSION_VERSION` / `< SESSION_VERSION`
- feature catalog expansion
- `files` と `webFiles` authority validation
- Linear comparison draft migration
- CLI-only replay
- current writer validation message

v40/schema5 を load して Generate せず Save した場合も、session 41/schema6 を書く。schema5 record へ `display: {isCircular: null, startCoordinate: null}` を加える pure canonical promotion を行い、v41/schema5 を偶発的に生成しない。これは geometryを変えず、再ロード時に embedded source metadata から topology を検出できる。

current writer fixtureとは別に、`main` writerで作った v40/schema5 fixture を immutable positive control として追加する。render requestだけでなく、record order、editor feature catalog、layout、comparison authority、resource bindingを検証する。

### 5.5 Package-root Python API

`gbdraw.draw_circular()` と `gbdraw.draw_linear()` に同じ keyword-only parameter を追加する。

```python
record_displays: Sequence[RecordDisplayOptions] | None = None
```

- `None` は全 record default。
- sequence length は normalized record count と完全一致させる。
- string/bytes、wrong type、余剰/不足を `ValidationError` にする。
- adapter は各 in-memory `RecordInput.display` へ渡すだけで、topology/file validationを再実装しない。
- `RecordDisplayOptions` は package root からも export する。

signature、`gbdraw.__all__`、`tests/fixtures/public_contract.json`、Python reference と実行可能 exampleを意図的に更新する。

### 5.6 CLI

正式な per-record surface は両 mode 共通の `--records_table` とする。optional columns を追加する。

| Column | Values | Default |
| --- | --- | --- |
| `topology` | blank/`auto`、`circular`、`linear` | `auto` |
| `display_start` | blank または positive 1-based source coordinate | blank / no shift |

table parse は enum、positive integer、明示 `linear + display_start`、`region + display_start` を row/column付き error で検証する。upper bound と auto topology は record materialization 後に検証する。

単一 displayed record のため、両 subcommand に同じ convenience flags を追加する。

```text
--record_topology {auto,linear,circular}
--display_start_coordinate INT
```

- direct `--gbk` / `--gff --fasta` input が exactly one record に解決される場合だけ許可する。
- multi-record解決では明示エラーを出し、`--records_table` を案内する。
- records table との併用は拒否する。
- `--circular_top_coordinate` や mode-specific alias は追加しない。
- session replayに fresh flag を混ぜる既存禁止規則を維持する。

parser、namespace normalization、public CLI contract snapshot、generated help、command-line referenceを両 mode 同時に更新する。

### 5.7 Web canonical/draft contract

committed meaning と editable draft を分ける。

```text
renderRequest.records[].display
    = last successful render の意味

config.recordDisplayDrafts[]
    = 未生成を含む user editable state

webFiles
    = resource bytes と source binding のみ
```

draft entry は最低限次を持つ。

```json
{
  "scope": "circular",
  "sourceUid": "...",
  "selector": "#2",
  "recordId": "duplicate-or-public-id",
  "topologyOverride": null,
  "startCoordinate": null
}
```

`detectedTopology`、length、current default start は discovery result から導出し、persisted authority にしない。session load 中の表示に必要な一時 provenance として保持しても、source再発見結果を上書きしない。

`webFiles.bindings` に topology/start を複製しない。history snapshot と Reset Settings は draftを扱い、committed requestを変更しない。

## 6. Target architecture

### 6.1 Data flow

```text
source resource
  -> record discovery / typed source loading
  -> selector and cardinality resolution
  -> source topology detection + user override
  -> existing reverse-complement / affine source mapping
  -> per-record RecordDisplayTransform
  -> interval / series / alignment fragment projection
  -> Circular angle or Linear x geometry
  -> SVG parts with unchanged source identity metadata
```

fresh Web generation、session replay、CLI、package-root Python、typed API は全て typed planner 以降の同じ経路へ合流する。

### 6.2 Module ownership

| Owner | Responsibility |
| --- | --- |
| `gbdraw/api/requests.py` | raw per-record option と早期型検証 |
| `gbdraw/api/record_planning.py` | topology検出、override、length/region validation、provenance |
| `gbdraw/core/record_metadata.py` | existing source/local affine mapping。modular policyは持たない |
| new `gbdraw/layout/record_coordinates.py` | mode-neutral immutable transform、fragment、series、alignment breakpoint |
| `gbdraw/api/request_render.py` | plan-aligned transform tuple、materialization、batch projection |
| Circular/Linear assembly | record indexに対応する transform の中継 |
| drawers/groups/svg helpers | projection resultからgeometryを作る。topologyを再解決しない |
| `session_request_codec.py` | canonical schema 6 encode/decode |
| `session_io.py` / `api/session_compat.py` | session 41 と evidence-backed migration |
| `io/cli_tables.py` | CLI table shape とcell validation |
| Web `app/record-display-options.js` | discoveryとのreconcile、UI projection、feature shortcut |
| Web `services/session-request.js` | Web stateからcanonical requestへの唯一のprojection |
| Web `services/config.js` | draft/committed/session coordination |

transformを `features`、`labels`、`annotations`、各modeへ複製しない。`RecordDisplayTransform` は public API として exportせず、planner/render internal contract とする。

### 6.3 Transform API

実装時の正確な naming は型とtestを同時に更新する限り調整してよいが、責務を減らさない。

```python
@dataclass(frozen=True)
class RecordDisplayTransform:
    length: int
    source_base: int
    source_step: Literal[-1, 1]
    start_coordinate: int | None
    is_circular: bool

    def source_base_to_display_index(self, base: int) -> int: ...
    def source_boundary_to_display_offset(self, boundary: int) -> int: ...
    def display_index_to_source_base(self, index: int) -> int: ...
    def project_interval(self, interval: SourceInterval) -> tuple[DisplayFragment, ...]: ...
    def project_parts(self, parts: Sequence[SourceInterval]) -> tuple[DisplayFragment, ...]: ...
    def project_series(self, points: Sequence[SeriesPoint]) -> tuple[DisplaySeries, ...]: ...
    def alignment_cut_breakpoints(self, span: SourceInterval) -> tuple[float, ...]: ...
```

要件:

- constructor は length > 0、step ±1、start range、circular prerequisite を検証する。
- explicit start を source/local cutへ変換する処理は一度だけ行う。
- source↔display scalar mapping は round-trip propertyを持つ。
- interval length、orientation、biological terminalを失わない。
- full-length spanをzero-lengthに畳まない。
- identity transform は既存 raw coordinatesを返す fast pathを持つ。
- geometry scalarやSVG objectへ依存しない。

### 6.4 Renderer consumer matrix

| Consumer | Circular | Linear | Required behavior |
| --- | --- | --- | --- |
| feature/underlay | angle/path projection | seam fragment | source identity、terminal arrow、part ordering維持 |
| horizontal/radial label | display midpointとleader | continuous fragment選択 | label text/filterは不変 |
| tick/ruler | source label、mapped angle | source label、mapped x | 1/L jump、RCを正しく表示 |
| record definition | length/source range | wrapped coordinate summary | 誤解を招く単一 `start-end` を避ける |
| GC/content/skew | mapped closed path | reordered/broken series | seam chord/diagonalなし |
| depth | mapped samples | reordered/broken series | logical track index不変 |
| annotation | mapped span/parts | seam fragment | artificial capなし、topology-aware wrap |
| Circular comparison ring | mapped reference span | N/A | raw HSP metadata維持 |
| Linear pairwise/collinear | N/A | common-t fragments | query/subject対応、match identity維持 |
| Similarity group alignment | N/A | mapped member center | record placement前に適用 |
| interactive SVG | all partsを同じsource identityへbinding | 同左 | popup/downloadはsource coordinate |

## 7. 実装フェーズ

### Phase 0: baseline、fixture、座標oracleを固定する

Status: pending

#### 作業

1. branch、HEAD、dirty fileを記録し、production、test、documentation、generated artifactの既存差分を分ける。
2. `main` first-parent の schema 5/session 40 constants と代表 session を記録し、compatibility fixture の由来を `tests/fixtures/sessions/README.md` に残す。
3. 既存 reference SVG の status と checksumを記録する。通常testで `tests/reference_outputs/` を書かない。
4. coordinate-domain oracle を test dataとして作る。
   - length 10
   - forward/reverse
   - start `None`, 1, 5, 10
   - base 1/10
   - boundary 0/10
   - 1 bp、ordinary、cut-crossing、full-length interval
5. 一つの source に circular/linear/unknown recordを含む multi-record fixtureを用意する。duplicate record ID の別fixtureまたは同fixtureを用意する。
6. seam 前後に ordinary、multipart、source-origin-spanning、±strand featureを置く fixtureを用意する。
7. unsampled cutを持つ GC/depth series、query-only/subject-only/both-cut の plus/minus HSP fixtureを用意する。
8. focused baseline testを変更前に実行し、pre-existing failureと所要時間を記録する。

#### Baseline commands

```bash
git status --short --untracked-files=all
git diff --stat
git rev-parse --short=12 HEAD
python -m pytest tests/test_api_requests.py tests/test_record_planning.py tests/test_api_request_render.py -q
python -m pytest tests/test_session_request_codec.py tests/test_session_io.py tests/test_session_compat.py -q
python -m pytest tests/test_output_comparison.py::TestOutputComparison -v
TMPDIR=/tmp node --test tests/web/session-request.test.mjs tests/web/session-authority.test.mjs
```

#### 完了条件

- 数式の期待値を実装から独立したoracleとしてreviewできる。
- v40/schema5 positive controlがbranch writerで上書きされない。
- default referenceの変更前baselineが得られる。

### Phase 1: mode-neutral coordinate transform を実装する

Status: pending

#### 対象owner

- new `gbdraw/layout/record_coordinates.py`
- 必要な affine read helper のみ `gbdraw/core/record_metadata.py`
- new `tests/test_record_display_transform.py`

#### 作業

1. source base、source boundary、record-local、display indexの型/小dataclassを定義する。単なる `int` を異なるdomain間で無検証に渡すpublic helperを作らない。
2. `RecordDisplayTransform` と `DisplayFragment`、`DisplaySeries` を frozen contractとして実装する。
3. identity fast path、explicit forward、explicit reverseを実装する。
4. scalar inverse、interval/part projection、series seam interpolation、alignment cut breakpointを実装する。
5. full-length span、zero-length、boundary 0/L、modulo重複を別々に扱う。
6. invalid length、step、source base、start、non-circular startを typed errorへ変換できるようにする。
7. transform moduleから Circular/Linear/SVG/Biopython objectをimportしない。

#### Unit/property tests

- D7 の全oracle値
- `display_to_source(source_to_display(x)) == x`
- forward/reverseで全base permutationが重複しない
- projected fragment length合計がsource span lengthと一致する
- seamなし intervalは1 fragment
- seamあり intervalは2 fragment
- full-length intervalはlength Lのまま
- multipart order、terminal、artificial seam flags
- sampled/unsampled cutのseries endpoint補間
- t breakpointの0/1除外、sort、deduplicate
- identity callが入力interval/seriesを不要に分割しない

#### 完了条件

- Circular/Linear codeをimportせず全座標oracleが通る。
- rendererが独自modulo計算を追加せずに必要なfragmentを得られる。

### Phase 2: typed model、record resolution、plan threadingを実装する

Status: pending

Depends on: Phase 1

#### 対象owner

- `gbdraw/api/requests.py`
- `gbdraw/api/record_planning.py`
- `gbdraw/api/request_render.py`
- `gbdraw/api/__init__.py`
- `gbdraw/api/prepared.py` または既存prepared contract
- `tests/test_api_requests.py`
- `tests/test_record_planning.py`
- `tests/test_api_request_render.py`

#### 作業

1. `RecordDisplayOptions` を追加し、`RecordInput.display` を末尾へ追加する。
2. public export と `__all__` を更新する。
3. topology normalizerを一つ作り、raw annotationの大小文字/unknownを決定的に処理する。
4. selectionされた各 raw recordについて detected、override、effective circularを解決する。
5. existing affine metadataから current source start と orientation stepを取得し、transformを構築する。
6. start upper bound、circular prerequisite、record-level/collection-level region conflictをplannerで検証する。errorには input index、resolved record ID、record keyを含める。
7. `ALL` default inference、uniform fan-out、heterogeneous exact-oneの契約を実装する。
8. `ResolvedRecordProvenance`、Circular/Linear request plan、prepared inputへrecord-aligned display/transformを保持する。
9. `_coerce_resolved_collection()`、`_materialized_record_inputs()`、`_is_materialized_exact_one_request()`、Circular batch item projectionのexplicit copyを更新する。
10. Circular topology warningはdiagram mode warningとして残してよいが、effective topology resolverを再実装しない。explicit overrideを考慮した文面にする。
11. schema番号を埋め込んだdocstring/errorを `current canonical request` へ直す。

#### Tests

- strict model type/range/default/positional compatibility
- detected circular/linear/unknown と True/False override
- one mixed source + `ALL` + auto
- uniform override fan-outの全record validation
- selector付きexact-oneでheterogeneous設定
- reverse complement current start と explicit source anchor
- record region / collection region rejection
- resolved/materialized requestがdisplayを失わない
- Circular single/grid/batch、batch item、Linear planのrecord-transform alignment
- default planのdiagram output parity

#### 完了条件

- planの各 recordにちょうど一つのtransformが対応する。
- batch、resolve→plan、already-materialized short pathで値が消えない。

### Phase 3: schema 6、session 41、Web writerを原子的に更新する

Status: pending

Depends on: Phase 2

このPhaseはPython codecだけ、またはWeb constantsだけを先にlandしない。

#### Python canonical/session owner

- `gbdraw/session_request_codec.py`
- `gbdraw/session_io.py`
- `gbdraw/session.py`
- `gbdraw/api/session_compat.py`
- `gbdraw/cli_utils/session.py`
- related session fixtures

#### Web canonical/session owner

- `gbdraw/web/js/services/session-request.js`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/services/session-authority.js`
- `gbdraw/web/js/services/session-resources.js` verification-only unless binding changes are found

#### 作業

1. writer schemaを6、supported schemasを1/2/5/6へ更新する。
2. record `display` exact encode/decodeとschema5 default projectionを追加する。
3. grouping/top-level decodeを `schema >= 5` のformat-era predicateへ変える。
4. current sessionを41、supportedを27–33/39/40/41へ更新する。
5. current-authority behaviorを `version >= 40`、record-display draft availabilityを `version >= 41` で判定する。
6. `< CURRENT_SESSION_VERSION` によるv40誤migrationを全件分類し、historical thresholdへ置換する。
7. v40/schema5 load→save without Generateのpure promotionを実装する。
8. Python writer、Web writer、session authority inventoryのversion/schemaを同期する。
9. v40 main-era fixture、v41 current fixture、v39/v33 regression fixture、rejected 34–38を別々にtestする。
10. version/schemaをhard-codeするdocumentation/capture/packaging testを更新する。

#### Required session matrix

| Input | Expected read | Expected save/replay |
| --- | --- | --- |
| v41/schema6 explicit | current | same semantics |
| v41/schema6 auto/null | current | source metadataでresolve |
| v40/schema5 | positive current-authority input | v41/schema6、no shift |
| v39/schema5 | supported historical Web | direct supported migration |
| v33/schema1/2 | supported typed bridge | current requestへmigration |
| v27–30 | CLI replay only | existing boundary維持 |
| v34–38 | reject | reject |
| future v42/schema7 | reject | reject |

#### Tests

- schema6 exact fields、missing、unknown、wrong type、range
- schema5 single/grid/batch grouping positive control
- schema1/2 positive、3/4 negative
- Web request build/project round-trip
- v40のrecord order、catalog、editor、layout、comparison、resource authority
- load v40→Save before Generateがv41/schema6になる
- current sessionにlegacy `files` があればv40/v41とも拒否
- v40にold option/CLI/comparison migrationを再適用しない

#### 完了条件

- PythonとWeb writerが41/6で一致する。
- v40/schema5がlegacy誤判定されず、図とeditable authorityを維持する。
- current writerからv41/schema5が生成されない。

### Phase 4: Circular consumerをshared transformへ移行する

Status: pending

Depends on: Phase 1–3

#### Assembly/context

- `gbdraw/layout/circular.py`
- `gbdraw/diagrams/circular/assemble.py`
- `gbdraw/diagrams/circular/builders.py`
- Circular render context/configurator boundaries

各 record contextへtransformを一度だけ渡す。multi-record grid/batchも同じplan tupleからrecord indexで取得する。

#### Features、underlays、labels

- `gbdraw/render/drawers/circular/features.py`
- `gbdraw/svg/circular_features.py`
- `gbdraw/features/coordinates.py` / `objects.py` のdomain adapter
- `gbdraw/annotations/feature_underlays.py`
- `gbdraw/labels/circular_candidates.py`
- `gbdraw/labels/circular.py`
- `gbdraw/labels/circular_radial.py`
- `gbdraw/render/drawers/circular/labels.py`

source-domain part coalescingとfeature identityを先に確定し、geometryだけをprojectする。candidateに source midpoint と display midpointを分ける。horizontal/radial solver、leader anchor、embedded判定、collision obstacleがdisplay geometryを使う。

#### Ticks、definitions、numeric tracks

- `gbdraw/render/groups/circular/ticks.py`
- `gbdraw/svg/circular_ticks.py`
- `gbdraw/render/groups/circular/definition.py`
- `gbdraw/render/drawers/circular/gc_content.py`
- `gbdraw/render/drawers/circular/gc_skew.py`
- `gbdraw/render/drawers/circular/depth.py`
- `gbdraw/svg/circular_tracks.py`

tick valueとangleを分ける。display index 0を12時の `-90 degree` とする。series pathのstart/end、baseline closure、percent/deviation axisを監査する。value axis geometryは回転しない。

#### Annotations、comparison rings、interactive

- `gbdraw/annotations/resolve.py`
- `gbdraw/render/drawers/circular/annotations.py`
- `gbdraw/analysis/conservation.py`
- `gbdraw/render/drawers/circular/conservation.py`
- `gbdraw/svg/circular_conservation.py`
- interactive context/metadata builders

raw annotation/HSP metadataを維持し、display spanだけをprojectする。source-origin spanとdisplay seam spanを区別する。

#### Focused tests

- `tests/test_circular_feature_paths.py`
- `tests/test_feature_underlay_rendering.py`
- `tests/test_circular_ticks.py`
- `tests/test_circular_label_placement.py`
- `tests/test_circular_radial_labels.py`
- `tests/test_circular_annotation_tracks.py`
- `tests/test_circular_conservation.py`
- `tests/test_depth_track.py`
- Circular multi-canvas、SVG ID integrity、interactive metadata tests

#### 完了条件

- features、tracks、labels、ticks、annotations、comparison ringsが同じanchorを示す。
- radial/horizontal labelとleaderに衝突・逆向き・全円chordがない。
- default Circular referenceはbyte-identical。

### Phase 5: Linear record、feature、ruler、track、annotationを移行する

Status: pending

Depends on: Phase 1–3

#### Assembly/context

- `gbdraw/layout/linear.py`
- `gbdraw/diagrams/linear/assemble.py`
- `gbdraw/diagrams/linear/builders.py`
- `gbdraw/diagrams/linear/precalc.py`
- `gbdraw/diagrams/linear/positioning.py`

Linear global contextへ一つのtransformを置かず、record index→transform tupleをassemblyが所有する。row scale、normalize length、record gap、center alignmentはlengthを使い続け、display startで変えない。

#### Features、labels、underlays

- `gbdraw/render/drawers/linear/features.py`
- `gbdraw/svg/linear_features.py`
- `gbdraw/labels/linear.py`
- `gbdraw/render/drawers/linear/labels.py`
- `gbdraw/annotations/feature_underlays.py`

cut-crossing partをfragment化し、artificial seamにterminal arrow/capを付けない。同じ feature stable IDを維持し、fragment DOM IDだけを一意にする。embedded labelは最長continuous fragmentを選び、同長はdisplay start、part order、stable IDの順で決定する。

#### Ruler、definition、tracks

- `gbdraw/render/groups/linear/seq_record.py`
- `gbdraw/render/groups/linear/definition.py`
- `gbdraw/render/groups/linear/length_bar.py`
- `gbdraw/render/drawers/linear/gc_content.py`
- `gbdraw/render/drawers/linear/gc_skew.py`
- `gbdraw/render/drawers/linear/depth.py`
- `gbdraw/svg/linear_tracks.py`

explicit display startを持つrecordのon-axis rulerは、multi-record rowでもmapped source coordinate labelを使う。global bottom scaleは物理距離のままにする。definition coordinate lineはforwardなら `[K..L], [1..K-1]`、reverseなら `[K..1], [L..K+1]` のwrapped notationを使い、空区間は表示せず、`hide_length` の既存visibilityに従う。したがって `K=1`、`K=L` でも範囲外の `0` や `L+1` を文字列として出さない。

seriesはdisplay順へ並べ、seamで別pathにする。cutにsampleがない場合は既存interpolation policyに従う新しいendpointを作る。logical depth track identity、legend、shared y-axisは変えない。

#### Annotations

- `gbdraw/annotations/io.py`
- `gbdraw/annotations/resolve.py`
- `gbdraw/render/drawers/linear/annotations.py`
- `gbdraw/annotations/layout.py`

mode-only `wraps_origin` rejectionをparseから外し、resolved record topology/cropで検証する。bracket/line/band/highlightをfragment化し、人工端にcap/arrowを置かない。lane packingとlabelはprojected boundsを使う。

#### Focused tests

- `tests/test_linear_feature_paths.py`
- `tests/test_linear_label_placement.py`
- `tests/test_linear_selectors.py`
- `tests/test_linear_definition_alignment.py`
- `tests/test_linear_track_layout.py`
- `tests/test_linear_annotation_tracks.py`
- `tests/test_annotations.py`
- `tests/test_depth_track.py`
- multi-record row / duplicate ID / vertical layout tests

#### 完了条件

- `K..L,1..K-1` geometryに負幅、全幅diagonal、余分なarrow/capがない。
- source tick/popupとdisplay xが正しく対応する。
- default Linear referenceはbyte-identical。

### Phase 6: Linear comparison、orthogroup alignment、fragment interactionを移行する

Status: pending

Depends on: Phase 5

#### 対象owner

- `gbdraw/render/groups/linear/pairwise_match.py`
- `gbdraw/diagrams/linear/orthogroup_alignment.py`
- `gbdraw/diagrams/linear/assemble.py`
- `gbdraw/diagrams/linear/builders.py`
- `gbdraw/web_support/feature_catalog.py`
- `gbdraw/render/interactive_svg.py`
- match popup/standalone interactive modules

#### 作業

1. query/subject oriented endpointsを共通 `t` parameterへ正規化する。
2. 両 record transformのcut breakpointを合併し、各t区間からribbon/curve fragmentを作る。
3. plus/plus、plus/minus、minus/plus、minus/minusを同じparameter規則で扱う。
4. gapを持つHSPはendpoint linear interpolationであることをmetadata/docsに明記する。
5. fragmentは元matchのsource q/s endpoints、score、identity、raw evidence identityを保持する。
6. common match ID + fragment indexをSVGへ付ける。
7. catalog/interactive collectorをmatch ID→element listへ変更し、metadata不一致のduplicateは明示エラーにする。
8. popupとsequence downloadは一つのsource HSPとして動作し、fragmentだけのsequenceを誤exportしない。
9. Similarity-group/Collinear member centerをtransform後にxへ変換し、その後にrecord translation/alignmentを適用する。
10. presentation-only start変更がLOSAT/cache keyへ入らないことをcounter付きtestで確認する。

#### Tests

- query-only、subject-only、both-cut
- cut exactly endpoint、cut inside span、no cut
- plus/minus全組合せ
- gapなし exact linear case、gapあり documented approximation
- ribbon/curve、nucleotide/protein、pairwise/collinear
- one source match→multiple paths→one catalog entry/popup
- duplicate metadata mismatch rejection
- fragment DOM ID uniquenessとsemantic common ID
- source sequence download、reference/comparison/both actions
- orthogroup anchorがdisplay transform後に揃う
- display start変更でLOSAT executor/cache probe 0

#### 完了条件

- comparison linkがcanvasを誤横断せず、query/subject対応を失わない。
- interactive match identityがfragment数で増殖しない。

### Phase 7: Web record discovery、draft state、UI、feature shortcutを実装する

Status: pending

Depends on: Phase 2–6。schema/session projectionはPhase 3のcontractを使う。

#### Discovery/state owner

- `gbdraw/web/js/app/record-discovery.js`
- `gbdraw/web/js/app/python-helpers.js`
- `gbdraw/web/js/app/linear-record-selector.js`
- `gbdraw/web/js/app/annotations/record-catalog.js`
- new `gbdraw/web/js/app/record-display-options.js`
- `gbdraw/web/js/state.js`
- `gbdraw/web/js/app/app-setup.js`

#### 作業

1. JS fast GenBank parser と Python/Biopython discovery が selector、record ID、length、normalized topologyを同じshapeで返すようにする。
2. GFF3/FASTA discovery は topology unknown とする。
3. source UID + selector/local indexからmode-scoped logical record keyを作る。record IDをkeyにしない。
4. `record-display-options.js` が discovery reconcile、override、inactive start、range validation、feature shortcut計算、stale entry purgeを所有する。
5. `linearSeqs[]` はsource cardだけを所有し、per-record display stateを追加しない。
6. file replace時だけ古いsource UIDのdraftを削除する。async rediscovery、record order、selector切替では対応するdraftを維持する。
7. Reset Settingsはtopology override/startをdefaultへ戻す。undo/redo snapshotはdraftを含む。

#### UI owner

- `gbdraw/web/index.html`
- focused CSS/template
- `gbdraw/web/js/app/app-setup.js`

Circular uploader直下に、single/grid/batchの別なく discovery済み record を表示する `Records` sectionを置く。Multi-Record Canvasがoffでも表示する。

Linearは各source cardの `Record options` 内にrecord rowsを置く。

- explicit selector: 選択recordの一行
- Automaticかつmulti-record: 発見recordごとの入れ子行
- source card共通checkboxは作らない

各行のvisible contract:

```text
<record ID>  <selector>  <length> bp  <Detected: circular|linear|unknown>
[x] Circular record                         [Reset to detected]
Display start [          ]  Current start: <source coordinate>
[Use selected feature 5′ end] [Use selected feature midpoint] [Reset start]
```

mode-specific label:

- Circular: `Coordinate at 12 o'clock`
- Linear: `Coordinate at left edge`

同じstate fieldへprojectし、modeごとに別名fieldを作らない。

#### UI規則

- checkboxはresolved effective valueを表示する。
- detected unknownはunchecked。ユーザーがcheckするとexplicit circular overrideになる。
- user toggle後はReset actionを表示し、autoへ戻せる。
- non-circular、crop、length不明ではstart inputとshortcutをdisabledにし、理由を行内表示する。
- start inputの空欄は追加shiftなし。`1 <= K <= length` を行内validateする。
- checkboxをoffにしても以前のstartはinactive draftとして保持してよいが、canonical requestには `startCoordinate:null` を出す。再checkでdraftを復元する。
- duplicate record IDにはselectorを必ず併記する。
- accessible nameにrecord label/selectorを含める。
- errorはGenerateを止めるが、last successful previewを消さない。

#### Feature shortcut

- last committed feature catalogのselectionを使う。
- exactly one selected featureで、record keyが対象rowと一致するときだけ有効。
- `location_parts` とstrandからD10の5′ end/midpointを計算する。
- stale committed catalog、different source UID、multiple selection、unknown strandではdisabled。
- shortcut actionはdraft/historyだけを更新し、LOSAT、file read、preview mutationを開始しない。

#### Canonical/session projection

- `gbdraw/web/js/services/session-request.js`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/services/reset.js`
- `gbdraw/web/js/services/history-snapshot.js`
- `gbdraw/web/js/services/session-authority.js`

Automatic Linear expansion時、一つのsource resourceから複数exact-one canonical recordsを作り、logical keyで各display optionを対応させる。resource bytesは一つのままにする。

session loadは committed canonical displayを先に復元し、同じlogical recordにだけeditable draftをoverlayする。未生成draftをcommitted previewへ混ぜない。`webFiles` は変更しない。

#### JS unit tests

- new `tests/web/record-display-options.test.mjs`
- `tests/web/record-selector.test.mjs`
- `tests/web/session-request.test.mjs`
- `tests/web/session-authority.test.mjs`
- `tests/web/session-draft-authority.test.mjs`
- `tests/web/session-resources.test.mjs`
- `tests/web/history.test.mjs`
- `tests/web/feature-selection.test.mjs`

cover:

- JS/Python discovery parity
- circular/linear/unknown、override/reset
- duplicate ID、Automatic/explicit selector、record reorder
- file replace purge、async rediscovery retention
- inactive start round-trip
- one resource→multiple canonical records
- committed vs draft authority
- feature shortcutの±strand、multipart、origin-spanning、even length

#### Browser tests

- extend `tests/web/linear-multi-record.playwright.spec.js`
- add or extend Circular multi-record spec
- add focused feature shortcut scenario

browser acceptance:

1. mixed-topology GenBankを一つのLinear cardへloadし、recordごとの行を確認する。
2. circular rowだけstartを編集し、Generate後に左端feature/tick/popupを確認する。
3. Circular single/batch/gridの全 groupingで同じrecord controlsを確認する。
4. save/reload後もsource card count、record rows、override、inactive draft、committed diagramが一致する。
5. selector切替とduplicate IDで値が漏れない。
6. feature選択→5′ end/midpoint→Generateでsource popup座標が不変。
7. 1280 x 720と390 x 844で横overflow、row label、focus-visible、disabled reasonを確認する。

#### 完了条件

- file/card単位stateを経由せず、全record設定が正しいcanonical recordへ届く。
- draftとcommitted resultのauthorityが維持される。
- keyboardだけでcheckbox、number input、Reset、feature shortcutsを操作できる。

### Phase 8: CLI、package-root API、public contractを接続する

Status: pending

Depends on: Phase 2–6

#### CLI owner

- `gbdraw/io/cli_tables.py`
- `gbdraw/api/record_planning.py` の manifest adapter
- `gbdraw/circular.py`
- `gbdraw/linear.py`
- shared CLI normalization helper
- `tests/test_cli_tables.py`
- Circular/Linear selector/integration/session tests

#### CLI作業

1. `RecordsTableRow` 末尾へnullable topology/start fieldを追加し、positional compatibilityを維持する。
2. `_RECORDS_COLUMNS`、reader、cell error、manifest projectionを更新する。
3. `topology` はblank/auto/circular/linear、`display_start` はblank/positive intを受ける。
4. table rowのsame-file selectorsが異なる値を保持し、`order` sorting後もalignmentを維持する。
5. 両 parserへsingle-record convenience flagsを追加し、共通validatorでexact-oneを要求する。
6. records table/direct flags/session replayのmutual exclusionを明示する。
7. CLI requestとsession sidecarがschema6 display fieldを保持することをtestする。
8. parser helpをsourceから生成し、Circular/Linear option setの差をreviewする。

#### Package-root/public API owner

- `gbdraw/interface.py`
- package-root export module
- `gbdraw/api/__init__.py`
- `tests/test_python_interface.py`
- `tests/test_public_contract.py`
- `tests/fixtures/public_contract.json`
- `tests/test_api_library_usage.py`

#### Python作業

1. `RecordDisplayOptions` をpackage root/typed namespaceから同じclassとしてexportする。
2. `record_displays` keywordを両 draw functionへ追加し、length/type validationを共有helperに置く。
3. root adapterはin-memory `RecordInput.display`へ写すだけにする。
4. typed API exact-one/mixed-source recipeとroot in-memory recipeを実行testへ追加する。
5. public snapshotはhelperで再生成し、意図したsymbol/signatureだけを適用する。
6. default root callのSVG bytesが変わらないことをtestする。

#### 完了条件

- Web、CLI table、CLI single-record、typed API、root APIが同じtransform結果を返す。
- API/CLI adapterにtopology inference、modulo、fragment処理の複製がない。
- expected `ValidationError` の型と文面にrecord contextがある。

### Phase 9: documentation、Gallery session、visual evidenceを更新する

Status: pending

Depends on: Phase 3–8

#### Technical documentation owner

- `docs/REFERENCE/web-app.md`
- `docs/REFERENCE/command-line.md`
- `docs/CLI_Reference.md` generated help
- `docs/REFERENCE/input-formats-and-tsv-schemas.md`
- `docs/REFERENCE/python-api.md`
- `docs/REFERENCE/typed-requests.md`
- `docs/SESSION_COMPATIBILITY.md`
- `docs/REFERENCE/session-and-request-compatibility.md`
- `docs/SVG_SEMANTIC_HOOKS.md`
- `docs/RECIPES.md`
- `docs/RELEASE_NOTES_0.14.0b0.md`

#### Documentation decisions

1. public termは `display start`、UIの具体labelは12時/左端とする。`origin` をpublic identifierにしない。
2. source coordinatesが変わらないこと、SnapGeneとの相違、crop制約、RC default、HSP interpolationを明記する。
3. records tableの同一multi-record file例、single-record CLI flags、typed/root API recipeを実行する。
4. Linear ruler/definitionのwrapped notationを説明する。
5. semantic match hookが一つのmatchに複数fragmentを許すことを更新する。
6. session 41/request 6 と40/5 migration historyを正確に書く。

#### Tutorials/Gallery

- `docs/TUTORIALS/GUI/first-circular-genome-diagram.md` にdisplay-start操作を既存workflowを壊さないoptional stepとして追加する。
- Linear anchoringは新しいpublic pageを作らず、`docs/RECIPES.md` とWeb referenceへ統合する。
- control追加で既存operation cropがずれるGallery tutorialをinventoryする。
- current-writer asset契約に従い、全Gallery sessionをowner toolで41/6へ再生成する。
- tutorial screenshotを変更する場合は `web-gallery-screenshot-maintenance` とcapture recipeを使い、bitmapを手編集しない。
- owner-maintained `examples/gbdraw_social_preview.png` は変更しない。

#### Commands/evidence

```bash
python tools/update_cli_reference_help.py --check
python -m pytest tests/test_documentation_reference_contracts.py tests/test_documentation_contracts.py -q
python -m pytest tests/test_gallery_session_semantics.py tests/test_refresh_gallery_sessions.py tests/test_web_packaging.py -q
python tools/capture_gallery_tutorial_screenshots.py --all --check
```

対象tutorialを再撮影した場合は、そのexampleのclean capture、Gallery Playwright、desktop/mobile目視を追加する。

#### 完了条件

- docsに旧 `circular_top_coordinate` / `top_coordinate`、default=1、Circular-only semanticsが残らない。
- command、Python code、table recipeがclean directoryで実行できる。
- Gallery sessionsは41/6で、diagramとeditor stateを維持する。

### Phase 10: final verification、performance、diff audit

Status: pending

Depends on: all prior phases

#### Focused Python gate

```bash
python -m pytest tests/test_record_display_transform.py tests/test_api_requests.py tests/test_record_planning.py tests/test_api_request_render.py -q
python -m pytest tests/test_session_request_codec.py tests/test_api_session.py tests/test_session_compat.py tests/test_session_io.py -q
python -m pytest tests/test_cli_tables.py tests/test_python_interface.py tests/test_public_contract.py tests/test_api_library_usage.py -q
python -m pytest tests/test_circular_feature_paths.py tests/test_circular_ticks.py tests/test_circular_label_placement.py tests/test_circular_radial_labels.py -q
python -m pytest tests/test_linear_feature_paths.py tests/test_linear_label_placement.py tests/test_linear_definition_alignment.py tests/test_linear_track_layout.py -q
python -m pytest tests/test_annotations.py tests/test_circular_annotation_tracks.py tests/test_linear_annotation_tracks.py tests/test_feature_underlay_rendering.py -q
python -m pytest tests/test_depth_track.py tests/test_circular_conservation.py tests/test_linear_multi_record_comparisons.py tests/test_comparisons.py tests/test_collinearity.py tests/test_protein_colinearity.py -q
python -m pytest tests/test_interactive_svg_cli_format.py tests/test_web_feature_catalog.py tests/test_web_feature_metadata.py -q
```

実際に存在するtest filenameへPhase 0でcommandを調整する。新規testを別ownerへ統合した場合は重複commandを残さない。

#### Focused Web gate

```bash
TMPDIR=/tmp node --test tests/web/record-display-options.test.mjs tests/web/record-selector.test.mjs
TMPDIR=/tmp node --test tests/web/session-request.test.mjs tests/web/session-authority.test.mjs tests/web/session-draft-authority.test.mjs tests/web/session-resources.test.mjs
TMPDIR=/tmp node --test tests/web/history.test.mjs tests/web/feature-selection.test.mjs tests/web/match-sequences.test.mjs
npx playwright test tests/web/linear-multi-record.playwright.spec.js --project=chromium --workers=1
```

Circular/feature-shortcut用specを新設した場合は同じChromium projectで追加する。Node PlaywrightがなければPython Playwrightで同一scenarioを実行し、browser sandbox failureは権限付きで同じcheckを再実行する。

#### Broad gate

```bash
ruff check gbdraw/
TMPDIR=/tmp node --test tests/web/*.test.mjs
python -m pytest tests/ -v -m "not slow"
python -m pytest tests/test_output_comparison.py::TestOutputComparison -v
python tools/verify_gui_offline.py
python -m build
```

browser wheelが必要なtestだけ、事前に `python tools/prepare_browser_wheel.py` を実行する。deployable bundleを作る別作業でない限りcache-bust tokenを更新しない。

#### Performance check

1. default identity pathでtransform object以外のper-feature fragment allocationが増えていないことをprofileする。
2. representative Circular multi-record、Linear multi-record comparison、large depthを各3回測定する。
3. no-rotation defaultに10%を超えるmedian regressionがあれば原因を調べ、単に許容しない。
4. SeqRecord sequence bytesやfeature listのrotation copyを作っていないことをmemory/profileで確認する。

#### Visual check

- nondefault Circular: anchor featureが12時、radial/horizontal labels、ticks、GC/skew/depth、comparison ringが一致
- nondefault Linear:複数circular recordsを同じgeneへanchorし、features/tracks/ribbonsがseamを誤横断しない
- plus/minus arrow、multipart、annotation cap
- interactive feature/match popupがsource座標と一つのidentityを示す
- 1280 x 720 / 390 x 844 Web controls
- PNG/PDF等のraster/vector exportでclippingがない

#### Diff audit

1. production diffでtransform owner、topology resolver、Web projectionが一つずつであることを確認する。
2. testsで既存assertionを削った箇所とreplacement contractを確認する。
3. schema/session fixture、Gallery generated assets、user docsをproduction/test diffと分ける。
4. default reference SVGに差分がないことを確認する。
5. nondefault reference追加はsource recipeとvisual reviewを持つことを確認する。
6. public contract snapshotの変更が新class/export/signature/CLI optionsだけであることを確認する。
7. generated wheel、`dist/`、`gbdraw.egg-info/` がcommit対象に入っていないことを確認する。

#### 完了条件

- focused/broad/browser/documentation gateが通る。
- pre-existing failureは再現commandと本work packageとの差を記録する。
- reference、Gallery、session fixtureの由来と変更理由を説明できる。

## 8. Cross-surface matrix

`Supported` はそのsurfaceから非default値を表現し、同じtyped transformまで到達できることを意味する。

| Surface | Circular | Linear | Per-record mixed source | Persistence / replay | Final consumer |
| --- | --- | --- | --- | --- | --- |
| typed `gbdraw.api` | Supported | Supported | selector付きexact-one / uniform ALL | canonical/session | request plan transform |
| package-root Python | Supported | Supported | in-memory record順の`record_displays` | session builder使用時 | same request plan |
| CLI direct | exactly-one only | exactly-one only | No、tableを案内 | CLI session sidecar | same request plan |
| CLI records table | Supported | Supported | Supported、same source repeat可 | embedded table dependencies/session | same request plan |
| Web fresh | single/grid/batch | single/multi-row/comparison | Supported、Automatic expansion | editable draft + committed request | worker typed renderer |
| Web session replay | Supported | Supported | Supported | 41/6、40/5 migration | worker typed renderer |

このmatrixの一行だけ独自のtopology inference、coordinate math、fragment logicを持ってはならない。

## 9. Test scenario matrix

### 9.1 Topology/start/transform

| Case | Expected |
| --- | --- |
| detected circular + auto + start None | identity、control enabled |
| detected circular + explicit K | Kが12時/左端 |
| detected linear + auto + K | validation error |
| unknown + auto + K | validation error、explicit circularを案内 |
| unknown + override circular + K | accepted |
| detected circular + override linear + inactive K | requestはstart null、draftは保持可 |
| forward full record + start 1 | identity |
| reverse full record + start None | existing RC identity、先頭source base L |
| reverse full record + start 1 | source base 1が12時/左端 |
| region + start K | validation error |
| ALL default | per-record detection |
| ALL uniform override | every recordへfan-out、one invalidならcontext付きerror |
| heterogeneous same file | selector付きrecordsが各値を保持 |

### 9.2 Geometry consumers

| Consumer | Seam before | Seam inside | Seam at endpoint | RC |
| --- | --- | --- | --- | --- |
| 1 bp feature | one part | N/A | no zero width | source strand/display direction確認 |
| ordinary arrow | one part | two fragments、one biological head | artificial headなし | terminal correct |
| multipart/origin feature | order維持 |必要数へsplit | duplicate fragmentなし | stable identity |
| horizontal/radial label | same text | projected anchor/leader | collision stable | upright/placement |
| annotation line/bracket/band/highlight | one mark | fragments | cap重複なし | source label |
| GC/skew/depth | continuous path | two continuous paths | duplicate diagonalなし | display order |
| Circular HSP ring | one arc | split arc/path | zero arcなし | original endpoints |
| Linear HSP | one ribbon | common-t fragments | duplicate fragmentなし | plus/minus全組合せ |

### 9.3 Identity and metadata

- feature stable IDはstart変更前後で同じ。
- rendered feature fragmentはunique DOM IDを持ち、common biological identityへgroupされる。
- source HSPはfragment数に関係なく一つのmatch catalog entryを持つ。
- feature/match popupはsource location/endpointsを表示する。
- feature/match sequence downloadはsource spanを使う。
- duplicate record IDはrecord index/keyで別instanceとして扱う。
- annotation table、depth table、comparison tableのsource valuesをwriterが書き換えない。

### 9.4 Persistence and authority

- current v41/schema6 fresh save/load/regenerate
- ungenerated Web draft と committed render の分離
- v40/schema5 load/replay/save promotion
- v39/v33 supported migration
- 34–38 rejection
- schema5 Circular grouping single/grid/batch保持
- resource alias/rewrite後もper-record display alignment保持
- Reset、undo、redo、mode switch、selector switch
- failed/cancelled/stale Generate後もprevious preview/export維持

## 10. Risk register

| Risk | Detection | Mitigation / stop rule |
| --- | --- | --- |
| base/boundaryの1 bpずれ | D7 oracle、1 bp/cut endpoint/RC property test | local `+1/-1` patchを禁止しadapterへ戻す |
| RC default outputが変わる | start Noneのreference/RC fixture | defaultを1へ正規化しない |
| hidden direct coordinate mathが残る | `rg`で`360 *`、base/step、position/lengthをaudit | consumerをtransformへ移すまでfeatureをcompleteにしない |
| SeqRecordを回転してidentity/cacheが変わる | sequence object/content、feature hash、cache counter | physical rotationを禁止 |
| schema 5をlegacy groupingとしてdecodeする | schema5 single/grid/batch fixture | equalityでなくhistorical threshold |
| session 40へold migrationを再適用する | v40 authority fixture | `>=40` era predicate、current constant比較を禁止 |
| v40 load→Saveが41/5になる | no-Generate save test | pure schema promotionをwriter boundaryで必須化 |
| one source cardの設定が別recordへ漏れる | mixed/duplicate/selector-switch browser test | source UID + selector key、card stateへ置かない |
| HSP splitがquery/subject対応を壊す | common-t oracle、plus/minus cases | independent split/cartesian productを禁止 |
| gapありHSPをexactと誤解させる | docs/popup contract review | endpoint interpolationと明記。exact表現をしない |
| match fragmentがduplicate ID errorになる | catalog/interactive multi-element test | common match identity + fragment attributeを同時land |
| quantitative pathがseamを横切る | unsampled cut visual/structural test | endpoint補間してpathを分離 |
| annotation mode checkが残る | Linear circular record wraps-origin test | parse/resolve責務を分ける |
| identity fast pathが遅くなる | Phase 0/10 median、allocation profile | no-shift fast path、no sequence copy |
| UI追加でGallery cropがずれる | capture `--check` と visual diff | affected recipeから再撮影、手編集禁止 |
| reference更新で欠陥を隠す | default reference diff audit | default差分は実装修正。再生成で受け入れない |

## 11. Dependency and landing rules

```text
Phase 0
  -> Phase 1 transform
      -> Phase 2 typed planning
          -> Phase 3 schema/session/Web writer (atomic)
          -> Phase 4 Circular consumers
          -> Phase 5 Linear consumers
              -> Phase 6 comparisons + interactive grouping (atomic)
                  -> Phase 7 Web UI
                  -> Phase 8 CLI/Python surfaces
                      -> Phase 9 docs/Gallery
                          -> Phase 10 final gate
```

Phase 4とPhase 5の実装作業はPhase 2/3後に並行できるが、同じtransform contractを変更する場合は先にPhase 1 testsを更新する。

次の変更は分離してlandしない。

1. schema 6 Python codec、Web writer、session 41 constants/reader。
2. comparison t-space fragment、multi-element match catalog、interactive SVG grouping。
3. public root signature、export、contract snapshot、executable documentation。
4. current session writer bump、Gallery current-writer assets、packaging contract。

開発途中でCircularだけを内部testすることはできるが、Linear consumer、session、Web/CLI/API parityが未完の状態でuser-visible controlをreleaseしない。比較を含むLinearでsilent wrong geometryを返すfallbackは作らない。

## 12. Acceptance criteria

### Semantic contract

- [ ] `start_coordinate=None` は全既存入力で追加shiftなし。
- [ ] 明示Kは1-based source baseであり、Circularの12時/Linearの左端に同じ意味で配置される。
- [ ] source sequence、feature/annotation/depth/comparison coordinatesを再採番・書換えしない。
- [ ] diagram mode、detected topology、override、display startが別ownerである。
- [ ] unknown topologyはmodeからcircularと推定されない。
- [ ] crop + explicit startは一貫したerrorになる。
- [ ] reverse complement + explicit source anchorがD7の数式に従う。

### Per-record ownership

- [ ] 一つのfile内のcircular/linear/unknown recordsを別々に設定できる。
- [ ] Web Automatic expansionで一つのsource cardから複数canonical recordsへ正しく投影される。
- [ ] duplicate record ID、reorder、selector switchで設定が漏れない。
- [ ] Circular single/grid/batchとLinear single/multi-rowでrecord-transform alignmentが保たれる。

### Rendering

- [ ] feature、underlay、label、leader、tick、definition、GC/skew、depth、annotation、comparisonが同じtransformを使う。
- [ ] cut-crossing feature/annotation/seriesに負幅、全幅chord/diagonal、人工arrow/capがない。
- [ ] Circular comparison ringとLinear ribbon/curveがcutで正しくfragment化される。
- [ ] Linear HSP fragmentはcommon-t endpoint interpolation規則に従う。
- [ ] orthogroup/collinear alignmentはprojected feature centerを使う。
- [ ] default reference SVGは変更されない。

### Identity and interaction

- [ ] stable feature identityはdisplay startで変わらない。
- [ ] 一つのsource matchは複数pathでも一つのmatch identity/popupを持つ。
- [ ] popup、tick、annotation metadata、sequence actionはsource coordinatesを使う。
- [ ] SVG DOM IDは一意で、document-local referencesが全て解決する。
- [ ] semantic hooksはdocumented common ID + fragment contractに従う。

### Surfaces

- [ ] typed API、package-root Python、CLI direct、CLI records table、Webで同じ結果を作れる。
- [ ] Web checkboxはdetected/overrideを区別し、start controlを正しくlock/unlockする。
- [ ] feature 5′ end/midpoint shortcutがstrand/multipart/originを正しく扱う。
- [ ] CLI multi-record ambiguityはrecords tableを案内するerrorになる。
- [ ] public export/signature/help/contract snapshotが実装と一致する。

### Persistence

- [ ] current writerはsession 41/request 6のみを書く。
- [ ] session 40/schema5をpositive inputとして完全にreplayできる。
- [ ] v40 load→Save without Generateがv41/schema6になる。
- [ ] v40をpre-40 legacy migrationへ通さない。
- [ ] committed request、editable draft、webFiles authorityが混ざらない。
- [ ] Gallery current sessionsが41/6へowner toolで再生成される。

### Quality and release

- [ ] display start変更でLOSAT/cache probe/file readを増やさない。
- [ ] focused Python/Node/browser gatesが通る。
- [ ] non-slow suite、ruff、output comparison、offline check、buildが通る。
- [ ] desktop/mobileとnondefault final SVG/PNGを目視確認する。
- [ ] technical docs、recipes、session history、release notesがcurrent contractを説明する。
- [ ] generated artifacts、reference outputs、Gallery assetsのdiffを別々にreviewする。

## 13. Completion ledger

| Phase | Status | Evidence required before completion |
| --- | --- | --- |
| 0 Baseline/oracle | pending | baseline commands、fixture inventory、v40 provenance |
| 1 Shared transform | pending | property/oracle tests |
| 2 Typed planning | pending | plan/materialization tests |
| 3 Schema/session | pending | 41/6 writer、40/5 matrix、Web parity |
| 4 Circular | pending | focused tests、default parity、visual artifact |
| 5 Linear | pending | focused tests、wrapped ruler/series visual |
| 6 Comparison/interactive | pending | common-t/multi-element identity tests |
| 7 Web | pending | JS/browser/session/feature shortcut evidence |
| 8 CLI/Python | pending | public contract、table/direct/API parity |
| 9 Docs/Gallery | pending | executable docs、generated sessions/captures |
| 10 Final gate | pending | broad commands、performance、diff audit |

statusを `completed` にするには、そのPhaseのproduction changeだけでなく、列挙したevidenceが必要である。reference regeneration、manual screenshot、session version定数の変更だけを完了証拠にしない。

## 14. Definition of done

Work package Cは、complete circular recordを任意のsource baseへアンカーし、Circularでは12時、Linearでは左端に表示でき、その結果が全tracks/comparisons/interactionsで一致するときに完了する。

同時に、default出力、source coordinate、stable identity、LOSAT cache、v40/schema5 replay、mixed-record ownershipが維持されなければならない。一つでもconsumerが旧座標式を使う、または一つのsurfaceだけが別の意味を保存する状態では完了としない。
