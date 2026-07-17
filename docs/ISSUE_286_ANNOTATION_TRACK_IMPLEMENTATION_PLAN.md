# Issue #286 Annotation Track Implementation Plan

- 作成日: 2026-07-17
- 対象バージョン: `0.14.0b0`
- 設計スナップショット: branch `cli_tsv_inputs`, HEAD `4f10f3c` の作業ツリー
- 対象 Issue: [#286 Region label](https://github.com/satoshikawato/gbdraw/issues/286)
- 状態: 設計合意済み、未実装
- 目的: 複数 feature または任意のゲノム領域を、line、bracket、band、hatch と任意ラベルで示す annotation track を Circular/Linear の第一級機能として実装する

## 1. 結論

Issue #286 は既存 feature label の拡張として実装しない。`RegionAnnotation` を独立した
domain object、`AnnotationSet` をデータ集合、既存の `CircularTrackSlot` / `LinearTrackSlot` を
配置契約とし、両 mode に `renderer="annotations"` を追加する。

最終構造は次のとおりとする。

```text
coordinate range / selected features / imported table
                         |
                         v
                 RegionAnnotation
                         |
                         v
                   AnnotationSet
                         |
        TrackSlot references one AnnotationSet
                         |
                         v
          target resolution and lane layout
                         |
            +------------+------------+
            |                         |
            v                         v
   Circular annotation track   Linear annotation track
            |                         |
            +------------+------------+
                         v
               static / interactive SVG
```

次の要件を初回実装の必須範囲とする。

1. annotation track を他の track と任意順に配置できる。
2. 複数の独立した annotation ring / row を同一図に配置できる。
3. annotation track 固有の `width` / `height` / `spacing` を TrackSlot で指定できる。
4. 座標範囲と複数 feature 選択を同じ描画経路へ正規化できる。
5. Circular、Linear、origin spanning、crop、reverse complement、multi-record を扱える。
6. line、bracket、solid band、hatch band と任意ラベルを扱える。
7. CLI、Python API、Web、session、interactive SVG が同じ contract を保持する。

この文書の phase はレビュー可能な実装順序を示す。Phase 途中の機能を「Issue #286 完了」とは
扱わず、全完了条件を満たした時点で一つの完成機能とする。

## 2. 設計上の判断

### 2.1 `FeatureObject` の subclass にしない

`FeatureObject` は strand、directionality、qualifier、feature track ID、block/intron shape を持つ。
Region annotation は feature と置換可能ではなく、feature factory、feature overlap、feature label、
feature legend に参加させる必要もない。ダミー strand や feature type を与える設計は採用しない。

### 2.2 既存 feature label placement に混ぜない

Region label は region mark に付随する明示的な文字列であり、feature label filtering、
embedded/external 判定、leader-line routing の対象ではない。共用するのは font metrics、text
measurement、text path などの低位 primitive に限定する。

### 2.3 独立 overlay ではなく正式な track renderer にする

annotation の任意順、複数 ring、専用 width/height は track layout の責務である。別の annotation
専用 positioning system は作らず、既存 slot normalization、axis boundary、radial/vertical layout、
track geometry metadata を拡張する。

### 2.4 renderer 名は `annotations` とする

公開 renderer 名は既存の `features`、`ticks`、`depth` と同じ粒度で `annotations` とする。
初回に受理する item は `RegionAnnotation` のみとし、point annotation や任意 SVG object は非目標とする。
将来 item type を追加しても renderer 名を変更しなくてよい一方、今回のために汎用 plugin framework は
導入しない。

### 2.5 data、style、placement、paint order を分離する

- `RegionAnnotation`: 対象領域、label、mark、item-level style override。
- `AnnotationSet`: annotation の集合と set-level default style / legend metadata。
- `TrackSlot`: set の参照、side、width/height、spacing、paint `z`。
- annotation layout: lane assignment と必要 extent。
- mode-specific drawer: 配置済み geometry から SVG を作る。

TrackSlot の並び順は既存 axis-index contract に従う物理順序、`z` は SVG paint order のみとする。
物理順序を `z` で変更しない。

### 2.6 非目標

- feature内容からoperon、gene cluster、phage moduleを自動推定しない。
- arbitrary SVG、freehand drawing、画像埋込をannotation itemとして受け付けない。
- GenBank/GFF3 featureの意味や保存形式をannotation編集のために変更しない。
- annotation一機能のために全track rendererをplugin/registry architectureへ置換しない。
- annotation labelとfeature callout labelを一つの汎用placement engineへ統合しない。

## 3. 公開データ契約

### 3.1 record selector

annotation target は既存 record selection と同じ ID / index の意味を使う。新しい selector 文法は
作らない。重複 ID を含む場合に index で一意に指定できることを必須とする。

```python
@dataclass(frozen=True)
class AnnotationRecordSelector:
    record_id: str | None = None
    record_index: int | None = None
```

既存 `RecordSelector` を循環依存なしに再利用できる場合は新型を作らず再利用する。再利用できない場合も
validation と正規化を共通 helper に置き、同じ判定を複製しない。

### 3.2 target

```python
@dataclass(frozen=True)
class CoordinateSpan:
    record: RecordSelector | None
    start: int
    end: int
    coordinate_space: Literal["source", "local"] = "source"
    wraps_origin: bool = False
    out_of_bounds: Literal["clip", "skip", "error"] = "clip"


@dataclass(frozen=True)
class FeatureSpan:
    record: RecordSelector | None
    selectors: tuple[FeatureSelector, ...]
    envelope: Literal["outer_bounds", "segments"] = "outer_bounds"
    circular_path: Literal["shortest", "forward", "reverse"] = "shortest"


RegionTarget: TypeAlias = CoordinateSpan | FeatureSpan
```

公開座標は1-based inclusiveとする。resolver 境界で0-based half-open segmentへ変換し、rendererに
公開座標を渡さない。

- `source`: 元 input の座標。crop / reverse complement後も既存 coordinate mapを通して表示位置へ写す。
- `local`: materialized record の `1..N`。
- `wraps_origin=False`: `start <= end` を必須とし、reverse complement上でも両値は方向ではなく範囲境界を表す。
- `wraps_origin=True`: `start > end` を必須としてCircularでoriginをまたぐ。Linearでは validation error。
- `outer_bounds`: 選択 feature 全体を囲む連続領域。
- `segments`: 非連続部分を保持する。同一 annotation ID の複数 markとして描く。

`FeatureSelector` は既存 feature selector value、stable feature ID、qualifier matchingのownerを再利用する。
Webで選択した featureだけのために別のfeature検索実装を作らない。

### 3.3 style

```python
@dataclass(frozen=True)
class HatchStyle:
    angle: float = 45.0
    spacing: float = 5.0
    color: str = "#666666"
    width: float = 1.0
    cross: bool = False


@dataclass(frozen=True)
class RegionAnnotationStyle:
    stroke: str = "#404040"
    stroke_width: float = 1.5
    stroke_dasharray: tuple[float, ...] = ()
    line_cap: Literal["none", "tick", "arrow"] = "tick"
    fill: str | None = None
    fill_opacity: float = 0.2
    hatch: HatchStyle | None = None
    label_color: str = "#202020"
    label_font_size: float | None = None
    label_orientation: Literal[
        "auto", "horizontal", "tangent", "radial", "arc"
    ] = "auto"
    label_position: Literal["center", "start", "end"] = "center"
    label_offset: float = 4.0
```

色は既存 color normalization、数値は有限値、opacity は `0..1`、spacing / width は正値として
入力境界で検証する。SVG attributeへ渡す前に再解釈しない。

### 3.4 annotationとset

```python
@dataclass(frozen=True)
class RegionAnnotation:
    id: str
    target: RegionTarget
    label: str = ""
    mark: Literal["line", "bracket", "band"] = "bracket"
    lane: int | None = None
    style: RegionAnnotationStyle | None = None
    legend_label: str | None = None
    metadata: Mapping[str, str] = field(default_factory=dict)


@dataclass(frozen=True)
class AnnotationSet:
    id: str
    annotations: tuple[RegionAnnotation, ...]
    default_style: RegionAnnotationStyle = field(
        default_factory=RegionAnnotationStyle
    )
    legend_label: str | None = None
```

同一 set 内の annotation ID は一意とする。異なる set で同じ annotation IDは許容するが、SVG DOM IDは
`record + slot + set + annotation` から生成して衝突を防ぐ。

style優先順位は次で固定する。

```text
application theme
    < AnnotationSet.default_style
    < annotation TrackSlot params override
    < RegionAnnotation.style
```

`None` はinherit、空文字列は明示的な無効値として扱わない。無効化が必要なfieldは型付き値で表す。

### 3.5 DiagramOptions

`DiagramOptions` のfield増加を一つに抑える。

```python
@dataclass(frozen=True)
class AnnotationOptions:
    sets: tuple[AnnotationSet, ...] = ()
    table: DataFrame | None = None
    table_file: str | None = None


@dataclass(frozen=True)
class DiagramOptions:
    ...
    annotations: AnnotationOptions | None = None
```

`sets` と table入力を同時に与えた場合のmergeは行わず validation errorとする。CLIはtableを
`AnnotationSet` へ変換してからtyped requestへ渡す。

## 4. TrackSlot契約

### 4.1 renderer追加

両 mode のsupported rendererへ `annotations` を追加する。

```text
Circular: features, ticks, numeric tracks, depth, conservation, annotations, spacer
Linear:   features, numeric tracks, depth, annotations, spacer
```

Annotation用の別TrackSpecは作らない。既存slotの共通fieldをそのまま使う。

- Circular: `side`, `radius`, `width`, `spacing`, `inner_gap_px`, `outer_gap_px`, `z`。
- Linear: `side`, `height`, `spacing`, `z`。

### 4.2 typed renderer params

raw `params` はnormalization時に専用contractへ変換する。

```python
@dataclass(frozen=True)
class AnnotationTrackParams:
    set_id: str
    lane_gap_px: float = 3.0
    padding_px: float = 2.0
    overflow: Literal["error", "compress", "clip"] = "error"
    show_labels: bool = True
    style_override: RegionAnnotationStyle | None = None
    anchor_slot: str | None = None
    layer: Literal["underlay", "foreground"] = "foreground"
```

- `set_id` は必須。存在しないset参照はerror。
- 同じsetを複数slotから参照してよい。
- 明示slot modeで参照されないsetはwarningとし、描画しない。
- `overflow=error` をdefaultとし、明示size不足をsilent clipしない。
- `compress` はlane gapとlabel offsetだけを下限まで縮小し、font sizeは変更しない。
- `style_override` は同じsetを異なるtrack viewで描く場合に使い、slot parserでは型付きstyleへ正規化する。

### 4.3 axis boundaryとside

既存axis-index contractを保持する。

- Circular: axis indexより前はoutside、以後はinside。
- Linear: axis indexより前はabove、以後はbelow。
- 同一side内の相対slot順を保持する。
- `z` は物理配置に影響しない。

現状のaxis-derived normalizationはCircularで非-tick overlay、Linearで非-feature overlayを拒否する。
`annotations` についてのみ次を追加する。

- adjacent annotation trackはaxis境界からsideを導出する。
- `side=overlay` は明示的な `anchor_slot` を必須とし、axis-indexによるside導出の対象外とする。
- overlay annotation slotはlayout extentを予約しない。
- `anchor_slot` は存在し、同じrecord layoutに属し、描画可能なbandを持たなければならない。
- `layer=underlay` はanchorより小さいpaint order、`foreground` は大きいpaint orderを要求する。
  矛盾する `z` はerrorにしてsilent補正しない。

### 4.4 default slot synthesis

annotation inputがありcustom track slot未指定の場合は、AnnotationSet順にdefault slotを合成する。

- Circular: feature/tick layoutと衝突しない既定のoutside slot。
- Linear: feature axisのabove slot。
- width/heightはauto。
- bandであってもdefaultではadjacent trackに描き、feature overlayは明示slotだけで有効にする。

既存annotationなしのdefault layoutとSVGを変更しない。

## 5. 入力と解決

### 5.1 annotation table

両 mode共通のTSV readerを追加する。CLI optionは `--annotation_table PATH` と
`--annotation-table PATH` を同じdestinationへ受ける。

必須column:

```text
set_id  id  mark
```

targetは各行で次のどちらか一方を必須とする。

```text
record  start  end  coordinate_space  wraps_origin
```

または

```text
record  feature_selector  envelope  circular_path
```

optional columns:

```text
label lane legend_label stroke stroke_width stroke_dasharray line_cap
fill fill_opacity hatch_angle hatch_spacing hatch_color hatch_width hatch_cross
label_color label_font_size label_orientation label_position label_offset
```

同じtable内の `set_id` でrowをgroup化し、first appearance順を保持する。unknown column、空ID、
重複ID、target混在、mode非互換値、非有限数値をrow/column付き `ValidationError` にする。

### 5.2 resolver

resolverはI/Oとrenderingの間に一度だけ実行する。

```python
@dataclass(frozen=True)
class ResolvedRegionAnnotation:
    id: str
    set_id: str
    record_index: int
    segments: tuple[tuple[int, int], ...]
    midpoint_bp: float
    span_bp: int
    label: str
    mark: str
    lane: int | None
    style: RegionAnnotationStyle
```

解決順序:

1. materialized recordへrecord selectorをbindingする。
2. `source` coordinateを既存coordinate mapでlocalへ変換する。
3. `FeatureSpan`をvisible featureに限定せず、materialized recordのfeature集合から解決する。
4. crop範囲外を `clip` / `skip` / `error` policyで処理する。
5. reverse complement後の方向とsegment順を正規化する。
6. Circular origin spanを最大2segmentに分割する。
7. empty spanを除外またはerrorにし、warningを構造化して返す。

feature visibilityやcolor設定によってRegionAnnotationのtargetが変わらないことを保証する。

## 6. Layout

### 6.1 lane assignment

lane計算はmode-independentなinterval packingをownerとし、pixel footprint計算だけをmode adapterへ渡す。

occupied intervalには次を含める。

- genomic mark span。
- bracket capの余白。
- measured label widthをbpへ換算した範囲。
- annotation間padding。

自動laneは決定的なgreedy interval coloringとする。sort keyはrecord、normalized start、normalized
end、set ID、annotation IDで固定し、入力containerのhash順に依存させない。明示laneは固定し、同laneの
衝突を `overflow` policyに従って扱う。

origin-spanning区間はCircular interval setとして衝突判定する。Linearでorigin spanを受け入れない。

### 6.2 required extent

```text
required extent
  = outer padding
  + lane count * lane extent
  + lane gaps
  + inner padding
```

- Circular `width=None`: required radial widthを使用。
- Linear `height=None`: required vertical heightを使用。
- 明示sizeが不足し `overflow=error`: slot ID、required、availableを含むerror。
- `compress`:許容最小gapまで縮小して再計算。
- `clip`:slot clipPathを作り、interactive metadataにclipped stateを残す。

### 6.3 layout contract

```python
@dataclass(frozen=True)
class PlacedRegionAnnotation:
    annotation: ResolvedRegionAnnotation
    lane: int
    mark_geometry: ...
    label_geometry: ...
    occupied_bounds: ...


@dataclass(frozen=True)
class ResolvedAnnotationTrack:
    slot_id: str
    set_id: str
    placements: tuple[PlacedRegionAnnotation, ...]
    required_extent_px: float
    occupied_bounds: ...
```

Circular radial layoutとLinear track layoutは `required_extent_px` を通常slotと同じpackingへ渡す。
annotation track用の第二layout engineをassembly側に作らない。

## 7. Rendering

### 7.1 mode-specific geometry

共通resolver/layoutの後、幾何だけを分ける。

- Linear: bpからx、laneからyへ変換し、line/bracket/rect bandを生成する。
- Circular: bpからangle、laneからradiusへ変換し、arc/bracket/annular sectorを生成する。

drawerはDataFrame、SeqFeature、TrackSlot parser、config dictを読まない。配置済みcontractとSVG groupだけを
受け取る。

### 7.2 mark

- `line`: spanに沿ったline/arc。
- `bracket`: lineにstart/end capを付ける。`line_cap=arrow` は方向が定義されたtargetだけ許容する。
- `band`: track lane全幅またはoverlay anchor band内の矩形/annular sector。
- hatch:共通pattern registryがSVG `<defs>` を一意に生成し、同styleをdeduplicateする。

originをまたぐline/bandは視覚的には連続するが、SVG pathは安全なsegmentへ分割してよい。

### 7.3 label

- Linear `auto`: horizontal、mark中央。
- Circular `auto`:十分なarc lengthならarc text、短い場合はtangent、反転が必要な半円ではpath方向を補正。
- label幅をlane layoutへ含め、描画後に別のcollision補正をしない。
- annotation labelはfeature `show_labels`とは独立し、trackの `show_labels` で制御する。
- labelなしannotationを正常系として扱う。

### 7.4 layerとSVG ID

adjacent trackは通常track paint orderに従う。overlayは次のphaseを持つ。

```text
annotation underlay
anchor track
annotation foreground
external feature labels
```

各mark/text groupへ次を付ける。

```text
data-gbdraw-annotation-id
data-gbdraw-annotation-set-id
data-gbdraw-annotation-track-id
data-gbdraw-record-id
data-gbdraw-record-index
```

DOM IDは同一setの複数track描画でも一意にする。

## 8. CLIとPython API

### 8.1 CLI

共通argument ownerにannotation table optionを追加し、Circular/Linearの双方へ渡す。CLI validationは
table syntaxとoption conflictを扱い、生物学的target解決は共通resolverへ委譲する。

使用例:

```bash
gbdraw circular \
  --gbk genome.gbk \
  --annotation_table annotations.tsv \
  --circular_track_slot 'regions:annotations@set_id=functions,w=28px' \
  --circular_track_axis_index 1 \
  -o annotated-genome
```

```bash
gbdraw linear \
  --gbk genome.gbk \
  --annotation_table annotations.tsv \
  --linear_track_slot 'operons:annotations@set_id=operons,h=24px' \
  --linear_track_axis_index 1 \
  -o annotated-locus
```

最終option spellingとslot scalar syntaxは実装済みparserからCLI helpを生成して確定する。ドキュメントを
parserより先にtruth sourceにしない。

### 8.2 API

`AnnotationOptions`、target/style/set型、table readerを `gbdraw.api` から公開する。低位
`assemble_*` へannotation引数を個別に増殖させず、`DiagramOptions.annotations` をcanonical入口とする。
互換用assembly関数へ必要な場合も一つのresolved annotation bundleとして渡す。

`CircularDiagramRequest` / `LinearDiagramRequest` はmode-specific target値を早期検証するが、record長が
必要なvalidationはmaterialization後へ残す。

## 9. Sessionとinteractive SVG

### 9.1 canonical request codec

`AnnotationOptions` と全frozen dataclassを明示的にencode/decodeする。`asdict()` 任せでtype discriminatorを
失わない。

```json
{
  "diagramOptions": {
    "annotations": {
      "sets": [
        {
          "id": "functions",
          "annotations": [
            {
              "id": "replication",
              "target": {
                "kind": "coordinateSpan",
                "start": 12000,
                "end": 18500
              },
              "mark": "bracket",
              "label": "Replication"
            }
          ]
        }
      ]
    }
  }
}
```

- targetとstyleにstable discriminatorを持たせる。
- annotation table fileをresourceとして埋め込む経路と、materialized setを保存する経路を混在させない。
- session versionを更新し、旧versionはannotationなしとして読み込む。
- round-trip後のtyped request equalityと再生成SVG semanticsを検証する。

### 9.2 interactive metadata

feature metadataとは別のannotation collectionを追加する。interactive enrichmentはannotation DOM IDへ
metadataを結び付け、同じsetを複数trackで描いた場合もtrack IDで区別する。

## 10. Web UI

### 10.1 state

annotation dataとtrack配置stateを分ける。

```javascript
annotationSets       // domain data and styles
circularTrackSlots   // placement
linearTrackSlots     // placement
```

track rowへannotation itemを埋め込まない。history/config/sessionは両者を独立にclone/serializeする。

### 10.2 annotation editor

focused modulesを追加し、`app.js`や`run-analysis.js`へ編集ロジックを集中させない。

```text
gbdraw/web/js/app/annotations.js
gbdraw/web/js/app/annotations/state.js
gbdraw/web/js/app/annotations/target-actions.js
gbdraw/web/js/app/annotations/style-actions.js
gbdraw/web/js/app/annotations/table-codec.js
```

必要なworkflow:

1. coordinate rangeを直接追加する。
2. SVGまたはfeature editorで複数featureを選び、region annotationを作る。
3. AnnotationSetを作成、rename、duplicate、deleteする。
4. line/bracket/band/hatch、label、style、laneを編集する。
5. annotation trackをCircular/Linear track editorで追加し、setを選択する。
6. drag-and-dropで任意順に並べ、width/height/spacingを編集する。
7. generated annotationをクリックしてeditorの該当itemを選択する。

### 10.3 track editor

両Web track moduleのsupported renderer、labels、default ID、normalization、serialization、preview extent、
duplicate/remove/reorderへ `annotations` を追加する。

- annotation set selectorをrenderer固有controlとして表示する。
- Circularはwidth、Linearはheightをauto/manualで表示する。
- adjacent sideとoverlay anchor/layerを切り替える。
- unknown `set_id`、明示size不足、mode非互換をgeneration前に表示する。
- wheel capability probeへannotation option/renderer supportを追加する。

### 10.4 session/history/config

JSON clone可能なplain dataへ正規化し、File objectをannotation stateへ保持しない。annotation table uploadは
読み込み時にdomain stateへmaterializeし、history/sessionはmaterialized dataを保存する。

## 11. Legend

直接labelを持つbracketを自動でlegendへ重複追加しない。legend entryは次の明示値だけから生成する。

1. `RegionAnnotation.legend_label`
2. なければ `AnnotationSet.legend_label`
3. どちらもなければ追加しない

同一caption/styleはdeduplicateする。solid/hatch previewをCircular/Linearの既存legend groupへ追加し、
annotation専用legend layoutは作らない。

## 12. 実装構成

想定ownerは次のとおりとする。最終file名は既存package境界との整合をレビューして確定する。

```text
gbdraw/annotations/
    __init__.py
    models.py
    io.py
    resolve.py
    layout.py

gbdraw/render/
    patterns.py
    drawers/circular/annotations.py
    drawers/linear/annotations.py
    groups/circular/annotations.py
    groups/linear/annotations.py
```

既存変更面:

- `gbdraw/tracks/circular.py`, `gbdraw/tracks/linear.py`
- `gbdraw/diagrams/circular/radial_layout.py`
- `gbdraw/diagrams/circular/assemble.py`
- `gbdraw/diagrams/linear/track_slots.py`
- `gbdraw/diagrams/linear/assemble.py`
- `gbdraw/api/options.py`, `requests.py`, `request_render.py`, `diagram.py`, `__init__.py`
- `gbdraw/session_request_codec.py`, `session_io.py`, CLI session binding
- `gbdraw/circular.py`, `gbdraw/linear.py`, common CLI args
- `gbdraw/render/interactive_context.py`, `interactive_svg.py`
- Web state、track modules、annotation modules、config/history/session/run-analysis/template

`FeatureObject`、feature factory、feature label candidateへannotation固有branchを追加しないことをレビュー条件とする。

## 13. 実装順序

### Phase 0: contractとcharacterization

1. annotationなしのCircular/Linear track normalization、axis-index、SVG hashを固定する。
2. current overlay rejectionとslot orderをcharacterization testにする。
3. typed request/sessionのtrack round-tripを固定する。
4. Web track reorder、duplicate、axis movementの現行testを固定する。

完了条件:

- annotation contract追加前の既定SVGが変わらない。
- 既存rendererのslot semanticsを変更せずtestで説明できる。

### Phase 1: domain、table、resolver

1. immutable modelsとvalidationを実装する。
2. annotation table readerとDataFrame adapterを実装する。
3. record binding、source/local coordinate、feature target、crop/RC、originを解決する。
4. `AnnotationOptions` とAPI exportを追加する。

完了条件:

- target種別、record選択、座標変換、全error policyのunit testが通る。
- rendererを呼ばずにresolved contractを検証できる。

### Phase 2: track normalizationとlayout

1. 両modeへ `annotations` rendererとtyped paramsを追加する。
2. axis-derived adjacent sideと明示overlay/anchorを実装する。
3. deterministic lane packingとrequired extentを実装する。
4. Circular radial/Linear vertical layoutへ通常slotとして統合する。
5. auto size、explicit size、overflow policyを実装する。

完了条件:

- 複数annotation slotの任意順と独立extentがgeometry snapshotで確認できる。
- annotationなしのresolved track geometryが不変。

### Phase 3: SVG rendering

1. Linear line/bracket/band/labelを実装する。
2. Circular arc/bracket/annular band/labelを実装する。
3. hatch pattern registryと全export互換を実装する。
4. adjacent/overlay、underlay/foreground、origin spanを実装する。
5. legend entryとinteractive DOM metadataを実装する。

完了条件:

- 全mark/style/modeのtargeted geometry testが通る。
- PNG/PDF/EPS/PSでhatchとopacityを視認・自動検証できる。
- 参照SVG差分は新規annotation fixtureまたは意図したfixtureだけに限定される。

### Phase 4: CLI、typed request、session

1. 共通CLI optionとtable resource bindingを追加する。
2. `DiagramOptions.annotations` を両builderへ転送する。
3. canonical request codecとsession version migrationを実装する。
4. static/interactive exportを含むCLI/API/session round-tripを検証する。

完了条件:

- 同じannotation tableからCLIとAPIが意味的に同じSVGを生成する。
- session保存・読込後にannotation set、track order、style、targetが失われない。

### Phase 5: Web UI

1. annotation domain state/editor/table importを追加する。
2. 両track editorへannotations rendererとset bindingを追加する。
3. feature複数選択からregion targetを作るworkflowを追加する。
4. instant preview、click selection、history、config、sessionを接続する。
5. capability probeとPyodide wheel pathを更新する。

完了条件:

- coordinate/feature両targetをWebだけで作成、編集、再生成できる。
- 複数ring/rowをdrag-and-dropし、width/heightを変更できる。
- saved sessionをWebとCLIの双方で再生成できる。

### Phase 6: documentationとrelease

1. CLI helpと `docs/CLI_Reference.md` を同期する。
2. table schemaをTutorial 5、Circular配置をTutorial 1、Linear配置をTutorial 7へ追加する。
3. Python API example、recipe、FAQ、release noteを更新する。
4. Circular/Linearの代表例とWeb workflowをGalleryへ追加する。
5. screenshotが必要なWeb手順はgallery screenshot skillに従ってcaptureする。

完了条件:

- copy-paste可能な座標target、feature target、複数track、overlay band例がある。
- 文書のoption名、default、table schemaが実装と一致する。

## 14. テスト計画

### 14.1 domainとI/O

- valid/invalid target union、duplicate ID、unknown set。
- 1-based inclusiveから0-based half-openへの変換。
- source/local、crop、RC、partial/outside span。
- unique/duplicate record IDとrecord index。
- feature envelope/segments、missing selector、hidden feature independence。
- table row/column error、style normalization、non-finite values。

### 14.2 layout

- overlapping/non-overlapping/explicit lane。
- label幅がspanより長い場合。
- deterministic order。
- Circular origin overlap。
- auto/explicit extent、error/compress/clip。
- multiple slots、same set multiple views、unknown anchor。
- axis-index order、side conflict、overlay reservationなし。

### 14.3 rendering

- mode × mark × side × label orientation。
- origin span、short arc、reversed text path。
- underlay/foreground paint order。
- unique DOM IDsとinteractive attributes。
- hatch pattern deduplication。
- legend deduplication。

### 14.4 integration

- CLI/API/session equivalence。
- Circular single/multi-record、Linear multi-record。
- annotationあり/なしのstatic/interactive SVG。
- custom track slotsとdefault synthesis。
- all export formats。
- reference-output comparison。

### 14.5 Web

- annotation set CRUD。
- table import validation。
- feature selection to region。
- circular/linear slot creation、duplicate、remove、reorder、axis move。
- set binding、auto/manual size、overlay controls。
- config/history/session round-trip。
- generated SVG click selection。
- targeted browser testによる実際のPyodide generation。

## 15. 検証コマンド

実装中はfocused testを追加し、最終的に次を実行する。

```bash
pytest tests/test_annotations.py -v
pytest tests/test_circular_annotation_tracks.py -v
pytest tests/test_linear_annotation_tracks.py -v
pytest tests/test_circular_track_slots.py -v
pytest tests/test_linear_track_slots.py -v
pytest tests/test_api_requests.py tests/test_api_request_render.py -v
pytest tests/test_session_request_codec.py tests/test_api_session.py -v
pytest tests/test_interactive_svg_cli_format.py -v
pytest tests/test_output_comparison.py::TestOutputComparison -v
pytest tests/ -v -m "not slow"
ruff check gbdraw/
```

Web unit testはannotation modules、track modules、config/history/sessionを対象に追加する。browser検証が
必要な時はNode PlaywrightとPython Playwrightの両方の利用可能性を確認し、sandbox launch failureは
必要な権限で同じcheckを再実行する。

reference outputを更新する場合は通常比較でactual差分を確認した後、意図したgeometry変更だけを
`--update-reference-outputs` で更新し、SVG差分をレビューする。

## 16. リスクと抑制策

### 16.1 track分岐の増加

現行track layoutはrendererごとの明示分岐を持つ。annotation追加だけを理由に全renderer registryへ
全面改修しない。各分岐ではannotation ownerへ委譲し、assemblyへmark/style分岐を漏らさない。

### 16.2 coordinate contractの混乱

公開座標、materialized record座標、SVG座標を一つのobjectで兼用しない。resolver後は0-based half-openを
唯一の内部contractとし、type名で段階を区別する。

### 16.3 label engineとの二重実装

feature label placementを再利用しない一方、font metricsとtext primitiveは共有する。annotation layoutは
track内lane、feature label layoutはtrack外calloutという責務境界を保つ。

### 16.4 session payload増加

同じtable fileとmaterialized annotationsを重複保存しない。Webではmaterialized dataをcanonicalとし、
CLI file inputはresource referenceかmaterialized setの一方だけを保存する。

### 16.5 Web module肥大化

track managerにはplacement/set bindingだけを置き、annotation CRUD/style/target操作をfocused moduleへ
分ける。`run-analysis.js` はserializationされた入力を渡すだけにする。

### 16.6 export互換

hatch、opacity、text pathをSVGだけで完了扱いしない。CairoSVGとbrowser PDF/PNG経路でfixtureを検証し、
unsupported SVG constructを避ける。

## 17. SOLID、KISS、DRY

### SOLID

- SRP: target resolution、lane layout、track placement、SVG drawing、UI stateを分ける。
- OCP: feature/label codeを変更せずtrack rendererを追加する。将来のannotation item追加はannotation owner内に閉じる。
- LSP: RegionAnnotationをFeatureObjectとして偽装しない。
- ISP: drawerへDataFrame、SeqRecord、app state全体を渡さない。
- DIP: assemblyはresolved track contractへ依存し、TSV/Web/GenBank由来の違いを知らない。

### KISS

- 新しいannotation用track positioning systemを作らない。
- markごとのsubclass hierarchyを作らず、型付きenumと小さいdrawer dispatchを使う。
- annotation一機能のために全renderer/plugin registryを導入しない。
- default synthesisはcustom slot未指定時だけに限定する。

### DRY

- Circular/Linearでtarget resolution、style validation、lane packing、pattern registryを共有する。
- annotation table、API、Webを同じdomain modelへ変換する。
- record selector、coordinate map、feature selector、track scalar parser、font metricsを既存ownerから再利用する。
- session/CLI/Webで別々のannotation schemaを持たない。

## 18. 完了条件

Issue #286を完了とするには、次をすべて満たす。

1. `RegionAnnotation` と `AnnotationSet` が公開typed APIとして安定している。
2. Circular/Linearの `annotations` TrackSlotを任意順に複数配置できる。
3. autoおよび明示width/heightがlayout、canvas bounds、Web UIに反映される。
4. coordinate/feature target、origin、crop、RC、multi-recordが正しく解決される。
5. line、bracket、band、hatch、labelが全export形式で成立する。
6. overlay underlay/foregroundとadjacent trackが予測可能なpaint/space contractを持つ。
7. CLI、API、Web、session、interactive SVGが同じannotation semanticsを保持する。
8. annotationなしの既存SVG、track order、label behaviorに非意図差分がない。
9. targeted、non-slow、Ruff、reference comparison、browser検証が成功する。
10. CLI reference、tutorial、recipe、API docs、release note、Galleryが実装と同期する。
