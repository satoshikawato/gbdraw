# Linear Track Occupancy Layout Implementation Plan

- 作成日: 2026-07-19
- 対象バージョン: `0.14.0b0` 以降
- 設計スナップショット: branch `multi_records_per_line`, HEAD `9d0e386` の作業ツリー
- 状態: Phase 0–7実装済み
- 対象: Linear diagram の feature、GC content、GC skew、depth、annotation、comparison の縦方向配置

## 0. 実装状況（2026-07-19）

| Phase | 状態 | 実装範囲 |
|---|---|---|
| 0. Characterization | 実装済み | Middle + Separate Strands + Resolve、same-side custom slots、comparison、half-open interval、lane上限をfocused testで固定した。 |
| 1. Geometry primitives | 実装済み | `VerticalBand`、feature/label/numeric footprint測定、paint/reserve contractを追加した。 |
| 2. Per-record planner | 実装済み | Default/Customを共通slot intentへ正規化し、recordごとのmeasured footprintをordered cursorでpackする。 |
| 3. Renderer cutover | 実装済み | rendererはresolved originを使用し、track geometry metadata schema 2へpaint/reserve bandを出力する。 |
| 4. Record/comparison/canvas | 実装済み | record spacing、comparison endpoint、minimum corridor、canvas extentをrecord planから解決する。 |
| 5. Web/session migration | 実装済み | panel openとcustom enableを分離し、feature geometry schema v2、session version 33、v32 migration、canonical projection round-tripを実装した。legacy gallery sessionはv32 reader pathで検証する。 |
| 6. Legacy removal/docs | 実装済み | legacy offset ownerとsilent lane fallbackを削除し、Tutorial、CLI Reference、FAQ、release notes、Web helpを更新した。browser wheelはpackaging/browser検証用に再生成し、cache-bustだけをdeploy workflowに残す。 |
| 7. Constraint composition regression | 実装済み | X範囲付き`CollisionBand`とkind-pair policyで隣接rowの必要間隔を解き、body/comparison/definition制約を加算ではなく`max()`で合成する。comparisonは実際に境界を跨ぐ場合だけ有効で、single/multi-rowが同じsolverを使う。 |

Default/Custom、single/multi-record、comparisonは同じrecord-local vertical planを使用する。
欠損Depthはreserve bandだけを保持し、feature slotのheight/spacingはv2 contractで解決する。

2026-07-21の回帰修正では、canvas enclosureとrow spacingの責務を再度分離した。
左側definitionのようにplot bodyとXが交差しないbandは間隔を増やさない。一方、local headerのように
Xが交差するdefinition/body pairはdefinition clearanceを要求する。comparison exclusionは
body reserveに包含されるためcross-kindで再加算せず、activeなcomparison boundaryの
comparison/comparison制約だけが`comparison_height`の最小corridorを保証する。

## 1. 結論

この修正では、`Resolve Overlaps` に Middle 専用の Y 補正を追加しない。
Linear diagram の縦方向配置を、次の一方向 pipeline に統合する。

```text
feature lane assignment
        |
        v
horizontal record sizing
        |
        v
per-record footprint measurement
        |
        v
vertical track packing
        |
        v
record / comparison / canvas placement
        |
        v
SVG rendering from resolved coordinates
```

`Resolve Overlaps` の公開上の意味は、引き続き「genome coordinate 上で重なる feature を
追加 lane へ割り当てる」とする。ただし、lane 割り当て後の feature 全体が占有する縦方向の範囲を、
GC content、GC skew、depth、non-overlay annotation、comparison、record spacing、canvas sizing が
必ず参照する。

完成後は、次の配置式を複数箇所に持たない。

- feature lane 数から推測した record 高さ
- `cds_padding + vertical_padding` による固定 GC offset
- Middle / Above / Below ごとの comparison endpoint 特例
- custom track slots 用の後付け feature shift
- single-record-row と multi-record-row で異なる縦方向 extent 計算

すべての描画要素は、record ごとの immutable な vertical layout plan から解決済み座標を取得する。

## 2. 現状の問題

### 2.1 `Resolve Overlaps` の結果が下流 layout に伝わらない

feature-feature の genomic interval overlap は `gbdraw/features/tracks.py` で解決され、
各 feature に `feature_track_id` が割り当てられる。feature renderer はこの ID を使って追加 lane を描く。

一方、legacy Middle layout の GC content、GC skew、depth は、
`gbdraw/diagrams/linear/positioning.py` の固定式で軸から配置される。
この式は実際に増えた feature lane の下端を参照しない。

そのため、Middle、Separate Strands、Resolve Overlaps の組合せでは、reverse strand の2本目以降の
feature lane が GC content または GC skew の予約帯へ入る。

### 2.2 feature slot が高さ0として解決される

`gbdraw/diagrams/linear/track_slots.py` の `_slot_height()` は、features renderer に常に `0.0` を返す。
同じ side にある後続 slot の cursor は feature lane 数に応じて進まない。

assembly は `_custom_slot_track_offset_y()` と `_custom_record_top_extent()` /
`_custom_record_bottom_extent()` で一部を後付け補正しているが、次を一貫して扱えない。

- `side=overlay` の Middle feature
- feature と numeric track が同じ side にある任意の slot order
- record ごとに異なる feature lane 数
- feature label、leader line、ruler、stroke の突出
- comparison endpoint と custom slot の組合せ

### 2.3 `overlay` が二つの意味を持つ

現在の feature `side=overlay` は、実質的に次の二つを同時に意味する。

1. feature を record axis 上に配置する
2. vertical packing で領域を予約しない

Middle feature に必要なのは1だけである。軸上の feature lane は上下両側に paint されるため、
他の track に対しては領域を予約しなければならない。

annotation の明示的 overlay は、anchor slot と重ねること自体が意図された表現である。
この意味を feature の axis placement と同じ field に依存させない。

### 2.4 extent の owner が複数ある

現在は、次の値が別々の式で計算される。

- SVG group の Y transform
- record 間隔
- canvas height
- comparison ribbon の開始・終了 Y
- custom slot geometry metadata

見た目を修正しても別の式が古いまま残りやすく、clip、余分な空白、comparison overlap が再発する。

## 3. 設計上の不変条件

### 3.1 occupied-band invariant

各 renderer は、record axis を local `y=0` として次の二つの band を返す。

```python
@dataclass(frozen=True)
class VerticalBand:
    top_y: float
    bottom_y: float


@dataclass(frozen=True)
class LinearSlotFootprint:
    paint_band: VerticalBand
    reserve_band: VerticalBand
```

- `paint_band`: 実際に paint され得る全pixelの範囲。
- `reserve_band`: spacing と安全余白を含み、他要素が侵入してはいけない範囲。

`paint_band` には、nominal height だけでなく次を含める。

- fill path
- stroke の半幅
- line cap とmarkerの突出
- tick とtick label
- feature label とleader line
- annotation label
- ruler label

必須 invariant は次のとおりとする。

1. 同一 record の non-overlay `reserve_band` 同士は交差しない。
2. `z` は paint order のみに影響し、packing 結果を変えない。
3. `resolve_overlaps` は feature lane assignment のみに影響する。
4. lane assignment 後の feature footprint は、`resolve_overlaps` の値にかかわらず packing が消費する。
5. axis、feature、numeric track、annotation、comparison の座標は同じ plan から取得する。
6. slot `spacing` は reserve-band edge 間の距離として一度だけ適用する。
7. intentional overlap は明示的な overlay annotation とその anchor の間だけ許可する。
8. overlay annotation が anchor の外へ突出する場合、その突出を composite reserve band に含める。
9. canvas viewBox は全 `paint_band` と必要な title、legend、length bar、definition を包含する。
10. overlay annotationの`anchor_slot`はcomplete stack内の描画可能なnon-annotation slotへ解決する。unknown、non-drawing、annotation chainは拒否する。

### 3.2 placement とcollision policyを分離する

公開 input の互換性を維持しつつ、normalized model では次を分離する。

```text
anchor
  - axis
  - above
  - below
  - slot:<id>

collision policy
  - reserve
  - overlay
```

既存 input は次のように正規化する。

| Existing input | Internal anchor | Collision policy |
|---|---|---|
| Middle `features@side=overlay` | `axis` | `reserve` |
| `features@side=above` | `above` | `reserve` |
| `features@side=below` | `below` | `reserve` |
| non-overlay annotation | `above` / `below` | `reserve` |
| overlay annotation | `slot:<anchor_slot>` | `overlay` |

公開 CLI に新しい option を追加する必要はない。まず内部 contract を正しく分離する。

### 3.3 測定と描画で同じgeometry sourceを使う

feature lane の Y 座標は、`calculate_feature_position_factors_linear()` と同じ規則から一度だけ解決する。
renderer が同じ factor を独自に再計算しない。

推奨する内部型は次のとおりとする。

```python
@dataclass(frozen=True)
class LinearFeatureLane:
    strand_pool: Literal["shared", "positive", "negative"]
    track_id: int
    band: VerticalBand


@dataclass(frozen=True)
class LinearFeatureLaneGeometry:
    lanes: tuple[LinearFeatureLane, ...]
    occupied_band: VerticalBand
```

純粋関数 `measure_linear_feature_lanes()` が feature dict、`cds_height`、strand separation、
track layout、axis gap からこの geometry を返す。

## 4. Renderer別のfootprint contract

### 4.1 Axis and ruler

- axis stroke の半幅を含める。
- ruler tick length を含める。
- ruler label のfont bboxとbaselineを含める。
- axis/rulerはvertical packerのstructural seedとする。

### 4.2 Features and labels

- 実際の `feature_track_id` ごとに top / bottom を測定する。
- feature stroke の半幅を含める。
- external label、rotated label、leader lineをunionする。
- embedded labelはfeature glyphの外へ出る場合だけbandを拡張する。
- Middleはaxis固定のcomposite seedとして上下両方を予約する。
- Above / Belowは動的高さを持つ通常slotとして扱う。

### 4.3 Dinucleotide content and skew

- deviation content: `[-height / 2, +height / 2]` を基本 footprint とする。
- percent content: `[0, height]` を基本 footprint とする。
- skew: `[-height / 2, +height / 2]` を基本 footprint とする。
- plot stroke、baseline stroke、axis/tickがある場合はその突出を含める。

### 4.4 Depth

- plot bodyは `[0, height]` を基本 footprint とする。
- border、quantitative axis、tick、tick textを含める。
- recordごとにdata rangeが異なっても、同じresolved heightを使う場合はfootprint contractを共有できる。

### 4.5 Annotation

- lane layoutが返すrequired extentを基礎にする。
- bar stroke、label bbox、label offsetを含める。
- `overflow=clip` はclip bandと外側strokeまでをpaint bandとする。
- overlay annotationはanchorとのcomposite footprintを作る。

### 4.6 Spacer

- paint bandを持たない。
- 指定された距離だけcursorを進める。
- spacerとslot spacingを二重加算しない。

### 4.7 Definition, title, legend, and length bar

これらはtrack slotではなくcanvas extent ownerとして扱う。
record bodyとX領域が異なるdefinitionを、常にtrack stackへunionして過大なrow間隔を作らない。

次を別々に保持する。

- `record_body_band`: feature、label、axis、numeric、annotation
- `comparison_exclusion_band`: comparisonが侵入してはいけないrecord範囲
- `canvas_band`: definitionなどを含むcanvas bounds用範囲

## 5. Per-record vertical packing

### 5.1 defaultとcustomを同じinputへ正規化する

custom slotsがない場合も、内部では `default_linear_track_slots()` をmaterializeする。

```text
simple options ---------+
                       +--> normalized slot intents --> per-record packer
custom track slots -----+
```

これにより、legacy numeric track pathとcustom slot pathの二重実装を廃止できる。

### 5.2 recordごとに解決する

slot structureはdiagram全体で共有できるが、最終offsetはrecordごとに異なる。

- feature lane数が異なる。
- label数とlabel placementが異なる。
- annotation lane数が異なる。
- record-local ruler labelが異なる。

したがって、globalな `LinearTrackLayout` だけではなく、各recordについて
`LinearRecordVerticalPlan` を生成する。

```python
@dataclass(frozen=True)
class LinearResolvedSlot:
    slot_id: str
    renderer: str
    origin_y: float
    paint_band: VerticalBand
    reserve_band: VerticalBand
    z: int


@dataclass(frozen=True)
class LinearRecordVerticalPlan:
    axis_band: VerticalBand
    slots: tuple[LinearResolvedSlot, ...]
    record_body_band: VerticalBand
    comparison_exclusion_band: VerticalBand
    canvas_band: VerticalBand
```

### 5.3 packing algorithm

1. axis/ruler structural bandをseedする。
2. Middle featureはaxis固定compositeへunionする。
3. Above側はaxisに近いslotから外側へ解決する。
4. Below側もaxisに近いslotから外側へ解決する。
5. non-overlay slotを配置するたびにoccupied edgeを更新する。
6. overlay annotationをanchor解決後に配置し、composite footprintを更新する。
7. 最終unionから `top_extent` と `bottom_extent` を一度だけ導出する。

Below slotの基本制約は次のとおりとする。

```text
resolved_slot_top >= occupied_bottom + gap
```

Above側は鏡像とする。

### 5.4 非衝突時の互換性

初回移行では、既存のpreferred originを候補として保持し、衝突時だけ外側へ押す。

```python
new_origin = max(
    legacy_preferred_origin,
    occupied_bottom + gap - local_paint_top,
)
```

これにより次を両立する。

- 既にclearな出力のY座標を維持する。
- feature laneが増えたrecordだけnumeric trackを必要量移動する。
- reference SVGの不要な差分を抑える。

bugを残すための長期的なlegacy layout flagは追加しない。

### 5.5 custom slot order

Above / Below featureをaxis-adjacent以外へ移動できる仕様を維持するなら、featureとlabel group全体へ
resolved Y shiftを渡し、slot orderを実際の描画へ反映する。

初回実装でfeature group shiftを安全に実装できない場合は、対応できない順序を
`ValidationError` にする。featureを高さ0として黙って重ねる挙動は残さない。

## 6. Horizontal and vertical solver separation

multi-record rowでは、label placementが最終sequence widthに依存する。
一方、共有 `px_per_bp` はsequence lengthとrecord gapから決定でき、vertical extentには依存しない。

layoutを次の二段階へ分ける。

1. Horizontal plan
   - `px_per_bp`
   - record sequence width
   - record X position
   - horizontal inset
2. Vertical plan
   - final widthを使ったlabel placement
   - per-record slot footprint
   - axis Y
   - row spacing
   - comparison corridor
   - canvas top / bottom

既存 `gbdraw/layout/linear_multi_record.py` の `LinearRecordMeasurement` と
`LinearRecordPlacement` はrenderer非依存のaggregate contractとして維持する。
feature固有fieldを追加せず、`LinearRecordVerticalPlan` のunionから
`top_extent`、`bottom_extent`、comparison extentを渡す。

通常の1 record per rowも同じplacement contractへ段階的に統合する。
ただし、`placement is not None` が現在record-local rulerなどのpresentation mode判定にも使われているため、
geometry placementとpresentation flagを先に分離する。

## 7. Comparison placement

### 7.1 推奨contract

comparison endpointは各record planから導出する。

```text
comparison_top_y    = exclusion_band.top    - endpoint_gap
comparison_bottom_y = exclusion_band.bottom + endpoint_gap
```

recordの上下関係に応じてquery / subjectのtopまたはbottom endpointを選択する。

`comparison_height` は固定offsetではなく、二つのendpoint間に必要な最小clear spanとする。
不足する場合はrow axis spacingを増やし、comparisonをrecord bodyへ押し込んで縮めない。

次を同じendpoint contractへ送る。

- legacy ordered BLAST inputs
- explicit `LinearComparison`
- one-record-per-row layout
- multiple-records-per-row layout
- ribbon style
- curve style

### 7.2 intentional underlay

既存Above / Belowのaxis-to-axis ribbonを意図的なunderlayとして残す場合は、
Middle / Above / Belowの暗黙分岐ではなく明示的なcomparison collision policyとして表す。

推奨defaultはcontent-edge endpointである。比較リボンもrecord bodyへ侵入させない。
underlayを残す場合は、visual compatibilityとして明示し、通常のnon-overlap保証の例外としてテストする。

## 8. Feature lane allocatorの隣接修正

vertical plannerが正しくても、feature lane allocatorがoccupied laneへsilent fallbackすると
`Resolve Overlaps` の保証は成立しない。

同じ実装期間に次を修正する。

- 固定候補数を使い切った場合に最初のlaneへsilent fallbackしない。
- 必要に応じてlaneを動的追加するか、明示的なresource limitで`ValidationError`にする。
- genomic intervalはBiopythonのlocation semanticsに合わせてhalf-open `[start, end)` として扱う。
- endpointが接するだけのfeatureをoverlapとして扱わない。
- compound featureは初回は `min(start)..max(end)` envelopeを維持し、segment-aware packingは別変更とする。

lane selection algorithm自体の大幅な視覚変更は、vertical layout cutoverと同時に行わない。
既存のlane preferenceを可能な限り維持する。

## 9. Rendering integration

rendererとbuilderは、解決済みoriginを受け取る。

- record axis group
- feature and label group
- dinucleotide content group
- dinucleotide skew group
- depth group
- annotation group
- definition group

rendererは次を行わない。

- `cds_padding + vertical_padding` の再加算
- feature lane数からのrecord高さ推定
- slot sideからの独自Y shift
- comparison countからのrecord endpoint推定
- canvas heightの更新

featureとaxisを別originへ配置する必要があるcustom orderでは、`SeqRecordGroup` 内のaxis/rulerと
feature/labelを分離する。axisはrecord axisに残し、feature/labelだけをresolved feature originへ移動する。

## 10. Web and session contract

### 10.1 UI wording

Web UIでは次の表記を使用する。

- `Resolve Feature Overlaps`
- Help: overlapping genomic features are assigned to additional lanes; other tracks are repacked automatically.
- Middle: `Features on axis`

`Resolve Overlaps` がGC、depth、annotation、comparison同士の汎用collision solverであるような説明はしない。

### 10.2 Feature slot height and spacing

v2 contractでは次の意味とする。

- feature slot `height`: 実測feature bandの最小予約高、`max(actual, configured)`
- `--feature_height`: feature glyph自体の厚さ
- feature slot `spacing`: 外側の隣接trackとのclearance

Webのtooltipは、feature slot heightがglyph heightを変更しないことを明示する。
Auto表示はrecordごとに実測値が異なることを示す。

### 10.3 Custom Track Slots state

現在はCustom Track Slots panelのcaret操作がcustom layoutのenable/disableも兼ね、
再有効化時にsimple controlsからstackを再生成する。

次を別stateにする。

- custom stackを描画に使用するか
- panelを開いているか

要件:

- panel close/openでslot内容とCLI引数を変更しない。
- disable/re-enableで既存stackを復元する。
- `Reset` 操作だけがsimple controlsからstackを再生成する。
- custom modeではslot listを描画設定の正本とする。
- simple track controlsはdisableするか、simple mode onlyと表示する。

### 10.4 Session migration

linear track slot geometry schemaをv1からv2へ更新する。

v1ではfeatureのheight / spacingが実質的に無効だったため、v1 session読込時はこの二値だけを除去し、
旧sessionのeffective renderingを維持する。他の値は保持する。

- slot order
- renderer
- enabled
- side / axis index
- numeric track height / spacing
- renderer params

canonical-only sessionも識別できるよう、session versionは33とする。
version migration後にcanonical projectionを行い、旧versionのcanonical-only sessionを読み込めるようにする。

canonical requestのbuild側が保持するlinear slots / axisを、projection側でも対称に復元する。
complete sessionのvalidation/migrationを先に行い、成功後だけlive stateへcommitする。不正なimportは現在のstateを変更しない。

## 11. 実装phase

### Phase 0: Characterization and failing tests

描画コードを変更する前に、現在の挙動を固定する。

- screenshot相当のMiddle + Separate Strands + Resolve + GC + skew collisionをsynthetic recordで再現する。
- reverse strand lane `-1`, `-2`, `-5` のfeature bandを固定する。
- 現在clearな代表出力のgroup Yを固定する。
- Above / Below comparisonの現行axis-to-axis contractをcharacterizationする。
- custom same-side feature + numeric trackの現行collisionを再現する。
- Below + centered GCの境界食い込みを再現する。

このphaseではreference SVGを更新しない。

### Phase 1: Pure geometry primitives and shadow measurement

- `VerticalBand` とvalidationを追加する。
- feature lane geometry measurementを追加する。
- GC、skew、depth、axis/ruler、annotation footprint helperを追加する。
- 旧height計算と新measurementを並行実行できるようにする。
- 既存描画は変更せず、代表caseで差をassertまたは記録する。

### Phase 2: Per-record track planner

- default optionsをnormalized slotsへmaterializeする。
- recordごとに`LinearRecordVerticalPlan`を生成する。
- Middle featureをaxis composite seedとして予約する。
- Above / Below packerを共通cursor algorithmへ統合する。
- preferred originを維持し、collision時だけ外側へ押す。
- custom slotsとlegacy numeric tracksを同じplannerへ送る。

このphaseの完了時点で、今回のfeature-GC/skew/depth collisionを解消する。

### Phase 3: Renderer cutover

- numeric/annotation builderへresolved originを渡す。
- feature/label groupへresolved originを渡す。
- axis/rulerとfeature/labelのgroup ownerを必要に応じて分離する。
- geometry metadataをplanからserializeする。
- metadata schema 2へpaint / reserve boundsをadditiveに追加する。

### Phase 4: Record, comparison, and canvas unification

- horizontal planとvertical planを分離する。
- ordinary one-record-per-rowも`LinearRecordMeasurement` / `LinearRecordPlacement`へ寄せる。
- record spacingをrecord planのextentから解決する。
- comparison endpointをrecord planから解決する。
- `comparison_height`をminimum corridorとして扱う。
- canvas top / bottomを最終planから計算する。
- title、legend、length bar、definitionのextent ownerを明示する。

### Phase 5: Web and session migration

- UI wordingとhelpを更新する。
- custom enable stateとpanel open stateを分離する。
- feature height / spacing v2 semanticsを実装する。
- linear slot schema v1からv2へのmigrationを追加する。
- session version migrationとcanonical projectionを修正する。
- saved session、gallery session、Web CLI argument generationを検証する。

### Phase 6: Legacy removal and documentation

parity確認後に次を削除する。

- `_precalculate_feature_track_heights()` の並行list contract
- `_custom_slot_track_offset_y()`
- `_custom_record_top_extent()` / `_custom_record_bottom_extent()`
- `_record_plot_track_stack_bottom_y()`
- fixed `position_gc_content_group()` / `position_gc_skew_group()` / legacy depth offset
- `if non_middle_layout` によるtrack spacingとcomparison特例
- canvas configurator内の描画位置ownerとしてのlegacy plot offset fields

実装完了後に更新する文書:

- `docs/TUTORIALS/7_Linear_Layout.md`
- `docs/TUTORIALS/2_Comparative_Genomics.md`
- `docs/CLI_Reference.md`
- `docs/FAQ.md`
- release notes
- Web help textとgallery tutorial

## 12. Test plan

### 12.1 Pure unit tests

新しい `tests/test_linear_vertical_layout.py` を推奨する。

- `VerticalBand` のunion、translate、expand、validation
- deviation content / skewの対称footprint
- percent content / depthの片側footprint
- stroke、tick、large tick fontの突出
- axis/ruler visible / hidden
- Middle separated-strand lane `-1`, `-2`, `-5`
- Above / Belowの鏡像性
- slot spacingのedge-to-edge保証
- spacer
- `z` を変更してもplacementが変わらないこと
- overlay annotationとanchor composite
- recordごとにfeature lane数が異なるcase
- 全non-overlay reserve bandの非交差property

Hypothesisなどの新しい依存は追加せず、resolved planの隣接bandをparametrized testで検査する。

### 12.2 Python integration and SVG tests

最低限、次をcoverする。

1. Middle、Separate Strands、Resolve、reverse CDS 2 lanes、GC、skew。
2. reverse lanes 5本、複数depth、GC、skew。
3. Below + deviation GCの境界。
4. Above / Below、separate true / false、resolve true / false。
5. percent GC。
6. external / rotated labelsとleader line。
7. ruler on axis。
8. annotation above / below / overlay、label、clip。
9. custom axis featureと上下numeric track。
10. custom same-side order。
11. recordごとにlane数が異なる2から4 records。
12. legacy BLAST、explicit comparison、multi-record comparison。
13. ribbon / curve style。
14. legend、title、length bar、definitionのcanvas enclosure。

planだけのtestではmeasurementとrenderingが同じ誤りを共有できるため、代表caseではSVG elementの
absolute bboxも検証する。

### 12.3 Browser tests

Playwrightで次を検証する。

- dense Middle + Separate Strands + Resolve + GC/skewの`getBBox()`非交差
- custom overlay feature + above/below numeric/annotation
- panel close/openでstackとCLI argsが不変
- disable/re-enableでstackが復元される
- Reset時だけstackが再生成される
- v1からv2 session migration
- v2 round-trip
- v32 canonical-only sessionをv33 readerで復元
- save、load、generate後も同じslot specsになること
- SVG、PNG、PDF exportで同じlayoutが使われること

Node `@playwright/test` がない環境ではPython Playwrightでtargeted checkを行う。

### 12.4 Reference output and performance

- targeted tests通過後に`TestOutputComparison`を実行する。
- reference outputはgeometry changeを承認してから更新する。
- diffは主にY transform、canvas height、comparison pathへ限定されることを確認する。
- X coordinate、bp mapping、feature identity、group ID、colorを不変にする。
- plannerはordered cursorで `O(records * slots + features)` を維持する。
- slot全pair比較による `O(slots^2)` packingは導入しない。
- dense recordと複数recordで既存比のmaterialな性能低下がないことを測定する。

## 13. Acceptance criteria

実装は次をすべて満たしたときに完了とする。

1. screenshot相当のMiddle collisionが解消される。
2. feature lane数が増えると、そのrecordの必要なtrackだけが外側へ移動する。
3. non-overlay reserve band同士が許容epsilon内で交差しない。
4. feature、label、GC、skew、depth、annotation、comparisonが同じrecord planを参照する。
5. defaultとcustom track slotsが同じvertical packerを使う。
6. Middle、Above、Belowで共通packing codeを使い、差はfeature footprint生成に限定する。
7. single-record-rowとmulti-record-rowが共通placement contractを使う。
8. comparison corridorがrecord exclusion bandへ侵入しない。
9. canvasが全painted contentをclipしない。
10. 既にclearなlegacy outputのY位置を可能な限り維持する。
11. feature lane allocatorがsilent overlapへfallbackしない。
12. Web session migrationで既存slot orderとnumeric geometryを失わない。
13. Custom Track Slots panelの開閉が描画設定を変更しない。
14. focused Python、Web、browser、reference comparison testが通る。
15. 旧offset ownerと後付け補正が削除され、geometryのownerが一つになる。

## 14. 主なriskとmitigation

| Risk | Mitigation |
|---|---|
| reference SVGの広範なY差 | preferred originを維持し、collision時だけ外側へ移動する |
| label measurementとrenderingのずれ | 同じfont metric、baseline、resolved label dataを共有する |
| feature custom orderでaxisまで移動する | axis/ruler groupとfeature/label groupを分離する |
| placement object導入でruler modeが変わる | geometry placementとpresentation flagsを分離する |
| comparison ribbonが過度に短くなる | `comparison_height`をminimum corridorとしてrow spacingを拡張する |
| diagram全体が高くなる | recordごとの実測bandとmove-only-on-collision policyを使う |
| overlay annotationの突出 | anchorとのcomposite reserve bandへunionする |
| font差とantialiasing | stroke/2と小さいsafety epsilonをreserve bandへ含める |
| Web sessionがslotを失う | versioned migrationとcanonical round-trip testを先に追加する |
| old/new layout ownerが併存する | shadow phase後に明示的cutoverし、legacy fieldsを削除する |

## 15. Non-goals

- Circular layoutを同時に移行しない。
- 汎用2D collision solverを導入しない。
- SVGを描画した後にDOM bboxを読み、位置を反復修正するruntime layoutにしない。
- 全rendererを継承させるplugin frameworkやscene graphを追加しない。
- feature compound locationのsegment-aware lane packingを同時に実装しない。
- label routing algorithmを全面的に置き換えない。
- comparison edge bundlingや自動record並べ替えを追加しない。
- bugを維持する長期legacy flagを追加しない。

Linear track間の問題はY方向の一次元packingとして解決できる。小さいmeasurement contractと
per-record planに限定し、汎用layout frameworkへ拡張しない。

## 16. ファイル別の変更計画

### Core geometry and layout

- `gbdraw/layout/linear.py`
  - axis-relative feature lane geometryのsingle sourceを追加する。
- `gbdraw/layout/linear_multi_record.py`
  - horizontal / vertical planning境界を整理する。
  - aggregate record extentとcomparison endpointを受け取る。
- `gbdraw/diagrams/linear/track_slots.py`
  - dynamic per-record footprintを受け取るpackerへ置き換える。
  - feature height 0 contractを廃止する。
- `gbdraw/diagrams/linear/precalc.py`
  - feature dict、lane geometry、label measurementをrecord dataへまとめる。
- `gbdraw/diagrams/linear/assemble.py`
  - plan生成とrenderer orchestrationに責務を限定する。
  - legacy/custom/single/multi/comparisonの重複offset式を削除する。

### Rendering

- `gbdraw/diagrams/linear/builders.py`
  - resolved origin / placementを各groupへ渡す。
- `gbdraw/diagrams/linear/positioning.py`
  - fixed GC/skew/depth Y式を削除する。
- `gbdraw/render/groups/linear/seq_record.py`
  - axis/rulerとfeature/labelの配置ownerを分離する。
- `gbdraw/render/drawers/linear/features.py`
  - measured lane geometryを使用する。
- linear label modules
  - measured feature bandとresolved label placementを共有する。
- linear numeric and annotation groups
  - resolved slot origin以外のlayout式を持たない。

### Configuration and public adapters

- `gbdraw/tracks/linear.py`
  - feature axis placementとreserve semanticsを正規化する。
  - feature slot height / spacing v2 contractをvalidationする。
- CLI / API adapters
  - default settingsを同じnormalized slot intentへ変換する。
  - 公開option名は初回変更しない。

### Web and sessions

- `gbdraw/web/index.html`
  - wording、tooltip、custom enable switch、independent caretを更新する。
- `gbdraw/web/js/app/linear-track-slots.js`
  - panel state、v2 feature geometry、reset behaviorを実装する。
- `gbdraw/web/js/app/run-analysis.js`
  - v2 slot specsをPythonへ渡す。
- `gbdraw/web/js/services/config.js`
  - session migrationを追加する。
- canonical request/session projection modules
  - linear slot listとaxisを対称的にrestoreする。

### Tests

- `tests/test_linear_vertical_layout.py`
- `tests/test_linear_track_layout.py`
- `tests/test_linear_track_slots.py`
- `tests/test_depth_track.py`
- linear annotation and comparison suites
- multi-record comparison suites
- `tests/web/` session and slot tests
- targeted Playwright specification

## 17. 実装判断

今回の根本修正で優先する判断は次の三点である。

1. defaultとcustomを常に同じvertical plannerへ通す。
2. recordごとの実測feature bandを他track、record spacing、comparison、canvasが共有する。
3. 描画コードはplanを消費するだけにし、Y座標を再計算しない。

Middleだけに `record_heights_below` を加算する修正は、今回の再現例だけを止めることはできる。
しかしcustom overlay、Belowのcentered GC、tick text、annotation、comparisonで同じ問題が残るため、
根本修正としては採用しない。
