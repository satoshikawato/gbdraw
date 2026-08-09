# Issue #289 Circular Radial Label Implementation Plan

- 作成日: 2026-07-16
- 文書区分: 内部実装計画
- 対象 Issue: [#289 Label arrangement for circular mode: radial arrangement](https://github.com/satoshikawato/gbdraw/issues/289)
- 対象: Circular mode の feature label 生成、配置、描画、canvas、CLI、Python API、session、Web UI、全出力形式
- 状態: 計画
- 原則: radial placement を短いラベルまたは outer label に限定しない。既存の label contract 全体を最初から扱う

## 1. 結論

Circular feature label に、現行の水平配置 `horizontal` と並ぶ `radial` 配置を追加する。
`radial` は外部ラベルの文字列の長手方向を円の半径方向へ向ける配置であり、単なる SVG text の
回転オプションとしては実装しない。候補生成、回転矩形の衝突判定、角度順序の維持、leader line、
feature/track clearance、legend、canvas bounds を含む独立した layout strategy とする。

初回実装から次をすべて対象にする。

1. outer label と inner label。
2. 短い gene symbol と長い product/override label。
3. `auto`、`embedded_only`、`external_only` の全 label rendering policy。
4. `tuckin`、`middle`、`spreadout`、custom circular track slots。
5. single record と multi-record circular canvas。
6. `resolve_overlaps`、`separate_strands`、origin-spanning/multipart feature。
7. side/top/bottom/corner legend、plot title、center definition、tick/numeric/depth/conservation tracks。
8. SVG、interactive SVG、PNG、PDF、EPS、PS、Web preview/export。
9. CLI、typed config、Python API config override、canonical session、Web UI、config save/load。

選択されたラベルを、文字数、ラベル数、密度を理由に黙って省略、短縮、切断しない。radial layout
には固定の label-count cutoff や text-length cutoff を設けない。必要な場合は配置半径と canvas を
拡張する。ユーザーが指定した固定 track geometry と必要領域が数学的に両立しない場合だけ、重なった
図を返すのではなく、原因と必要量を含む `ValidationError` を返す。

後方互換性のため default は `horizontal` とし、既存の horizontal reference SVG は変えない。

## 2. 成功条件と不変条件

### 2.1 機能上の不変条件

- 同じ入力、同じ設定から常に同じ配置を生成する。
- `radial` は label filtering、label table、qualifier priority、feature visibility の結果を変更しない。
- embedded 判定は placement strategy から独立させ、`auto` では収まるラベルを従来どおり feature 内へ描く。
- `embedded_only` は外部配置を行わず、`external_only` は長さに関係なく全候補を radial layout に渡す。
- `--labels out` は全 external label を外側へ配置する。
- `--labels both` は current strand semantics に従って outer/inner を分け、inner label を radial に描画する。
- genome 上の feature 順序と label の循環順序を保ち、leader line の交差を防ぐ。
- 0/360 度境界を通常の隣接関係として扱う。
- 左半円では text direction を反転し、文字を逆さにしない。
- label text、label bbox、leader line、feature、track、legend、center definition の意図しない衝突を残さない。
- layout が返す authoritative bounds を canvas、legend、multi-record packing が共通利用する。

### 2.2 品質上の不変条件

- radial 専用 geometry を Circular、Linear、canvas の三箇所へ重複実装しない。
- `gbdraw/labels/circular.py` に mode 判定と radial heuristic を追加し続ける構造にしない。
- renderer が feature semantics や衝突解決を再計算しない。
- Web 側に Python layout の複製を作らない。Web は option を渡し、生成 SVG を表示するだけにする。
- 新しい外部 dependency、Web build step、server API を追加しない。
- radial mode を理由に既存 horizontal solver の出力を変えない。

## 3. 現状と技術的課題

### 3.1 現在の Circular label pipeline

現在は `gbdraw/labels/circular.py` の `prepare_label_list()` が次の責務を同時に持つ。

- feature から label text を選ぶ。
- text bbox を計測する。
- embedded/external を判定する。
- outer/inner を分類する。
- horizontal label を楕円上へ配置する。
- label-label、label-feature、leader-label、leader-leader collision を解く。
- renderer 用の mutable dict を組み立てる。

`rearrange_labels_fc()` 以下には horizontal bbox を前提とした配置 heuristic が集中している。
`gbdraw/diagrams/circular/assemble.py` は別途 horizontal bbox を再計算し、legend collision と canvas
拡張を行う。`gbdraw/render/drawers/circular/labels.py` は external label を無回転 text として描く。

### 3.2 再利用すべき既存機能

- `gbdraw/core/text.py`: font-aware text width/height 計測。
- `gbdraw/layout/spatial.py`: AABB broad-phase index。
- `gbdraw/labels/linear.py`: rotated rectangle corners、SAT collision の既存実装。
- `gbdraw/diagrams/circular/radial_layout.py`: feature lane と全 circular track の authoritative radial geometry。
- `gbdraw/diagrams/circular/assemble.py`: canvas expansion、legend placement、single/multi-record assembly。
- `gbdraw/render/groups/circular/labels.py`: leader を feature の背面、text を feature の前面へ描く二段階 layering。
- `gbdraw/web/js/app/feature-editor/label-actions.js`: external label の直前にある二本の line から feature anchor を取得する contract。

### 3.3 radial で新たに必要なもの

horizontal AABB では、回転した text の実際の占有領域を表せない。radial mode には次が必要である。

- text anchor と rotation を含む oriented rectangle。
- oriented rectangle の AABB、SAT intersection、segment intersection。
- feature midpoint の preferred angle と label の placed angle の分離。
- radial text extent と tangential text extent の分離。
- inner/outer の available radial arena と center/track reservation。
- label layout から canvas/legend まで共有する正確な bounds。
- 過密時に label を落とさず、必要半径を単調に増やす収束手順。

## 4. Public contract

### 4.1 CLI

Circular subcommand に次を追加する。

```text
--label_placement {horizontal,radial}
```

- default: `horizontal`。
- `--labels none` の場合も値の parse/validation は行うが、描画結果には影響しない。
- Linear subcommand の既存 `--label_placement {auto,above_feature}` と同名でよい。subcommand parser が
  分離されているため user-facing contract は衝突しない。
- Circular CLI は `modify_config_dict(..., circular_label_placement=...)` に渡す。
- `--label_rendering`、`--label_font_size`、`--circular_label_spacing`、filter/table options はそのまま併用する。
- 現行 CLI の `--labels both` は GC/skew track を一律に suppress する。`horizontal` では互換性のため
  現行挙動を維持するが、`radial` ではこの suppress を行わず、inner solver が GC/skew を含む実際の
  track reservation と共存させる。mode 分岐は CLI、solver、assembler に重複させず、正規化済み
  placement と track request から一度だけ決定する。

使用例:

```bash
gbdraw circular \
  --gbk tests/test_inputs/NC_001879.gbk \
  --labels both \
  --label_placement radial \
  --qualifier_priority examples/qualifier_priority.tsv \
  -o tobacco_plastid_radial
```

### 4.2 Typed config

`labels.linear` と対称な typed section を追加する。

```toml
[labels.circular]
placement = "horizontal"
```

```python
@dataclass(frozen=True)
class LabelsCircularConfig:
    placement: Literal["horizontal", "radial"]
```

- 古い config/session に `[labels.circular]` がなくても `horizontal` として読む。
- normalization と validation は一箇所に置く。
- `LabelsConfig` に散在する既存の arc factor は horizontal strategy の設定として維持する。
- radial strategy は font size、spacing、stroke、feature clearance、track geometry を共有し、horizontal
  専用の ellipse center/arc-angle heuristic を暗黙に流用しない。

### 4.3 Python API

既存の API surface を増やし過ぎない。次の二経路を正式に対応させる。

```python
config_overrides={"circular_label_placement": "radial"}
```

```python
cfg.labels.circular.placement == "radial"
```

`modify_config_dict()` に `circular_label_placement` を追加し、既存の `label_placement` は Linear 用の
backward-compatible key として維持する。`DiagramOptions.config_overrides`、single-record builder、
multi-record builder、request renderer の全経路を test する。

### 4.4 Session と Web config

canonical `diagramOptions.configOverrides` に次を追加する。

```json
{
  "circular_label_placement": "radial"
}
```

- schema 1 の additive override とし、session version を placement 一項目だけのために上げない。
- field がない既存 session は `horizontal`。
- legacy CLI session import/export は Circular の `--label_placement` と相互変換する。
- Web state は Linear の `adv.label_placement` と混同せず、`adv.circular_label_placement` を所有する。

### 4.5 Web UI

Circular labels が有効なとき、Labels panel に native `<select>` を追加する。

```text
Circular label placement
  Horizontal
  Radial
```

- default は Horizontal。
- state、config save/load、canonical session、undo/history、run-analysis args の全経路へ接続する。
- Pyodide wheel capability detection に `--label_placement` の Circular support を追加し、古い wheel で
  Radial を選んだ場合は明確な rebuild/redeploy error を出す。
- horizontal 専用 X/Y arc offset controls は Radial 選択時に disabled/hidden とし、誤って効いているように見せない。
- radial layout の preview-only JavaScript 実装は作らない。

## 5. Target architecture

```text
feature_dict + label/filter/override policy
                 |
                 v
       CircularLabelCandidate[]
       (semantic data + text metrics)
                 |
                 v
      circular track/feature layout
                 |
                 v
    embedded/external + outer/inner projection
                 |
                 v
        placement strategy dispatch
          /                     \
 horizontal strategy       radial strategy
          \                     /
                 v
       CircularLabelLayoutResult
       - placed labels
       - oriented boxes/AABBs
       - leader segments
       - content bounds
       - required radius growth
                 |
                 v
     obstacle/canvas/legend resolution
                 |
                 v
       shared Circular label renderer
```

### 5.1 Module ownership

計画上の module 境界は次とする。実装時に既存 import graph を確認し、同じ責務を保てるなら細かな
filename は調整してよいが、責務を再結合しない。

| Module | 単一責務 |
|---|---|
| `gbdraw/labels/circular.py` | public/internal facade、strategy dispatch、既存 import compatibility |
| `gbdraw/labels/circular_candidates.py` | feature から semantic candidate を一度だけ生成 |
| `gbdraw/labels/circular_horizontal.py` | 現行 horizontal placement behavior |
| `gbdraw/labels/circular_radial.py` | radial placement、半径要求、order-preserving relaxation |
| `gbdraw/layout/text_geometry.py` | rotated box、AABB、SAT、segment/polygon intersection |
| `gbdraw/labels/circular_types.py` | small immutable contexts/results または TypedDict contract |
| `gbdraw/render/drawers/circular/labels.py` | placement result を SVG text へ変換 |
| `gbdraw/render/groups/circular/labels.py` | layering と二本の leader segment の描画 |
| `gbdraw/diagrams/circular/assemble.py` | track/label layout の orchestration、legend/canvas integration |

class hierarchy は作らない。二つの strategy は同じ callable contract を実装し、placement 名から
関数を引く小さな dispatch map を使う。

### 5.2 Candidate contract

Candidate は placement に依存しない情報だけを持つ。

- stable feature identity と input order。
- label text、font family、font size、measured width/height。
- feature type、strand、track id、directional flag。
- longest/origin-spanning segment の start/end/midpoint/span。
- current rendering policy と embedded eligibility に必要な値。
- feature coordinates。pixel radius は track layout 後に projection する。

Candidate 生成後に horizontal/radial それぞれで `get_label_text()` や font measurement を繰り返さない。

### 5.3 Placement result contract

各 external label は少なくとも次を持つ。

- `preferred_angle_deg` と `placed_angle_deg`。
- `text_x`、`text_y`、`rotation_deg`、`text_anchor`、`dominant_baseline`。
- `width_px`、`height_px`。
- oriented box corners と authoritative AABB。
- `feature_anchor_x/y`、`elbow_x/y`、`leader_start_x/y`。
- `is_inner`、`is_embedded=False`、feature/track metadata。

embedded label は従来の arc text path contract を保つ。renderer は placement strategy を推測せず、
result に含まれる text geometry をそのまま描く。

## 6. Shared text geometry

### 6.1 抽出対象

`gbdraw/labels/linear.py` にある次の処理を `gbdraw/layout/text_geometry.py` へ抽出する。

- anchor に対する未回転 rectangle corner。
- anchor 回りの corner rotation。
- corner union からの AABB。
- separating axis theorem による oriented rectangle intersection。
- point/segment と convex rectangle の intersection。
- padding を含む rectangle expansion。

Linear は shared helper を呼ぶように変更し、抽出前後で全 Linear label test の結果を一致させる。
Circular horizontal、Circular radial、canvas、legend は同じ geometry owner を使う。

### 6.2 座標と baseline

- layout 内部の box は font metric の width/height を使い、SVG text の baseline 差を一箇所で吸収する。
- text anchor は `start`、`middle`、`end` の三種類を扱う。
- corner order は一貫した clockwise/counter-clockwise とし、SAT と leader contact の双方で使う。
- geometry は DOM の `getBBox()` に依存しない。CLI、Pyodide、test で同じ結果にする。
- `calculate_bbox_dimensions()` の既存 CLI/Pyodide font fallback contract を維持する。

## 7. Radial placement specification

### 7.1 Coordinate system

Feature midpoint の preferred angle を次で定義する。

```text
theta = 2π * middle_bp / genome_length - π/2
```

`theta=−π/2` が 12 時、`theta=0` が 3 時である。feature anchor は actual feature lane の local
inner/center/outer radius から求める。base canvas radius の推測値を使わない。

### 7.2 Text orientation

- text baseline は placed angle の radial axis と平行にする。
- 右半円では読み方向が通常となる outward/inward anchor を選ぶ。
- 左半円では rotation を 180 度反転し、最終 rotation を読みやすい半平面へ正規化する。
- anchor を反転しても rectangle が outer label では feature から外側、inner label では中心側へ伸びるようにする。
- 12 時/6 時付近は epsilon 付きの一貫した tie-break を使い、微小な入力差で 180 度反転しない。
- orientation と box corner は同じ pure function から生成する。

### 7.3 Initial radial arenas

Outer label の最小開始半径は次の最大値とする。

- configured label anchor radius。
- その genome interval に存在する feature の最大 outer radius。
- feature/text clearance と leader safety margin。
- 外側 circular track の reserved radius。

Inner label の最大開始半径は次の最小値とする。

- configured inner label anchor radius。
- その feature の最小 inner radiusから clearance を引いた値。
- inner track/definition reservation が許す outer edge。

固定の一枚の円として扱わず、feature interval ごとの local track radius を使う。`resolve_overlaps` で
feature が複数 lane に分かれた場合も、leader は実際の lane edge から始まる。
inner side に GC/skew/depth/comparison/custom track の完全な annulus がある場合は、reservation を跨ぐ
一つの空間とみなさず、feature edge から中心へ向かう連続した available interval を列挙する。label text
は一つの interval 内に収め、足りなければ 7.7 節の radius preflight でその interval を広げる。

### 7.4 Cyclic order

1. preferred angle で stable sort する。同角度は stable feature identity/input order で決める。
2. 最大の angular gap を seam に選び、循環列を一時的に線形化する。
3. label の順序を変えずに forward/backward relaxation する。
4. 最後に seam pair も通常の neighbor として検査する。
5. placed angle の cyclic order と feature midpoint の cyclic order が一致することを assertion/test で保証する。

この順序保存を leader crossing 防止の第一条件とする。pairwise intersection test は最終検証として使う。

### 7.5 Collision-free packing

各 side について次の決定論的手順を使う。

1. candidate を preferred angle、最小 radial offset に置く。
2. oriented box の AABB を broad phase index へ登録する。
3. 衝突した contiguous angular cluster だけを取り出す。
4. cluster 内で角度を最小 spacing まで forward/backward に緩和する。
5. 角度だけで解けない場合は、cluster の radial start を外側（inner は中心側）へ移す。
6. 全体 circumference が不足する場合は必要半径を実測値から計算し、blind fixed-step ではなく不足量だけ拡張する。
7. SAT で text-text collision、segment/polygon で leader-text collision、segment intersection で leader-leader collision を再検査する。
8. collision-free な範囲で angle drift、radial displacement、leader length を縮める polishing pass を一度行う。

候補評価の優先順位は固定 tuple とする。

```text
(
  text_collision_count,
  leader_text_collision_count,
  leader_crossing_count,
  order_violation_count,
  total_angular_drift,
  max_angular_drift,
  total_radial_displacement,
  total_leader_length,
  stable_tie_break,
)
```

randomness、set iteration order、wall-clock timeout を選択結果に使わない。

### 7.6 Radial tiers

添付例のように dense cluster だけを feature から遠ざけ、疎な label は近くに保つため、start radius
は label ごとに持つ。ただし任意の二次元 bin packing にはしない。

- 同じ angular cluster 内で必要なときだけ離散 tier を増やす。
- tier gap は text spacing と stroke safety から導く。
- 新しい tier を追加したら cluster 全体の order と leader collision を再検査する。
- outer side の tier 数に固定上限を置かない。
- 最終 fallback は cluster 共通半径を必要量まで広げるため、全て有限サイズの label に配置可能解がある。

### 7.7 Inner labels と必要 radius の収束

Inner label は中心側の有限領域を使うため、長文と dense set を正しく扱うには track layout と協調する。

1. semantic candidates を track layout 前に生成する。
2. provisional circular radial layout を解く。
3. candidate を feature lanes へ projection し、embedded/external を確定する。
4. inner radial solver が center reservation、inner tracks、text width、spacing から必要な axis radius を返す。
5. 不足があれば `canvas_config.radius` と dependent radial layout を exact deficit だけ拡張し、projection/placement を再実行する。
6. radius は単調増加させ、同じまたは小さい要求へ戻さない。
7. requirements が変化しなくなった時点で geometry を freeze し、その後に SVG group を描く。

この preflight は SVG element を追加する前に完了させる。描画後の全 group scaling で辻褄を合わせない。

明示的な絶対 radius/center reservation が移動不可で、必要な inner interval と矛盾する場合は、
該当 slot ID、available px、required px を含む `ValidationError` を返す。inner label を黙って outer
へ移す、重ねる、削除する挙動は採用しない。

### 7.8 Long labels

- text を省略記号で短縮しない。
- 自動 wrap、改行挿入、font-size 縮小を行わない。
- finite な一行 text は計測幅をそのまま radial extent とする。
- outer label では canvas を必要量まで拡張する。
- inner label では前節の preflight で circle/track radius を広げる。
- Unicode、italic を含まない現在の plain external label contract、空白、記号を既存 text metric で扱う。
- multiline/rich-text label 自体は現行 feature-label contract に存在しないため、#289 で別形式を導入しない。

### 7.9 Leader lines

既存 Web editor contract を維持するため、external label ごとに二本の SVG `line` を保つ。

```text
feature anchor -> elbow -> text-facing box edge
```

- elbow は placed angle の radial ray 上へ置く。
- text contact は oriented box の feature-facing edge 上で、corner を避ける。
- outer/inner と text direction flip の全組合せで同じ geometry helper を使う。
- feature-to-elbow segment は feature annulus を横切らない。
- elbow-to-text segment は他 label の oriented box interior を横切らない。
- zero-length segment は許容するが、NaN/inf は許容しない。
- layering は leader phase -> feature/track phase -> text phase を維持する。

## 8. Canvas、legend、multi-record integration

### 8.1 Authoritative bounds

`CircularLabelLayoutResult` は全 external label の oriented AABB と leader endpoints の union を local
coordinates で返す。次がこれを共有する。

- label-label/legend collision。
- `_expand_canvas_to_fit_external_labels()`。
- top/bottom legend の content bounds。
- side/corner legend obstacle resolution。
- multi-record subcanvas inset/visible-content packing。
- clipping regression tests。

`assemble.py` と `labels/circular*.py` が別々に horizontal/radial bbox を再計算する状態をなくす。

### 8.2 Legend

- top/bottom legend は final label content bounds の外へ、既存 minimum gap を保って配置する。
- side/corner legend は local legend AABB を forbidden obstacle として strategy-specific reflow に渡す。
- radial reflow は collided angular cluster だけを移動し、rotation、box、leader、bounds を同時更新する。
- label を移動できない場合は legend/canvas を外へ拡張する既存 fallback を使う。
- legend collision 修正後に古い leader anchor を残さない。

### 8.3 Multi-record canvas

single-record subcanvas は final label bounds を drawing metadata として持つ。multi-record assembler は
SVG text の `rotate(...)` を不完全に推測せず、その metadata を既存 content estimator と merge する。

- mixed record lengths と harmonized label font size。
- per-record scaling mode `auto`、`linear`、`equal`。
- row/column gap ratio と explicit row placement。
- shared legend と plot title。
- record ごとに異なる radial extent。

これらで label 同士または隣接 record が重ならないことを検証する。

### 8.4 Export

External radial text は plain SVG `<text>` に anchor point 回りの `rotate(angle x y)` transform を付ける。
`textPath` は embedded label にのみ使う。

- DOMPurify が `transform` と必要 data attributes を保持すること。
- CairoSVG が SVG/PNG/PDF/EPS/PS で同じ方向を描くこと。
- Browser canvas PNG と svg2pdf.js PDF が clipping せず text rotation を維持すること。
- interactive SVG の feature popup/label editor が rotated text の CTM を正しく扱うこと。

## 9. SOLID、KISS、DRY の適用

### 9.1 SOLID

- SRP: candidate extraction、placement、geometry、rendering、canvas orchestration を分離する。
- OCP: placement strategy は共通 input/result contract に追加する。assembler/renderers に mode ごとの
  field inference を増やさない。
- LSP: horizontal と radial の result は同じ consumer から扱え、bounds/leader contract を満たす。
- ISP: solver に app 全体の config dict や canvas object を渡さず、必要な immutable layout context、
  candidates、obstacles だけを渡す。
- DIP: assembly は具体的な radial heuristic ではなく placement callable に依存する。class hierarchy や
  dependency-injection framework は導入しない。

### 9.2 KISS

- public switch は `horizontal|radial` の一つだけから開始する。
- 汎用 constraint solver、random optimizer、physics simulation を導入しない。
- cyclic sort、cluster relaxation、radial tier、exact bounds expansion の決定論的手順に限定する。
- existing font metrics、spatial index、track layout、two-line leader layering を使う。
- config knobs を speculative に増やさない。spacing/font/stroke/filter/track options は既存値を使う。
- hidden fallback を作らない。解けない explicit constraint は actionable error にする。

### 9.3 DRY

- rotated rectangle math は Linear/Circular で共有する。
- `get_label_text()` と bbox measurement は candidate builder の一回だけ。
- label bounds は collision、canvas、legend、multi-record で一つの result を共有する。
- placement normalization は config/CLI/session/Web で同じ canonical values を使う。
- feature anchor と actual feature lane の解決は radial track layout の一箇所を authority とする。
- leader contact と text anchor flip を renderer、solver、Web へ重複実装しない。

### 9.4 Review gate

PR review では次を reject 条件とする。

- `if placement == "radial"` が candidate、assembler、renderer、canvas に無秩序に散在する。
- radial mode だけの bbox formula が複数ファイルにある。
- label count/length による silent skip。
- collision を warning だけ出して残す。
- nondeterministic iteration/randomness。
- Web JavaScript で label geometry を再実装する。
- existing horizontal references の理由のない更新。

## 10. Implementation phases

### Phase 0: Characterization and contracts

実装前に現行 behavior を test で固定する。

1. `prepare_label_list()` の candidate text、embedded flag、outer/inner assignment。
2. horizontal label bbox、leader anchors、legend fallback、canvas expansion。
3. `NC_001879.gbk`、`MjeNMV.gbk`、HmmtDNA、`NC_001454.1.gbk` の現行 horizontal output invariants。
4. Linear rotated-label geometry の current corner/bounds/collision behavior。
5. config/session に placement field がない場合の current defaults。

完了条件:

- existing reference SVG に差分がない。
- geometry extraction 前後を比較できる characterization test がある。

### Phase 1: Shared geometry and typed placement option

1. rotated text geometry を `layout/text_geometry.py` へ抽出する。
2. Linear label placement を shared helper へ切り替える。
3. `LabelsCircularConfig` と placement normalization を追加する。
4. `config.toml` default、`modify_config_dict(circular_label_placement=...)`、API override を追加する。
5. invalid placement を CLI/API/config で同じ意味の validation error にする。

完了条件:

- Linear/Circular horizontal の全 test/reference が変更なしで通る。
- geometry owner が一箇所である。

### Phase 2: Candidate/layout boundary refactor

1. semantic candidate builder を `circular.py` から分離する。
2. embedded 判定を feature-layout projection 後の共通処理にする。
3. horizontal solver を strategy contract の adapter で包む。
4. common placement result と authoritative bounds を導入する。
5. renderer と canvas consumer を common result に切り替える。
6. 既存の未使用/重複 layout type が新 contract と競合しないか確認し、`CircularFeatureLabelLayout` は
   authority として採用するか削除するかを決め、二つの型定義を残さない。

完了条件:

- horizontal output の byte/semantic reference が維持される。
- candidate text measurement が feature ごとに一度だけである。
- assembly が mutable raw label dict の内部 heuristic に依存しない。

### Phase 3: Full radial solver

1. quadrant orientation と anchor flip を実装する。
2. outer cyclic packing、seam handling、radial tiers を実装する。
3. inner cyclic packing と center/track reservation を実装する。
4. long-label radius requirement と monotonic preflight を実装する。
5. text-text、leader-text、leader-leader、feature/track collision の final validator を実装する。
6. deterministic polishing を実装する。
7. sparse/dense、short/long、outer/inner の synthetic unit tests を追加する。

完了条件:

- rendering policy が solver へ渡した external label 数と rendered external label 数が一致する。
- final collision counters が全 category で 0。
- 同じ入力の repeated run が同じ geometry を返す。
- label-count/text-length cutoff が存在しない。

### Phase 4: Assembly, canvas, legend, multi-record

1. track/label preflight を SVG draw 前へ組み込む。
2. inner deficit に応じた radius re-resolution を実装する。
3. authoritative bounds で canvas を拡張する。
4. radial-aware legend obstacle resolution を実装する。
5. multi-record drawing metadata と packing を更新する。
6. custom track slots、resolve-overlaps lanes、definition reservation と統合する。

完了条件:

- label、leader、legend、隣接 record の clipping/overlap がない。
- top/bottom/side/corner legend と plot title の minimum gap が維持される。
- custom fixed geometry の矛盾が actionable error になる。

### Phase 5: Rendering and all output formats

1. external label text の rotation/anchor rendering を追加する。
2. two-line leader sibling order と group IDs を維持する。
3. interactive SVG metadata/label editor を確認する。
4. CairoSVG の PNG/PDF/EPS/PS conversion test を追加する。
5. Browser PNG/PDF/SVG export の targeted test を追加する。

完了条件:

- 全 format で orientation と bounds が一致する。
- Web label edit/reflow 後も radial placement が保持される。

### Phase 6: CLI, API, session, Web UI

1. Circular CLI parser/help と validation を追加する。
2. radial では `--labels both` による GC/skew の blanket suppression を解除し、horizontal の現行挙動は
   保持する。
3. API config override/request builder を接続する。
4. Python session import/export mapping を追加する。
5. Web state、select、run args、capability check を追加する。
6. canonical session round-trip と config save/load を追加する。
7. browser wheel を準備し、実 browser flow を確認する。

完了条件:

- CLI/API/Web が同じ config value と SVG geometryを生成する。
- old session/config は horizontal のまま開く。
- radial session は save/load 後も radial を保持する。

### Phase 7: Documentation and references

1. `docs/CLI_Reference.md` に Circular placement option を追加する。
2. Circular labels tutorial に horizontal/radial の使い分けと full example を追加する。
3. `NC_001879.gbk` と gene qualifier priority/color table を使った radial example を追加する。
4. `tests/reference_outputs/` に radial 専用 reference を新規追加する。
5. existing horizontal references は更新しない。差分が出た場合は regression として調査する。

## 11. Test plan

### 11.1 Shared geometry unit tests

- 0、90、180、270 度と任意角の corner/AABB。
- `start|middle|end` anchor。
- padding と edge-touch semantics。
- rotated rectangle SAT: separate、touch、overlap、containment。
- segment vs rectangle: edge touch、corner touch、interior crossing。
- Linear の既存 helper と shared helper の characterization equality。

### 11.2 Radial solver unit tests

- four quadrants の readable rotation と text extent direction。
- 12/6 時 tie-break stability。
- sparse labels が preferred angle に留まる。
- dense contiguous cluster が overlap なしで広がる。
- 0/360 度を跨ぐ cluster。
- 同一 midpoint の複数 feature。
- long outer label が canvas requirement を増やす。
- long inner label が axis-radius requirement を増やす。
- outer/inner mixed labels。
- `auto|embedded_only|external_only`。
- `out|both`。
- radial の `both` と GC/skew/depth/comparison/custom inner tracks の共存。
- `resolve_overlaps` の複数 feature lane。
- `separate_strands`。
- origin-spanning と multipart feature。
- Unicode、空白、記号を含む label override。
- no text-text、leader-text、leader-leader collision。
- cyclic order preservation。
- repeated-run determinism。
- 500 以上の synthetic labels でも skip せず、broad-phase candidate checks が全 pair scan を避ける。

wall-clock の fragile な秒数だけを performance test にしない。candidate pair 数、measurement call 数、
solver pass 数の上限を fixture ごとに検証し、必要に応じて slow benchmark を併用する。

### 11.3 Dataset integration tests

最低限、次を使う。

| Dataset | 目的 |
|---|---|
| `tests/test_inputs/NC_001879.gbk` | 約156 kb plastid、短い gene labels、高密度、Issue の代表例 |
| `tests/test_inputs/MjeNMV.gbk` | dense viral labels、long product labels、inner labels |
| HmmtDNA fixture | circular origin、tRNA、large font、resolve-overlaps clearance |
| `tests/test_inputs/NC_001454.1.gbk` | multipart/origin/inner regression |
| synthetic records | 同角度、極端な長文、極端な密度、明示的 track conflict |

各 dataset で `tuckin|middle|spreadout`、labels `out|both`、rendering `auto|external_only` を
必要な pairwise matrix に絞って実行する。ただし機能対応自体は全組合せとし、冗長な全直積を reference
SVG にしない。

### 11.4 Canvas and multi-record tests

- single record の全 label box が viewBox 内。
- top/bottom/left/right/corner legend と label bounds の gap。
- plot title/definition と inner labels の clearance。
- numeric/depth/conservation/tick tracks との clearance。
- custom absolute/relative circular track slots。
- mixed-length multi-record、各 scaling mode、row arrangement。
- record content bounds が rotated text を含む。
- adjacent record labels が column/row gap を侵食しない。

### 11.5 CLI/API/session tests

- Circular CLI default は `horizontal`。
- `--label_placement radial` parse と config forwarding。
- radial の `--labels both` が GC/skew request を保持し、horizontal は既存 suppression contract を保つ。
- invalid value rejection。
- Linear の同名 option semantics に回帰がない。
- `modify_config_dict(circular_label_placement="radial")`。
- typed config missing/new field。
- single/multi-record public API config override。
- canonical request codec/session round-trip。
- legacy session import/export。
- old config/session missing field -> horizontal。

### 11.6 Web tests

- Circular label placement select の表示条件と default。
- state -> args -> Python wrapper の `--label_placement radial`。
- wheel capability error。
- config save/load、session save/load、undo/history。
- radial SVG の rotated external text と二本の leader line。
- feature label editor で label text を変更し、reflow 後も radial。
- SVG/PNG/PDF export が clip しない。
- Node Playwright がなければ Python Playwright で同じ targeted flow を実行する。

### 11.7 Reference output policy

- `circular_radial_labels` の新規 reference を追加する。
- reference は short gene labels だけでなく、少なくとも一つの long override、outer/inner、0/360 seam を含む。
- intentional radial geometry は新規 reference で review する。
- existing horizontal reference の更新を radial feature の実装に混ぜない。

### 11.8 Suggested commands

```bash
python -m pytest tests/test_circular_radial_labels.py -q
python -m pytest tests/test_circular_label_placement.py tests/test_circular_radial_layout.py -q
python -m pytest tests/test_circular_multi_canvas.py -q
python -m pytest tests/test_api_library_usage.py tests/test_api_request_render.py tests/test_api_session.py -q
python -m pytest tests/test_session_io.py tests/test_session_request_codec.py -q
node tests/web/session-request.test.mjs
python -m pytest tests/test_web_packaging.py -k "label or session or index" -q
python -m pytest tests/test_output_comparison.py::TestOutputComparison -v
python -m pytest tests/ -q -m "not slow"
ruff check gbdraw/
python -m build
```

CairoSVG と browser verification は環境 availability を両経路で確認し、sandbox launch error の場合は
同じ check を必要な権限で再実行する。

## 12. Acceptance criteria

すべて満たしたとき Issue #289 を完了とする。

1. Circular external labels に `horizontal|radial` を選べる。
2. default horizontal の既存出力と既存 session が変わらない。
3. radial が outer/inner の両方を描ける。
4. radial が短い/長い label を省略、短縮、wrap、font shrink せず描ける。
5. `auto|embedded_only|external_only` の semantics が維持される。
6. `out|both`、全 preset、custom tracks、resolve-overlaps、separate-strands で動作する。
7. rendering policy が radial solver へ渡した external label count と rendered external label count が一致する。
8. final text-text、leader-text、leader-leader、label-feature/track collisions が 0。
9. genome と label の cyclic order が維持される。
10. left-half labels が逆さにならず、axis tie-break が安定する。
11. inner labels の必要空間が preflight で track/canvas geometry に反映され、radial の `labels=both` が
    GC/skew/depth/comparison/custom inner tracks を黙って suppress しない。
12. 全 label/leader が canvas/viewBox 内に収まる。
13. legend、title、definition、multi-record peers と重ならない。
14. SVG、interactive SVG、PNG、PDF、EPS、PS で orientation と clipping が一致する。
15. CLI、Python API、session、Web UI が同じ canonical setting を使う。
16. Web feature-label edit/reflow と export が radial で動作する。
17. radial solver に label-count/text-length cutoff、randomness、silent-overlap fallback がない。
18. rotated geometry、candidate extraction、bounds computation が重複していない。
19. targeted tests、fast suite、reference comparison、ruff、build、browser smoke が通る。
20. 実行した test、未実行 test と理由、intentional new reference を PR に記録する。

## 13. Risks and mitigations

| Risk | Impact | Mitigation |
|---|---|---|
| candidate extraction refactor が horizontal を変える | 大きな reference diff | Phase 0 characterization、horizontal adapter、既存 reference 不変 gate |
| rotated bbox と SVG baseline がずれる | clipping/見かけ上の overlap | shared anchor/baseline geometry、export/browser getBBox spot check |
| dense labels で計算量が増える | CLI/Pyodide freeze | AABB broad phase、contiguous cluster、measurement cache、pair/pass counters |
| inner label と track が循環依存する | layout 非収束 | semantic prepass、radius 単調増加、exact deficit、stable convergence assertion |
| explicit radius が inner requirement と矛盾する | silent corruption | slot/available/required を含む ValidationError |
| legend reflow 後に geometry が stale | leader/text mismatch | immutable result replacement、bounds/leader を同時再計算 |
| multi-record estimator が rotate を無視する | record 間 overlap | authoritative drawing metadata を estimator と merge |
| CairoSVG/browser で rotation 差 | format ごとの崩れ | plain text rotate、全 export targeted tests |
| Web label editor が sibling order に依存 | feature association failure | existing two-line + text order を維持、browser test |
| config/session naming が Linear と衝突 | wrong option restore | CLI 名は mode-local、canonical key は `circular_label_placement` |

## 14. Explicit non-goals

この Issue は radial arrangement を全 existing feature-label contract へ実装する。次は別機能であり、
radial 対応を短縮するための除外ではない。

- multiline/rich-text feature label format の新設。
- user が label を一つずつ drag して座標を永続化する manual editor。
- radial text を文字単位に曲げる新しい curved-text style。
- Linear mode の新しい radial 相当 layout。
- label filtering/qualifier semantics の変更。
- default placement を horizontal から radial へ変更すること。

## 15. Definition of done

- Acceptance criteria をすべて満たす。
- facade、candidate builder、horizontal/radial strategy、shared geometry、renderer、assembly の責務が分離されている。
- radial mode の full dataset matrix に outer/inner、short/long、dense、all rendering policies が含まれる。
- selected label を隠す cutoff または unresolved collision がない。
- horizontal reference output が不変である。
- radial の新規 reference SVG を目視 review し、Issue 添付図の読み方向、密度、leader behavior を満たす。
- CLI/API/session/Web docs と user-facing help が同じ terminology を使う。
- PR description に SOLID/KISS/DRY review gate の確認結果を記録する。
