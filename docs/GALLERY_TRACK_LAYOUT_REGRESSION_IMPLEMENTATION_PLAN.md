# Gallery Track Layout Regression Implementation Plan

- 作成日: 2026-07-28
- 調査基準: `73bf033` (`Remove branch-only session compatibility layers`)
- 状態: **次セッション向け実装計画。本セッションでは実装・生成物更新を行わない**
- 対象:
  - `HmmtDNA_ATskew`: `middle` が canonical session replay で `tuckin` 相当になる回帰
  - `hepatoplasmataceae_collinear` / `hepatoplasmataceae_orthogroup`:
    Feature と GC content の間隔が過大になる回帰
- 関連設計:
  - [`LINEAR_TRACK_OCCUPANCY_LAYOUT_IMPLEMENTATION_PLAN.md`](LINEAR_TRACK_OCCUPANCY_LAYOUT_IMPLEMENTATION_PLAN.md)
  - [`LOSAT_CACHE_IDENTITY_AND_LINEAR_SPACING_REGRESSION_IMPLEMENTATION_PLAN.md`](LOSAT_CACHE_IDENTITY_AND_LINEAR_SPACING_REGRESSION_IMPLEMENTATION_PLAN.md)

> **実装時の受入条件更新 (2026-07-28):** 画面確認に基づき、Linear の最終 target は
> Feature→GC content、GC content→GC skew、GC skew→comparison の3境界を共通の
> `track_spacing = 0 px` にする。以下の 8 px は回帰調査時の旧 contract / baseline を示し、
> 最終受入値ではない。

## 0. 結論

今回の二つの表示不具合は、同じ Gallery 更新で表面化したが、根因は独立している。

1. **Workstream A — HmmtDNA circular preset**
   - UI state は `track_type=middle`、Feature は `side=overlay`,
     `lane_direction=split` で正しい。
   - canonical request を作る際に、`split` が「middle preset の既定値」として省略され、
     `features:features` だけが保存された。
   - Python の custom slot parser は、意味が省略された Feature を正しく `inside` と解釈する。
     そのため canonical session replay だけが `tuckin` 相当になる。
   - live Generate と canonical session writer が異なる serializer option を使っていることが
     直接原因である。
2. **Workstream B — Hepatoplasmataceae linear spacing**
   - Default と明示 Linear Track Slots は同じ意味の設定なのに、Default にだけ
     legacy preferred origin の上書きが入る。
   - per-record packer は、その axis 基準 preferred gap を、すでに Feature footprint を含む
     occupied edge の外側から再加算する。
   - 結果として Feature の下側 extent が実質的に二重計上され、宣言上 8 px の間隔が
     reserve band 間で約 23.6 px になる。

二つの修正が共有するのは、Gallery の再生成と geometry contract の最終検証だけとする。
Circular preset の修正で Linear solver を変更せず、Linear spacing の修正で session preset の
解釈を変更しない。

採用しない対症療法は次のとおり。

- Gallery example ID、accession、比較方式に対する条件分岐を追加する。
- HmmtDNA の半径や Hepatoplasmataceae の Y 座標を固定値で補正する。
- `vertical_padding` を小さくして過大 gap を隠す。
- `track_type=middle` を Python の全 custom slot に暗黙適用し、既存の bare spec の意味を変える。
- session、source SVG、Gallery SVG、thumbnail を手作業で編集する。
- 意図を確認せず reference SVG を一括更新する。
- Track Slot editor の human-readable な compact 表示まで transport format に合わせて変更する。

## 1. 再現結果と固定する baseline

### 1.1 HmmtDNA_ATskew

現行 session の同一 document 内に、相反する二つの表現が保存されている。

| 観測箇所 | 現行値 | 期待値 |
|---|---|---|
| `config.form.track_type` | `middle` | `middle` |
| `config.adv.circular_track_slots[0]` | `side=overlay`, `lane_direction=split` | 同左 |
| canonical `circularTrackAxisIndex` | `0` | `0` |
| canonical Feature spec | `features:features` | Feature が `split` であることを明示 |

SVG geometry も canonical request の意味欠落と一致する。

| Geometry | 現行 | 正常な middle |
|---|---:|---:|
| Axis radius | 390.0 px | 390.0 px |
| Feature radial band | 349.05–386.10 px | 371.475–408.525 px |
| Axis と Feature の関係 | Feature が全て内側 | `inner < axis < outer` |

現行 session を CLI の `--session` で再描画しても誤配置を再現するため、enrichment や
browser viewer のみの問題ではない。`lane_direction=split` を明示した同値 request は正常な
middle geometry を生成する。

履歴上は次の順に表面化した。

- `4fcc37c` (2026-07-15): canonical session writer が live Generate と異なる
  serialization path を持ち、潜在的不整合が入った。
- `520ba22` (2026-07-20): schema-3 Gallery session の promotion / replay で見える回帰になった。
- `aa555a0` (2026-07-24): tracked source SVG も誤った内側配置へ更新された。

### 1.2 Hepatoplasmataceae

`collinear` と `orthogroup` は比較データが異なるが、各 record の track Y 座標は同一である。
したがって比較方式や LOSAT result は根因ではない。

先頭 record の代表値は次のとおり。

| Geometry | 旧 Gallery | 現行 Gallery |
|---|---:|---:|
| Feature / Axis group Y | 67.0 px | 67.0 px |
| GC content group Y | 92.0 px | 113.5 px |
| center delta | 25.0 px | 46.5 px |
| Feature→GC の実描画 gap | 約 3 px | 約 24.5 px |
| Feature→GC の reserve gap | 未記録 | 約 23.6 px |
| 宣言された既定 spacing | — | 8.0 px |

同値な明示 Linear Track Slots では、調査時の geometry metadata 上で GC origin は約 30.9 px、
Feature→GC reserve gap は 8.0 px になる。Default path だけが約 46.5 px へ膨らむ。
これらの小数値は原因確認の参考値であり、受入条件の固定値ではない。

修正後の target は旧版の center delta 25 px ではない。旧版は Feature と GC の実描画間隔が
約 3 px しかなく、現行の 8 px spacing contract を満たさない。受入条件は次の関係式とする。

```text
gc.reserve.top - feature_composite.reserve.bottom
  ≈ configured spacing within the declared geometry tolerance
```

現在の Feature / repeat-region underlay / stroke を含む footprint では、GC center は概ね
30.5–31 px 付近になる見込みだが、この値自体を test の正本にしない。

履歴上は次の順に表面化した。

- `d94d14b` (2026-07-20): per-record occupancy planner と default preferred origin を導入した。
- `38ca401` (2026-07-20): sparse-depth compaction を追加したが、同じ preferred-gap 構造を維持した。
- `520ba22` (2026-07-20): Gallery render に過大 gap が現れた。
- `aa555a0` (2026-07-24): tracked source SVG / thumbnail に反映された。
- `16762f3` (2026-07-26): Hepatoplasmataceae assets の再生成でも回帰が維持された。

## 2. 根因

### 2.1 Circular: preset-relative compression を canonical transport に使っている

現行 Web path は概ね次の連鎖を持つ。

```text
UI slot
  side=overlay, lane_direction=split, preset=middle
        |
        v
normalizeCircularTrackSlot()
  「preset の既定形」と判定して side/lane を省略
        |
        v
buildCircularTrackSlotSpec()
  canonical writer は forceSplitLane=false
        |
        v
renderRequest: "features:features", axisIndex=0
        |
        v
Python circular_track_slots_with_axis_side()
  slot index 0 は Axis の内側 -> lane_direction=inside
```

該当する主な箇所は次のとおり。

- `gbdraw/web/js/app/circular-track-slots.js`
  - preset default shape の検出
  - slot normalization
  - `buildCircularTrackSlotSpec()`
- `gbdraw/web/js/services/session-request.js`
  - canonical `buildTracks()`
- `gbdraw/web/js/app/run-analysis.js`
  - live CLI args の生成
- `gbdraw/tracks/circular.py`
  - bare Feature spec と Axis boundary の context-free 解決

live Generate は `forceSplitLane: true` を渡すため、
`features:features@lane_direction=split` を生成して正常に動く。一方、canonical writer は同じ
UI state を `features:features` に縮約する。これは一つの意味に二つの serializer policy がある
DRY 違反である。

`track_type` は preset から default slots を作るための入力であり、custom slots の意味を後から
補完する authority ではない。canonical custom slot は、それ単独と Axis boundary から
同じ配置へ戻せる必要がある。

### 2.2 Linear: absolute preferred origin を relative gap として再適用している

Default Linear path では、`resolve_linear_track_layout()` の後にだけ次の処理が入る。

- `gbdraw/diagrams/linear/assemble.py`
  - `_layout_with_preferred_origins()`
  - `_default_preferred_origins()`
  - 二つの default-only call site

代表的な Middle case では概ね次の値になる。

```text
base = cds_padding 17 + vertical_padding 8 = 25
nominal GC origin = 18
axis_clearance = 8
overridden preferred GC origin = 25 + 18 - 8 = 35
```

その後、`gbdraw/diagrams/linear/track_slots.py::_pack_side_slots()` は、
preferred band と Axis の間を約 24 px の `preferred_gap` に変換する。その gap を、
すでに Feature / underlay reserve band の下端約 12.5 px から再び足す。

```text
occupied edge already contains Feature extent
  + axis-relative preferred gap
  + GC local top extent
  = GC origin 約 46.5 px
```

これは [`LINEAR_TRACK_OCCUPANCY_LAYOUT_IMPLEMENTATION_PLAN.md`](LINEAR_TRACK_OCCUPANCY_LAYOUT_IMPLEMENTATION_PLAN.md)
で定めた次の式とも一致しない。

```python
resolved_origin = max(
    absolute_preferred_origin,
    occupied_bottom + spacing - local_reserve_top,
)
```

Default-only override を外した同値 explicit path は 8 px gap を生成するため、比較、record 数、
palette、repeat-region underlay は主原因ではない。underlay は Feature composite reserve band を
約 0.5 px 広げるが、その突出分だけを反映するのが正しい。

## 3. 修正後の contract

### 3.1 Circular semantic serialization

1. preset expansion 後の Feature placement は、canonical / CLI transport で省略しない。
2. `middle` Feature は少なくとも `lane_direction=split` を必ず保持する。
3. live Generate、Save Session、Gallery promotion が同じ pure serializer を使う。
4. canonical custom slot の解釈は `track_type` に依存しない。
5. Python の既存 contract を維持し、bare `features:features` を一律 middle へ読み替えない。
6. explicit `inside` / `outside` / `split` は preset にかかわらず round-trip 後も同じ意味を持つ。
7. current schema-3 session の `renderRequest` は通常どおり authority とし、load 時の
   implicit re-projection は追加しない。

### 3.2 Linear occupancy and spacing

1. Default options と意味的に同値な explicit slots は同じ normalized intent と geometry を持つ。
2. `spacing_after_px` は隣接 reserve-band edge 間に一度だけ適用する。
3. `preferred origin` を残す場合、それは Axis からの絶対座標であり、gap ではない。
4. Below slot の origin は次の式で解く。

   ```python
   collision_origin = occupied_bottom + incoming_spacing - local_reserve_top
   resolved_origin = max(preferred_origin, collision_origin)
   ```

5. Above slot は同じ helper の鏡像として解く。

   ```python
   collision_origin = occupied_top - incoming_spacing - local_reserve_bottom
   resolved_origin = min(preferred_origin, collision_origin)
   ```

6. Feature と overlay underlay は composite footprint とし、underlay が Feature band 内なら
   numeric track を動かさない。突出時は突出量だけ動かす。
7. missing Depth の compaction contract は維持し、今回の変更に便乗して再定義しない。
8. comparison kind と record-row spacing は record-local Feature→GC geometry を変えない。
9. renderer は引き続き resolved origin のみを使い、独自の Y 補正を持たない。
10. Default の `track_spacing` は 0 px とし、Feature→GC、GC→skew、
    skew→comparison exclusion boundary を同じ接触条件にする。

### 3.3 SOLID / DRY / KISS の適用

| 原則 | 適用方針 |
|---|---|
| Single Responsibility | Web serializer は意味保存、Python parser は明示 spec の解釈、Linear packer は band constraint の解決だけを担当する。 |
| Open/Closed | Gallery ID 分岐ではなく、全 preset / 全 record に適用できる semantic・geometry contract を追加する。 |
| Liskov / equivalence | Default と同値 explicit slots、live と session replay を代替可能にする。 |
| Interface Segregation | preset、canonical transport、paint geometry の入力を混ぜず、必要な情報だけを各層へ渡す。 |
| Dependency Inversion | Gallery validation は SVG の固定 path ではなく、slot semantics と geometry metadata / band relation に依存する。 |
| DRY | live/session の serializer と Default/Custom の Linear layout path を各一つにする。 |
| KISS | 半径・Y 座標の補正表を作らず、明示 lane と `max` / `min` の constraint で解く。 |

## 4. 実装計画

### Phase 0 — Preflight と failing characterization

実装前に、現在の誤動作を小さい test で固定する。

1. 作業ツリーには本計画外の変更が存在し得るため、`git status` と
   `git diff --ignore-space-at-eol` で semantic diff を確認する。
2. 本計画の対象ファイルと既存変更が重なる場合は、既存変更を上書きせず scope を分離する。
3. HmmtDNA の巨大 SVG 全文 snapshot ではなく、次を抽出する helper を test 側に用意する。
   - Axis radius
   - Feature inner / outer radius
   - canonical Feature lane
4. Linear geometry metadata から次を抽出する共通 helper を用意する。
   - slot `origin_y`
   - `paint_band`
   - `reserve_band`
   - adjacent reserve gap
5. 最初に以下の failing assertions を追加する。
   - middle canonical request が explicit `split` を保持する。
   - live serializer と canonical serializer の semantic result が一致する。
   - default-only preferred-origin rewrite が Default と同値 explicit slots を不一致にする。
   - synthetic seed band に対し、`preferred_gap` 再加算が absolute preferred-origin contract と
     不一致になる。
   - Feature composite→GC reserve gap が `track_spacing` と許容誤差内で一致する。
6. no-feature、custom preferred origin、sparse Depth を修正前に characterization し、
   二つの Linear 変更を別々に適用しても意図した挙動を失わないことを確認する。

完了条件:

- 修正前に、HmmtDNA semantic test と Default Linear gap test が意図した理由で失敗する。
- unrelated cache、comparison、palette assertion は失敗しない。

### Phase 1 — Circular serializer を lossless に統合する

1. `gbdraw/web/js/app/circular-track-slots.js` で、effective Feature placement を一度だけ解決する
   pure helper を作る。
2. transport 用 slot serializer は、preset default と同じであっても Feature lane の意味を保持する。
   `inside` / `outside` / `split` の全てを明示できる形を優先する。
3. `forceSplitLane` のような呼び出し元ごとの policy 差を削除する。互換上すぐ削除できない場合も、
   public call sites ではなく共通 transport wrapper 内に閉じ込め、同一入力から同一 spec を返す。
4. `gbdraw/web/js/services/session-request.js` と
   `gbdraw/web/js/app/run-analysis.js` を同じ serializer へ切り替える。
5. Track Slot editor の compact formatter は UI 表示責務として維持してよい。ただし、
   effective placement を求める pure helper は transport serializer と共有する。
6. Python の `gbdraw/tracks/circular.py` には HmmtDNA や `track_type` fallback を追加しない。
   explicit split の既存解釈だけで正しい radial layout になることを test する。

完了条件:

- middle Feature は canonical / live とも `split`。
- tuckin は `inside`、spreadout は `outside`。
- project → canonical → project と save → CLI replay が冪等。
- explicit placement と Axis boundary の conflict validation は弱めない。

### Phase 2 — Linear Default/Custom layout と band solver を統合する

1. `gbdraw/diagrams/linear/assemble.py` の default-only preferred-origin rewrite と二つの
   call site を削除し、Default/explicit equivalence test だけを先に再実行する。
2. Default options も explicit slots も、次の一経路だけを通す。

   ```text
   slot inputs
     -> normalize
     -> resolve_linear_track_layout
     -> measure per-record footprints
     -> resolve_linear_record_vertical_plan
     -> render resolved origins
   ```

3. 1 の差分と test 結果を確認してから、`gbdraw/diagrams/linear/track_slots.py::_pack_side_slots()`
   から `preferred_inner_edge` / `preferred_gap` の再加算を除く。
4. Above / Below 共通の小さい pure origin resolver を作り、3.2 の `max` / `min` 制約を
   一箇所に実装する。
5. `initial_spacing` は Feature composite の `spacing_after_px`、以後は直前 slot の
   `spacing_after_px` として一度だけ消費する。
6. `compact_missing` は既存の sparse Depth test で保護し、欠損 slot が後続 track の
   不要な空間を残さない contract を維持する。
7. Default の旧固定中心を残す必要が test で判明した場合も、assembly caller 分岐を復活させない。
   normalized slot intent の明示 spacer / preferred origin として `track_slots.py` 内で一度だけ表現する。

完了条件:

- Middle + separate strands の Feature→GC reserve gap は 0 px。
- Default と同値 explicit slots が同じ geometry。
- Feature lane が増えた record だけ、必要量だけ numeric tracks が外側へ移動する。
- Above は Below の鏡像、Below は Feature と numeric が同じ側でも非衝突。
- GC→skew、Depth→GC、annotation→numeric の spacing は二重加算されない。

### Phase 3 — Session repair と Gallery 再生成

コードと focused tests が通った後にだけ生成物を更新する。

#### 3.1 HmmtDNA current session の再canonicalize

現行 HmmtDNA session はすでに schema 3 である。`refresh_gallery_sessions.py` は current
schema-3 `renderRequest` を authority として promoter を迂回するため、refresh だけでは
欠落した `split` を修復できない。

1. refresh 内部の `_promote_gallery_session()` ではなく、standalone promoter を明示的に実行し、
   既存 `config` から一時ファイルへ再canonicalizeする。

   ```bash
   gallery_tmp_dir="$(mktemp -d)"
   node tools/promote_gallery_session.mjs \
     gbdraw/web/gallery/sessions/HmmtDNA_ATskew.gbdraw-session.json \
     "$gallery_tmp_dir/HmmtDNA_ATskew.gbdraw-session.json"
   ```

2. 一時ファイルで次を検証してから tracked session と置換する。
   - session version / request schema は不変。
   - Feature spec は explicit split。
   - AT skew slot、palette、qualifier priority、resources、feature catalog は欠落しない。
3. semantic diff と focused migration test が通った後だけ、一時出力を tracked session に採用する。
   standalone promoter の出力であることを記録し、JSON を手編集しない。
4. current session load 時の一般的な auto-reprojection は追加しない。この一回の artifact migration と
   将来の writer fix を分離する。

#### 3.2 対象 session と assets の再生成

対象を明示して refresh する。

```bash
python tools/refresh_gallery_sessions.py \
  --session HmmtDNA_ATskew \
  --session hepatoplasmataceae_collinear \
  --session hepatoplasmataceae_orthogroup
```

この command は session を対象指定できるが、asset preparation は Gallery inventory 全体を処理し、
`sources/`、`examples/`、`thumbnails/`、`examples.json` を広く更新し得る。全差分を確認し、
今回と無関係な drift を混ぜない。

次を直接編集しない。

- `gbdraw/web/gallery/sessions/*.gbdraw-session.json[.gz]`
- `gbdraw/web/gallery/sources/*.svg`
- `gbdraw/web/gallery/examples/*.svg`
- `gbdraw/web/gallery/thumbnails/*.webp`
- `gbdraw/web/gallery/examples.json`

完了条件:

- HmmtDNA session result、source SVG、enriched Gallery SVG が同じ middle geometry。
- Hepatoplasmataceae の全 5 records で Feature→GC gap が contract を満たす。
- collinear / orthogroup の record-local track geometry は同じ。
- comparison count、identity、link metadata は更新前と同じ。
- enriched SVG の search、tooltip、zoom、drag metadata が維持される。

### Phase 4 — Regression gates

#### 4.1 Web unit / migration

- `tests/web/session-request.test.mjs`
  - middle / tuckin / spreadout の Feature lane serialization
  - canonical round-trip
  - live / session serializer equivalence
  - explicit lane が preset に上書きされないこと
- `tests/web/gallery-session-migration.test.mjs`
  - HmmtDNA promotion が `split` を保持すること
  - 再promotion が冪等であること
- `tests/test_web_packaging.py`
  - call-site 固有の `forceSplitLane` を期待する構造 test があれば、semantic assertion へ置き換える。

#### 4.2 Python geometry

- `tests/test_circular_track_slots.py`
  - explicit split + Axis index 0 が overlay/on-axis に解決されること
  - explicit inside / outside と conflict error が維持されること
- `tests/test_linear_vertical_layout.py`
  - Default / explicit equivalence
  - Middle / Above / Below
  - stranded / non-stranded
  - Feature lane 1 / multiple lanes
  - repeat-region underlay の有無と突出
  - GC-only、skew-only、Depth→GC→skew、GC percent mode
- `tests/test_linear_track_layout.py` または `tests/test_linear_track_slots.py`
  - pure `max` / `min` origin constraint
  - spacing once-only
  - sparse Depth compaction

assertion は `y == 67` や `radius == 390` のような example 固有値ではなく、次を使う。
Linear の float / stroke geometry は `pytest.approx` 相当の明示 tolerance で比較する。

```text
circular middle: feature_inner < axis < feature_outer
linear below: next.reserve.top - previous.reserve.bottom ≈ spacing
linear above: previous.reserve.top - next.reserve.bottom ≈ spacing
```

#### 4.3 Gallery semantic / browser

- `tests/test_gallery_session_semantics.py`
  - HmmtDNA canonical Feature が split。
  - HmmtDNA session result / source の Feature band が Axis を跨ぐ。
  - Hepatoplasmataceae の全 record で Feature→GC reserve gap が宣言値と一致。
  - two comparison modes の record-local track geometry が一致。
- `tests/test_refresh_gallery_sessions.py`
  - staged artifact validation が semantic geometry regression を拒否する。
  - transaction / schema / artifact sync の既存 contract を維持する。
- focused Playwright または Python Playwright
  - HmmtDNA と Hepatoplasmataceae を Gallery から開く。
  - SVG metadata / `getBBox()` を使う relational check を一つずつ行う。
  - viewer interactivity が壊れていないことを確認する。

巨大 interactive SVG の全文 snapshot は追加しない。共有 geometry helper と metadata を優先する。

### Phase 5 — Reference、full validation、目視

1. まず tracked references を更新せず比較する。
2. 意図した geometry 差だけであることを確認する。
3. reference 更新が必要な場合だけ、repository 指定の generator を使う。
4. 更新された SVG は、Feature/GC/Skew、labels、ruler、comparison、viewBox を目視 review する。
5. 最後に対象 Gallery の原寸表示と thumbnail を確認する。

## 5. 主な変更対象

| Workstream | 変更候補 | 役割 |
|---|---|---|
| Circular | `gbdraw/web/js/app/circular-track-slots.js` | effective placement と transport serializer の単一 owner |
| Circular | `gbdraw/web/js/services/session-request.js` | canonical request が共通 serializer を使用 |
| Circular | `gbdraw/web/js/app/run-analysis.js` | live path 固有 option を除去 |
| Circular tests | `tests/web/session-request.test.mjs` | semantic round-trip |
| Circular tests | `tests/web/gallery-session-migration.test.mjs` | Gallery promotion |
| Circular tests | `tests/test_circular_track_slots.py` | Python explicit lane contract |
| Linear | `gbdraw/diagrams/linear/assemble.py` | default-only preferred-origin owner を削除 |
| Linear | `gbdraw/diagrams/linear/track_slots.py` | band constraint による唯一の origin solver |
| Linear tests | `tests/test_linear_vertical_layout.py` | measured per-record geometry |
| Linear tests | `tests/test_linear_track_layout.py` / `tests/test_linear_track_slots.py` | pure packing contract |
| Gallery tests | `tests/test_gallery_session_semantics.py` | tracked session/source の意味と geometry |
| Gallery pipeline | `tools/refresh_gallery_sessions.py` | 必要な場合のみ generic staged validation を追加 |
| Generated | 対象 sessions / sources / examples / thumbnails / manifest | generator 経由でのみ更新 |

原則として変更不要な箇所:

- `gbdraw/tracks/circular.py`
  - explicit custom spec の context-free 解釈は正しい。
- `gbdraw/session_request_codec.py`
  - string を lossless に運ぶ責務であり、preset fallback を追加しない。
- LOSAT cache / comparison converter
  - Hepatoplasmataceae の record-local gap とは無関係。

## 6. Test matrix

| Case | Expected circular placement | Expected Linear adjacency |
|---|---|---|
| middle default/custom | Feature が Axis を跨ぐ | Feature composite→first numeric = declared spacing |
| tuckin | Feature outer ≤ Axis | 同上 |
| spreadout | Feature inner ≥ Axis | 同上 |
| explicit inside/outside/split | preset に依存せず指定どおり | — |
| stranded / non-stranded | — | measured Feature band に応じて必要量だけ移動 |
| overlap lanes 1 / N | — | gap は一定、origin だけ lane 増分移動 |
| repeat underlay in-band | — | underlay なしと同じ origin |
| repeat underlay protruding | — | 突出量だけ origin が移動 |
| missing Depth | — | 空 slot を compact、後続 spacing は保持 |
| no comparison / pairwise / orthogroup / collinear | — | record-local geometry は同一 |
| save/load/session replay | live と同じ | direct Generate と同じ |

## 7. 検証 command

Focused checks:

```bash
node tests/web/session-request.test.mjs
node tests/web/gallery-session-migration.test.mjs

python -m pytest \
  tests/test_circular_track_slots.py \
  tests/test_linear_track_slots.py \
  tests/test_linear_vertical_layout.py \
  tests/test_linear_track_layout.py \
  tests/test_gallery_session_semantics.py \
  tests/test_refresh_gallery_sessions.py \
  -q
```

Reference / broader checks:

```bash
python -m pytest tests/test_output_comparison.py::TestOutputComparison -v
python -m pytest tests/ -v -m "not slow"
ruff check gbdraw/
git diff --check
```

Python の tool / test file を変更した場合は、上記に加えて変更した `.py` file だけを
`ruff check` の対象にする。無関係な `tools/` / `tests/` 全体の既存 lint を今回の gate に混ぜない。

意図した reference geometry 更新が必要な場合のみ:

```bash
python -m pytest \
  tests/test_output_comparison.py::TestGenerateReferences \
  --update-reference-outputs \
  -v

python -m pytest tests/test_output_comparison.py::TestOutputComparison -v
```

browser verification が必要な場合は、root の `package.json` の有無だけで可否を判断しない。

```bash
command -v playwright && playwright --version
python -c "from playwright.sync_api import sync_playwright; print('python playwright ok')"
node -e "console.log(require.resolve('@playwright/test'))"
```

Chromium が sandbox error で起動できない場合は、同じ local browser check を必要な権限で再実行する。

## 8. Definition of Done

### HmmtDNA

- session の `track_type=middle` と canonical Feature placement が矛盾しない。
- live Generate、Save→Load→Generate、CLI `--session` が同じ slot semantics と radial geometry を持つ。
- Feature inner radius `<` Axis radius `<` Feature outer radius。
- GC content、GC skew、AT skew、ticks の順序と legend は維持される。
- HmmtDNA 固有の補正や Python parser heuristic がない。

### Hepatoplasmataceae

- 全 record の Feature composite reserve bottom→GC reserve top が既定 0 px。
- GC→skew と skew→comparison exclusion boundary も 0 px で接し、clipping しない。
- Default と同値 explicit slots の geometry が一致する。
- collinear / orthogroup の record-local geometry が一致する。
- comparison ribbon / link count / identity / metadata は変わらない。
- old center delta 25 px ではなく、band relation を受入基準にしている。

### Repository / Gallery

- focused、reference comparison、fast test、lint が成功する。
- session、source、enriched SVG、thumbnail、manifest が generator 出力として同期する。
- generator が変更した全 asset diff を review 済み。
- unrelated user changes、cache identity、comparison behavior、palette を巻き込まない。
- 実装、test、意図した generated artifacts を一つの reviewed commit としてまとめる。

## 9. Risk と mitigation

| Risk | Mitigation |
|---|---|
| current HmmtDNA schema-3 session が refresh だけでは自己修復しない | serializer 修正後に promoter で一度再canonicalizeし、内容を検証してから refresh |
| bare Feature spec の既存利用者を壊す | Python の bare-spec contract は変更せず、writer と tracked artifact を修復 |
| Linear gap を詰めすぎて feature / underlay / stroke と衝突する | paint ではなく composite reserve band と宣言 spacing で検証 |
| preferred-origin 削除で no-feature / sparse-depth layout が変わる | GC-only、Depth missing、Above/Below characterization を修正前に追加 |
| reference が既に回帰後の geometry を正としている | relational invariant を先に追加し、reference diff を人が review |
| Gallery helper が対象外 assets も再生成する | target session 指定後も全 asset diff を列挙し、無関係 drift を除外 |
| 大型 session / SVG test が遅く brittle になる | pure unit + metadata helper を主とし、tracked Gallery test は最小限 |
| dirty worktree の既存変更を上書きする | preflight で semantic diff を確認し、対象外変更を保存・維持 |

## 10. 次セッションの実行順

- [ ] Preflight で既存変更と対象ファイルの重なりを確認する。
- [ ] HmmtDNA semantic round-trip と Linear reserve-gap の failing tests を追加する。
- [ ] Circular transport serializer を一つに統合する。
- [ ] Default-only Linear preferred-origin rewrite を削除する。
- [ ] Linear packer を absolute preferred origin + collision constraint に単純化する。
- [ ] focused Web / Python tests を通す。
- [ ] HmmtDNA current session を promoter で一度だけ再canonicalizeする。
- [ ] 3 Gallery sessions と assets を generator で再生成する。
- [ ] generated diff、reference diff、comparison metadata を review する。
- [ ] fast suite、lint、browser check、原寸 SVG / thumbnail 目視を完了する。
- [ ] 一つの commit にまとめ、commit title / summary を handoff する。

推奨 commit title:

```text
Fix gallery track layout regressions
```

推奨 commit summary:

```text
- Preserve explicit circular feature placement across live and canonical serialization.
- Unify default and custom linear track packing around measured reserve-band spacing.
- Regenerate and validate the affected HmmtDNA and Hepatoplasmataceae gallery assets.
```
