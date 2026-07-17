# Issue #239 Linear Multi-Record Layout Implementation Plan

- 作成日: 2026-07-17
- 対象バージョン: `0.14.0b0` 以降
- 設計スナップショット: branch `main`, HEAD `0e4730c` の作業ツリー
- 対象 Issue: [#239 Multiple records on a single line; multi-multi comparison](https://github.com/satoshikawato/gbdraw/issues/239)
- 状態: 設計方針合意済み、未実装
- 目的: Linear mode で複数 record を同じ row に配置し、隣接 row 間の任意の record pair を共通スケールで比較できるようにする

## 1. 結論

Issue #239 は、既存の「1 record = 1 row」という Linear layout を、
「1 row = 1 個以上の record」へ拡張する機能として実装する。

新しい layout の基本形は次のとおりとする。

```text
Row 1   [record A]   [record B]
             \          \
              \          \
Row 2   [record C]   [record D]   [record E]
                  \          \
Row 3             [record F]   [record G]
```

初回実装では次を必須とする。

1. 同じ row に複数の chromosome、contig、plasmid、region を横並びにできる。
2. 全 record を一つの共通 bp/px スケールで描く。
3. record ごとに 0 から始まる小さな ruler を既存 ruler option で表示できる。
4. 隣接 row 間で、明示された任意の record pair に BLAST / protein comparison を描ける。
5. 既存の feature、label、GC、GC skew、depth、annotation、custom track slot を record ごとに維持する。
6. `--records_table` の既存 `row` / `column` を Linear mode でも有効にする。
7. CLI、Python API、Web、canonical request、session、interactive SVG が同じ endpoint contract を保持する。
8. 従来の「1 row に 1 record」と隣接比較の出力を変更しない。

この文書の phase はレビュー可能な実装順序を示す。途中 phase の機能を Issue #239 完了とは扱わず、
全受け入れ基準を満たした時点で完成とする。

## 2. 設計上の決定

### 2.1 用語

- `record`: 一つの表示対象 `SeqRecord`。crop / reverse complement 後の表示 record を含む。
- `row`: 同じ Y 軸上に配置する一つ以上の record の集合。
- `column`: 同じ row 内の左から右への順序。
- `placement`: record の row、column、X 座標、Y 座標、描画幅、上下 extent を解決した値。
- `comparison`: query record、subject record、match table を持つ明示的な比較 edge。

コード、CLI、Web、文書で `lane` と `row` を混在させず、公開用語は既存
`--records_table` に合わせて `row` とする。内部 helper 名では `lane_layout` を使用してもよい。

### 2.2 共通スケールだけを実装する

各 record を独立して同じ幅へ引き伸ばす `per_record` mode は追加しない。

MCScanX の公式 `dual_synteny_plotter` は左右の chromosome 総長から一つの `unit` を計算し、
`bar_plotter` も全 chromosome の最大長から一つの `unit` を計算している。

- [MCScanX `dual_synteny_plotter.java`](https://github.com/wyp1125/MCScanX/blob/master/downstream_analyses/dual_synteny_plotter.java)
- [MCScanX `bar_plotter.java`](https://github.com/wyp1125/MCScanX/blob/master/downstream_analyses/bar_plotter.java)

gbdraw でも全 row に共通する `px_per_bp` を一度だけ計算する。

```text
available_width = canvas alignment width
row_bp_total    = sum(displayed record lengths in the row)
row_gap_total   = record_gap_px * (record_count_in_row - 1)
row_fixed_width = row_gap_total + resolved decorative overhangs
row_scale       = (available_width - row_fixed_width) / row_bp_total
px_per_bp       = minimum row_scale across all rows
record_width    = displayed_record_length * px_per_bp
```

装飾テキストの左右 overhang は record の bp 幅に加えず、placement の外側 inset として扱う。
これによりラベル幅が配列長の視覚表現を変えない。

既存 `normalize_length=True` は record ごとに幅を正規化するため、同じ row に複数 record がある
layout とは意味が衝突する。初回実装では、いずれかの row に複数 record があり、かつ
`normalize_length=True` の場合は `ValidationError` とする。値を黙って無視または再解釈しない。

### 2.3 record-local ruler

新しい ruler mode は作らず、既存 `ruler_on_axis` と scale style を拡張する。

- ruler は各 record の表示開始を 0 とする。
- tick interval と bp/px は全 record で共通にする。
- 各 record の末端を越える tick は描かない。
- auto interval は最長の表示 record と共通スケールから一度だけ決め、短い record は収まる tick だけを描く。
- unit label は各 tick に繰り返さず、既存の表記規則に従う。
- ruler は任意表示とし、multi-record layout を有効にしただけでは強制表示しない。

単一の共通 length bar は multi-record row の個々の座標開始点を説明できないため、
multi-record layout で ruler を表示する場合は record-local ruler を authoritative scale とする。

### 2.4 初回に対応する comparison topology

初回実装は次だけを受理する。

- query と subject は異なる record。
- query row と subject row は隣接する。
- 一つの隣接 row pair 内で任意の N 対 M record pair を選べる。
- 同じ record pair に複数 input table がある場合、共通 filter 適用後に一つの comparison へ結合する。
- query が上、subject が下であることは要求しない。q/s の意味を保ったまま geometry の方向を決める。

次は初回スコープ外とする。

- 同一 row 内を弧で結ぶ comparison。
- row を飛び越える comparison。
- crossing を減らす record 自動並べ替え。
- edge bundling、density map、ribbon aggregation。
- comparison ごとの独立した filter、legend、color theme。

同一 row や非隣接 row の edge を黙って捨てず、record selector と row を含む validation error にする。

### 2.5 visual policy

- record axis、feature、label、ruler は comparison より前面に描く。
- comparison は既存 `ribbon` / `curve` style を再利用する。
- multi-record 専用の第三の path style は追加しない。
- wide ribbon を先に、小さい ribbon を後に描く既存 draw-order policy を維持する。
- identity、orientation、collinearity color の既存設定を維持する。
- Web の hover / click 強調は既存 interactive match metadata を使い、新しい selection engine を作らない。

### 2.6 非目標

- MCScanX の解析 algorithm や file format を再実装しない。
- 一般-purpose graph layout framework を導入しない。
- Linear record 一機能のために全 renderer を plugin architecture へ置き換えない。
- record ごとの独立スケール、broken axis、log scale を追加しない。
- 長い label を自動要約、翻訳、または任意位置へ route しない。
- N 対 M の全組合せを常に自動実行しない。利用者が選択または batch action で明示する。

## 3. 現状と変更境界

### 3.1 現状

- Linear assembly は record index ごとに一つの `record_offsets` Y 座標を作る。
- `record_offsets_x` は各 record を同じ alignment area 内で左寄せまたは中央寄せする値であり、
  同じ row 内の複数 record を表現できない。
- `PairWiseMatchGroup` は `comparison_count - 1` と `comparison_count` から q/s record を推定する。
- comparison file の並び順が record N と N+1 の対応を暗黙に表す。
- `--records_table` は `row` / `column` を読み込めるが、Linear mode は値を無視する。
- `RecordPresentation` と canonical request は `gridRow` / `gridColumn` をすでに保持できるが、
  `LinearDiagramRequest` は placement を拒否する。
- Web UI は sequence card 内の `Compare to next` に BLAST file を保持する。

### 3.2 変更する範囲

- Linear record placement と canvas dimension calculation。
- record ごとの coordinate mapping と ruler placement。
- explicit comparison endpoints と comparison path geometry。
- `--records_table`、新しい comparison manifest、CLI validation。
- typed Python request、canonical request schema、session migration。
- Web の record row assignment と comparison pair editor。
- interactive SVG の comparison endpoint metadata。
- unit、integration、reference output、Web unit、browser test。

### 3.3 変更しない範囲

- BLAST outfmt 6 / 7 の列意味。
- feature、label、track slot の公開 style contract。
- crop と reverse complement の適用順。
- pairwise identity / collinearity color calculation。
- output format と export path。
- Circular multi-record layout の見た目と option semantics。

## 4. Target architecture

```text
records + RecordPresentation(row/column)
                    |
                    v
       resolve record row/order selectors
                    |
                    v
       calculate LinearRecordPlacement[]
       - one shared px_per_bp
       - row widths and offsets
       - record-local coordinate mapper
       - top/bottom comparison anchors
                    |
          +---------+---------+
          |                   |
          v                   v
  record/track renderers   explicit comparisons
          |                   |
          |            resolve q/s placements
          |                   |
          |                   v
          |          PairWiseMatchGroup
          |          receives endpoints
          +---------+---------+
                    |
                    v
          static / interactive SVG
```

assembly は orchestration、layout module は pure geometry、render group は SVG construction を所有する。
file parsing、selector resolution、geometry、drawing を一つの class にまとめない。

## 5. Domain contracts

### 5.1 layout option

Circular option を拡張して mode 判定を増やすのではなく、Linear 専用 option を追加する。

```python
@dataclass(frozen=True)
class LinearMultiRecordOptions:
    record_gap_px: float = 24.0
    multi_record_positions: Sequence[str] | None = None
```

- `record_gap_px` は有限の 0 以上。
- `multi_record_positions` は既存 Circular と同じ `<selector>@<row>` 文法を使う。
- 同じ row 内の順序は position token の順序で決める。
- `--records_table` の `column` は token の順序へ正規化し、renderer へ column 値を直接渡さない。
- record を一度だけ含める、全 loaded record を含める、selector を一意に解決する、という既存検証を共有する。

`LinearDiagramRequest` に optional `layout: LinearMultiRecordOptions | None` を追加する。
`RecordPresentation.grid_row` / `grid_column` は Circular と Linear の双方で使える mode-neutral contract とし、
mode 固有 validation は request validator に置く。

直接 API の `assemble_linear_diagram_from_records()` には `layout` bundle を追加し、row list、column list、
gap list のような parallel arguments は追加しない。

### 5.2 calculated placement

内部 geometry は一つの immutable value object にまとめる。

```python
@dataclass(frozen=True)
class LinearRecordPlacement:
    record_index: int
    row: int
    column: int
    x: float
    axis_y: float
    sequence_width: float
    left_inset: float
    right_inset: float
    top_extent: float
    bottom_extent: float
    comparison_top_y: float
    comparison_bottom_y: float
    px_per_bp: float

    def x_for_position(self, position: float) -> float: ...
```

`x_for_position()` は crop / reverse complement 後の local displayed coordinates だけを受け取る。
公開 1-based 座標や source coordinates の変換を renderer に持ち込まない。

`left_inset` / `right_inset` は record definition、subtitle、端の ruler label の overhang を含む。
sequence width と分離し、文字列の長さで bp scale が変化しないようにする。

### 5.3 explicit comparison

order-dependent `DataFrame` list を drawing owner へ渡さず、endpoint を持つ型へ正規化する。

```python
@dataclass(frozen=True)
class LinearComparison:
    query_record_index: int
    subject_record_index: int
    matches: DataFrame
```

- index は input record sequence に対する 0-based index。
- row / column による表示順変更後も index の意味を変えない。
- `matches` は既存 comparison columns と metadata columns を保持する。
- filter は既存 `filter_comparison_dataframe()` を一度だけ適用する。
- renderer は DataFrame 内の qseqid / sseqid から placement を再推定しない。
- duplicate record ID があっても endpoint index で曖昧にならない。

`DiagramOptions` には `linear_comparisons: Sequence[LinearComparison] | None` を追加する。
既存 `blast_files` / `protein_comparisons` は後方互換 input として残し、boundary adapter で
`LinearComparison` へ変換する。

### 5.4 comparison table

CLI に `--comparisons_table TSV` を追加する。初回 schema は必要最小限にする。

```tsv
blast	query	subject
a_vs_c.tsv	#1	#3
b_vs_c.tsv	#2	#3
b_vs_d.tsv	#2	#4
```

必須 column:

```text
blast  query  subject
```

- path は table file からの相対 path を許可し、既存 table dependency machinery を使う。
- query / subject は unique record ID または `#index`。
- unknown column、missing file、unresolved selector、same-record、same-row、non-adjacent-row を
  table path、row number、column name 付き `ValidationError` にする。
- `--comparisons_table` と legacy `-b/--blast` は同時指定不可。
- label、color、filter override は初回 schema に追加しない。

legacy `-b` は次の条件で維持する。

- 各 row が一つの record だけなら、従来どおり file N を隣接 row N と N+1 へ割り当てる。
- いずれかの row に複数 record がある場合、endpoint が曖昧な `-b` を拒否し、
  `--comparisons_table` を案内する。

### 5.5 canonical request schema

canonical request は strict schema のため、comparison entry に endpoint field を追加する際に
schema 2 へ更新する。

```json
{
  "kind": "nucleotideBlast",
  "resourceId": "comparison-nucleotide-1",
  "queryRecordIndex": 0,
  "subjectRecordIndex": 2
}
```

precomputed protein comparison にも同じ endpoint fields を持たせる。generated comparison は
selected pair list を持てる形にする。

```json
{
  "kind": "generatedProteinComparison",
  "mode": "pairwise",
  "pairs": [
    {"queryRecordIndex": 0, "subjectRecordIndex": 2},
    {"queryRecordIndex": 1, "subjectRecordIndex": 2}
  ],
  "settings": {}
}
```

- schema 2 encoder は explicit endpoints を必ず書く。
- schema 1 decoder は従来の array order から adjacent endpoints を補う。
- schema 1 を復元して再保存した場合は schema 2 を出力する。
- unknown endpoint field を黙って無視しない strict policy を維持する。
- record reorder 後も endpoint が同じ logical record を参照することを round-trip test で固定する。

## 6. Geometry and rendering

### 6.1 row construction

1. record selector を input index へ解決する。
2. position token 順に `(row, column, record_index)` を作る。
3. row ごとに record を group 化する。
4. displayed sequence length、definition bbox、track extent、label extent を事前計算する。
5. 全 row から一つの `px_per_bp` を計算する。
6. row content width を計算する。
7. `align_center=True` は row 全体を alignment area の中央へ置く。
8. `align_center=False` は row 全体を alignment area の左へ置く。

既存 `align_center` を record 単位と row 単位の二重処理にしない。multi-record layout では
row group の alignment だけを意味する。

### 6.2 vertical extents

row の上下 extent は、その row に属する record の最大 extent とする。

```text
row_top_extent    = max(record top extents)
row_bottom_extent = max(record bottom extents)
next_axis_y       = current_axis_y
                    + row_bottom_extent
                    + comparison_height
                    + next_row_top_extent
```

definition、feature labels、custom slots、GC、skew、depth、annotation を extent に含める。
比較開始点と終了点は record ごとの `comparison_bottom_y` / `comparison_top_y` を使い、
row の最大 extent を全 endpoint に強制しない。

### 6.3 record renderers

`LinearCanvasConfigurator` は canvas 全体の style と available dimensions を所有し、
record ごとの幅や位置を mutable field に切り替えない。

各 record / quantitative track builder は必要な `LinearRecordPlacement` または最小の
`sequence_width` / coordinate mapper を受け取る。record ごとに canvas configurator を clone して
一時的に `alignment_width` を書き換える設計は採用しない。

次の renderer を同じ coordinate mapping path へ移す。

- sequence axis と feature
- feature label と leader line
- GC content / GC skew
- depth track
- annotation track
- axis ruler
- pairwise match endpoints

### 6.4 definition placement

multi-record row では record definition を左端の共通 column に置かない。各 record の上側に
record-local header として置く。

- label と subtitle の文字列生成、font、style は既存 `DefinitionGroup` を再利用する。
- position は record sequence span の中央を基準とする。
- header bbox の overhang を placement inset に含める。
- header 高さを row top extent に含める。
- one-record-per-row の definition placement は変更しない。

新しい definition text renderer を複製せず、既存 group に placement mode を追加する。

### 6.5 comparison rendering

`PairWiseMatchGroup` から `comparison_count` による record 推定を除く。constructor は
`LinearComparison` と query / subject placement を明示的に受け取る。

```text
query x   = query_placement.x_for_position(q coordinate)
subject x = subject_placement.x_for_position(s coordinate)
query y   = query placement comparison anchor
subject y = subject placement comparison anchor
```

SVG metadata の `data-query-record-index` / `data-subject-record-index` も同じ explicit endpoints から
設定する。metadata と geometry で別々の index 推定をしない。

## 7. SOLID、KISS、DRY の適用

### 7.1 SOLID

- SRP: table reader、selector resolver、lane geometry、comparison normalization、SVG drawing、Web state を
  分離する。`assemble.py` は処理順を組み立てるだけにする。
- OCP: renderer は `LinearRecordPlacement` を受け取ることで、row 数や同一 row の record 数に依存しない。
  新しい layout のたびに feature / depth / ruler の式を分岐させない。
- LSP: one-record-per-row の placement は既存位置、幅、DOM metadata と置換可能でなければならない。
  reference output の意図しない差分を許容しない。
- ISP: drawer に records table、request、app state 全体を渡さない。描画に必要な placement、record、style、
  match rows だけを渡す。
- DIP: assembly は pure layout function の出力に依存し、Web state、CLI namespace、file path を直接参照しない。
  このためだけの abstract base class や dependency container は導入しない。

### 7.2 KISS

- 共通スケール一つだけを実装する。
- 初回 comparison は隣接 row 間だけにする。
- row 内 order は利用者が決め、自動最適化しない。
- Web は pair matrix ではなく単純な comparison list を使う。
- comparison table の初回 column を三つに限定する。
- 既存 `curve` / `ribbon`、ruler、selector、records table、filter を再利用する。
- multi-record layout が無効なら新しい code path を最小限にし、既存 output を維持する。

### 7.3 DRY

- Circular と Linear で `<selector>@<row>` parser と selector resolution を共有する。
- CLI、API、session、Web は最終的に同じ `LinearComparison` endpoint contract へ正規化する。
- coordinate-to-X conversion は `LinearRecordPlacement.x_for_position()` の一箇所に置く。
- comparison endpoint index は geometry、metadata、popup で同じ値を使う。
- comparison filtering と concatenation は input owner で一度だけ行う。
- record header の text measurement と style は既存 `DefinitionGroup` を再利用する。
- legacy adjacency は boundary adapter 一箇所だけに残し、renderer に互換分岐を置かない。

## 8. CLI and Python API

### 8.1 CLI

Linear command に次を追加する。

```text
--multi_record_position SELECTOR@ROW   repeatable
--linear_record_gap PX
--comparisons_table TSV
```

`--records_table` が `row` / `column` を持つ場合は、同じ position contract へ変換する。

- `--records_table` と `--multi_record_position` は同時指定不可。
- position を一部 record だけに指定しない。
- `column` は同じ row 内で一意。
- row number は連続でなくても入力時に受理し、内部で 0-based contiguous row へ正規化する。
- comparison を持たない multi-record row も許可する。

### 8.2 Python API

- `LinearMultiRecordOptions` を `gbdraw.api` から export する。
- `LinearComparison` と comparison table reader を `gbdraw.api` から export する。
- typed `LinearDiagramRequest.layout` を authoritative request contract とする。
- long-form `assemble_linear_diagram_from_records()` は同じ option object を受ける。
- raw `blast_files` / `protein_comparisons` は deprecated にせず、legacy adjacency adapter を通す。

parallel list や dict の shape を複数追加しない。

## 9. Web UI

### 9.1 record layout editor

Linear input cardとは別に小さな `Record Layout` section を置く。

- 各 record row に表示名、row selector、同じ row 内の Up / Down を表示する。
- column number を直接編集させず、list order から導出する。
- record identity は array index ではなく既存 sequence `uid` を使う。
- `Auto` は従来どおり一 record 一 row。
- `Arrange in rows` を有効にしたときだけ row editor を表示する。
- Circular record order helper の selector表示と move logic を再利用できる部分は共通化する。

### 9.2 comparison editor

sequence card 内の `Compare to next` upload を、独立した comparison list へ移行する。

```text
From [record A]  To [record C]  Source [Upload / LOSAT]  [Remove]
```

初回 batch action:

- `Add adjacent pairs`: 各隣接 row の同じ順序の record pair を可能な範囲で追加する。
- `Add all pairs between adjacent rows`: 各隣接 row pair の N x M を追加する。
- `Clear comparisons`。

batch action 実行前に追加 pair 数を表示する。LOSAT の N x M job を暗黙に開始しない。

Web state は focused module に分離する。

```text
linear-record-layout.js     row assignment and order
linear-comparisons.js       endpoint pairs and source state
```

`run-analysis.js` は正規化済み layout / comparisons を args または canonical request へ渡すだけにする。
新しい CDN dependency や build step は追加しない。

### 9.3 legacy session migration

旧 `linearSeqs[index].blast` は record index `index` と `index + 1` の comparison へ変換する。
旧 session の record layout は一 record 一 row とする。migration 後は新しい comparison list と
canonical schema 2 で保存する。

### 9.4 interactive SVG

既存 popup は explicit endpoint metadata を読む。追加する SVG metadata は必要最小限にする。

- record group: row、column、record index。
- comparison group/path: query index、subject index、query row、subject row。

座標、identity、protein / orthogroup metadata は既存属性を維持する。

## 10. Implementation phases

### Phase 0: characterization tests

実装前に次を固定する。

1. one-record-per-row の record Y、X、width、definition、ruler geometry。
2. BLAST file N が record N / N+1 に接続される legacy metadata。
3. `ribbon` と `curve` の path geometry。
4. crop、reverse complement、align center、normalize length の既存出力。
5. Linear `RecordPresentation.grid_row` が現在は validation error になること。
6. schema 1 comparison round-trip と旧 Web session の `seq.blast` 復元。

完了条件:

- reference output を更新せず test が通る。
- 新しい layout contract 導入後に守るべき既存値が assertion になっている。

### Phase 1: shared placement contracts

1. selector / row resolver を Circular 固有 helper から mode-neutral helper へ抽出する。
2. `LinearMultiRecordOptions` と validation を追加する。
3. `LinearRecordPlacement` と pure row grouping / common scale calculation を追加する。
4. definition bbox と track extents を placement input に統合する。
5. one-record-per-row の既存 geometry を維持する compatibility test を通す。

完了条件:

- 2 + 3 record の二つの row について決定的な placement が得られる。
- record 順、row 順、gap、center / left alignment が unit test で説明できる。
- SVG rendering をまだ変更せず geometry module を単体検証できる。

### Phase 2: record and ruler rendering

1. record / track builder に placement または coordinate mapper を渡す。
2. axis、feature、label、quantitative track、annotation を record-local width で描く。
3. multi-record definition を record-local header として配置する。
4. `ruler_on_axis` を record-local ruler として描く。
5. canvas width、height、legend offset、plot title、bottom content の計算を更新する。

完了条件:

- comparison なしの 1xN、2xN、3-row layout が clip されない。
- shared ruler interval と proportional record width を SVG attribute / coordinates で検証できる。
- feature、GC、depth、annotation の X coordinate が同じ mapper を使う。
- one-record-per-row reference output に意図しない差分がない。

### Phase 3: explicit comparison endpoints

1. `LinearComparison` と normalization adapter を追加する。
2. `PairWiseMatchGroup` から `comparison_count` endpoint inference を除く。
3. q/s placement から path X/Y と metadata を作る。
4. same-row / non-adjacent-row validation を追加する。
5. `--comparisons_table` reader と relative dependency handling を追加する。
6. legacy `-b` adjacency adapter を追加する。

完了条件:

- 2x2、2x3、3-row の選択 pairだけが正しい recordへ接続される。
- duplicate record ID でも index metadata が正しい。
- crop / reverse complement 後の endpoint が表示 record の座標へ一致する。
- interactive match popup が正しい二つの record を表示する。

### Phase 4: generated protein comparisons

1. pairwise protein workflow が explicit pair list を受け取れるようにする。
2. multi-record row の default adjacent scope を、隣接 row 間の選択 pairとして定義する。
3. orthogroup / collinear result を explicit record index付き edgeへ変換する。
4. multi-record layout では `collinear_search_scope=all` と non-adjacent rendered edge の関係を
   validation し、初回は隣接 rowへ限定する。
5. LOSAT cache key に q/s record identity を含める。

完了条件:

- uploaded BLAST、LOSAT pairwise、orthogroup、collinear の全経路が同じ renderer input を生成する。
- N x M を選択しなければ不要な LOSAT jobを実行しない。
- cache が record reorder 後に別 pair の結果として再利用されない。

### Phase 5: CLI, typed request, and session schema

1. Linear CLI options と records/comparisons table validation を追加する。
2. `LinearDiagramRequest.layout` と request rendering を追加する。
3. canonical request schema 2 encoder / decoder と schema 1 migration を追加する。
4. session dependency collection に comparisons table と BLAST resources を含める。
5. public API exports と docstrings を追加する。

完了条件:

- CLI、typed request、canonical request が同じ SVG geometry を生成する。
- schema 1 session が復元でき、schema 2 として再保存できる。
- table path dependencies が session bundleへ埋め込まれ、復元時に書き換えられる。

### Phase 6: Web UI

1. record layout state / editor を追加する。
2. comparison list、upload、LOSAT pair、batch actions を追加する。
3. canonical request builder と legacy snapshot migration を更新する。
4. history / undo / redo で row move と comparison add/remove を一操作として扱う。
5. pair count、invalid topology、missing source を Generate 前に表示する。
6. hover / click popup と feature focus を explicit endpoints で検証する。

完了条件:

- 2x2 layout を作成し、三つの選択 comparisonを描ける。
- save / load、undo / redo、file reorder 後も endpoint が維持される。
- offline / local GUI で外部 dependencyなしに動作する。
- narrow viewport でも editor の selector と buttons が操作できる。

### Phase 7: verification and documentation

1. focused unit、integration、reference、Web、browser test を実行する。
2. intentional SVG geometry 用の reference output を一つ追加し、既存 reference は必要な場合だけ更新する。
3. Tutorials 2、4、5、7、CLI Reference、Recipes、Gallery を更新する。
4. records table と comparisons table の copy-pasteable example を追加する。
5. 実際の生成 SVG を使った multi-record example を追加する。

## 11. Test plan

### 11.1 pure layout tests

新規 `tests/test_linear_multi_record_layout.py` で検証する。

- one-record rows の legacy offsets。
- 2 records / 3 records の row grouping。
- common `px_per_bp` と proportional width。
- fixed pixel gap。
- center / left row alignment。
- definition / ruler overhang inset。
- row top / bottom extent。
- non-contiguous input row number normalization。
- duplicate / missing placement validation。
- multi-record row と `normalize_length=True` の rejection。

### 11.2 comparison tests

新規 `tests/test_linear_multi_record_comparisons.py` で検証する。

- explicit q/s index normalization。
- same row、same record、non-adjacent row rejection。
- multiple files for one pair の concatenation。
- identity / bitscore / evalue / length filter の一回適用。
- curve / ribbon endpoint geometry。
- reverse orientation。
- duplicate record ID。
- crop / reverse complement。
- interactive SVG endpoint metadata。

### 11.3 CLI and table tests

- records table row / column が Linear placementになる。
- comparisons table relative path resolution。
- query / subject ID と `#index`。
- table row / column付き error。
- `--comparisons_table` と `--blast` の排他。
- multi-record row に ambiguous legacy `--blast` を渡した場合の error。
- one-record rows の legacy `--blast` success。

### 11.4 API and session tests

- `LinearDiagramRequest.layout` validation。
- RecordPresentation row / column round-trip。
- `LinearComparison` typed request rendering。
- schema 2 encode / decode。
- schema 1 comparison migration。
- record reorder と endpoint identity。
- session embedded resource restoration。

### 11.5 Web unit tests

- uid-based row assignment。
- move within row / move across row。
- comparison pair add / remove / deduplicate。
- adjacent pair batch action。
- all-pairs count。
- invalid topology message。
- legacy `seq.blast` migration。
- history snapshot / canonical request round-trip。
- removed record に接続していた comparison cleanup。

### 11.6 browser tests

最低限、次の flow を実ブラウザで確認する。

1. 4 records を upload する。
2. #1 / #2 を row 1、#3 / #4 を row 2 に置く。
3. record-local ruler を有効にする。
4. #1 -> #3、#2 -> #3、#2 -> #4 を追加する。
5. BLAST upload または LOSAT を実行する。
6. curve ribbon が正しい record 区間へ接続される。
7. match をクリックし、q/s record と座標を確認する。
8. session を保存、復元し、layout と comparison pair が一致する。
9. SVG、interactive SVG、PNG、PDF を export する。

Node `@playwright/test` がなければ Python Playwright で同じ flow を検証する。

### 11.7 reference output

少なくとも一つの小さい deterministic fixture を追加する。

```text
Row 1: record A (10 kb), record B (5 kb)
Row 2: record C (8 kb), record D (7 kb)
Comparisons: A-C, B-C, B-D
```

fixture は shared scale、record gap、local ruler、curve path、metadata を一つの SVG で確認できるようにする。
既存 `tests/reference_outputs/` は read-onlyとして比較し、意図した geometry 変更だけを
`--update-reference-outputs` で更新して差分をレビューする。

## 12. Compatibility and migration

| Existing workflow | New behavior |
|---|---|
| One record per row, no comparison | 変更なし |
| One record per row with ordered `-b` files | 従来どおり adjacent comparison |
| `--records_table` without row / column | 変更なし |
| `--records_table` with row / column in Linear | placementとして有効化 |
| Multi-record row with raw `-b` files | ambiguity error、comparisons tableを案内 |
| Schema 1 canonical request | adjacent endpointsを補って読込 |
| Legacy Web `seq.blast` session | adjacent comparison listへmigration |
| `normalize_length=True`, one record per row | 変更なし |
| `normalize_length=True`, multi-record row | validation error |

Circular の `multi_record_size_mode`、gap ratio、multi-record canvas は変更しない。

## 13. Risks and mitigations

### 13.1 coordinate drift

Risk: feature、ruler、comparison が別々の width formula を使い、endpoint がずれる。

Mitigation: X conversion を `LinearRecordPlacement.x_for_position()` に集約し、全 renderer の同じ既知座標を
testする。

### 13.2 crop and reverse complement

Risk: BLAST coordinates と表示 record coordinates の空間が一致しない。

Mitigation: 既存と同じく比較 input は表示座標系に一致することを要求し、source-to-local transformationを
暗黙に推測しない。validation / warning 文言を docs と Web に表示する。

### 13.3 label overlap

Risk: 短い record の長い header が隣の record と重なる。

Mitigation: measured bbox overhang を left/right inset として row width に含める。font縮小や文字列切り捨てを
暗黙に行わない。

### 13.4 dense comparisons

Risk: N x M の ribbon 数と LOSAT job 数が急増する。

Mitigation: pairを明示し、batch action 前に件数を表示し、既存 filterを適用し、不要なpairを生成しない。
性能問題が計測された後にaggregationを別Issueとして検討する。

### 13.5 endpoint identity after reorder

Risk: array indexだけでWeb stateを保持すると、record移動後に別recordへ接続される。

Mitigation: Web編集中はuid、canonical/request/render境界ではstable input indexを使い、変換を一箇所に置く。

### 13.6 broad renderer changes

Risk: record-local width対応がfeature、GC、depth、annotationの既存出力へ波及する。

Mitigation: Phase 0 characterization、placementのoptional path、one-record legacy assertions、focused reference
diffを先に用意する。一度に全rendererを書き換えず、同じmapper contractへ順に移す。

## 14. Completion criteria

Issue #239 は次をすべて満たしたとき完了とする。

- [ ] Linear mode で同じ row に複数 record を配置できる。
- [ ] 全 record が一つの共通 bp/px スケールを使う。
- [ ] record-local ruler を既存 option で表示できる。
- [ ] selected N-to-M comparisons が隣接 row 間の正しいrecordへ接続される。
- [ ] feature、label、GC、skew、depth、annotation、custom track slot が配置に追従する。
- [ ] CLI records/comparisons table、Python API、Web が同じ結果を生成する。
- [ ] crop、reverse complement、duplicate ID が明示的なcontractとtestで扱われる。
- [ ] schema 1 / legacy Web sessionを復元できる。
- [ ] one-record-per-row の既存 SVG に意図しない差分がない。
- [ ] SVG、interactive SVG、PNG、PDF exportを確認している。
- [ ] focused tests、reference comparison、browser flowが通る。
- [ ] Tutorials、CLI Reference、Recipes、Galleryが実装と一致する。

## 15. 推奨する最初のPR境界

最初のPRは Phase 0 と Phase 1 に限定する。

```text
PR 1: characterization + placement contracts
PR 2: record/ruler rendering
PR 3: explicit comparisons + tables
PR 4: generated protein comparisons + request/session schema
PR 5: Web UI + browser verification + docs
```

各PRで後続phaseの仮UIや一時的な公開optionを追加しない。内部contractが安定する前に
Web stateとsession schemaを先行させず、rendererがorder-dependent DataFrame listを受ける期間を
長引かせない。
