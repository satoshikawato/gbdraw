## gbdraw リファクタリング提案（ライブラリ化・マルチトラック・カスタムトラック・インタラクティブ）

このドキュメントは、`gbdraw` を今後「パイプラインに組み込めるライブラリ」として拡張しつつ、Circosのような**任意トラックを任意位置へ配置できるマルチトラック（環状ゲノム含む）**と、**多様なカスタムデータトラック（SAM/BED/Wig/VCFなど）**、さらに**画面上で直接いじれるインタラクティブ化**へ発展させるための設計提案です。

（実装の“全部置き換え”ではなく、現行の `linear/circular` 描画ロジックを活かしながら段階的に到達する案を中心に書きます）

---

### 目標（ユーザー要望の再整理）

- **コード構成の合理化**
  - 凝集度を上げ、責務分離（IO / 設定 / レイアウト / 描画 / UI）を明確化
  - “拡張しやすい境界（API）”を作る
- **ライブラリとして使えるように**
  - CLIだけでなく、Python API（関数/クラス）で呼べる
  - 例：Snakemake/Nextflow、社内パイプライン、Webサービスに組み込み可能
- **環状ゲノムも含めたマルチトラック化**
  - Circos的な「任意のトラックを任意位置（半径/角度/順序）に配置」へ
  - ただし Circos のルール表記は複雑なので、**使いやすい設定DSL**へ
- **カスタムデータトラック**
  - BED/GFF/WIG/VCF/SAM(BAM)等の外部データを載せる
  - 大きいファイル（BAM/BigWig/VCF.gz）も扱える道筋（インデックス・ダウンサンプリング）
- **インタラクティブ化**
  - Web/GUI上で、トラック順序や位置、ラベル位置、色、フィルタなどを直接編集
  - 編集結果を設定として保存し、再現可能にする

---

### 現状（ざっくり俯瞰：どこがボトルネックか）

現行 `gbdraw` はすでに「CLI + Web UI（Pyodide）」を持っており、描画も十分リッチです。一方で、拡張（特にマルチトラックと外部データ）をするには“構造”が足りません。

- **エントリポイント**
  - CLI: `gbdraw/cli.py` → `gbdraw/circular.py` と `gbdraw/linear.py`
  - Web: `gbdraw/web/index.html`（Vue + Pyodide + micropip で gbdraw wheel を導入）
- **描画パイプライン（現状）**
  - `circular.py` / `linear.py` が「引数→設定ロード→入力ロード→Configurator生成→plot_*」まで担当
  - `circular_diagram_components.py` / `linear_diagram_components.py` が “決まった構成のトラック（features + GC + skew + …）” を順番にcanvasへ積む
  - 出力は主に SVG（変換はCLIで CairoSVG、Webはブラウザ側変換）
- **マルチトラック（環状）の前提が固定**
  - `config.toml` に「最大3トラック前提の track_ratio / track_dict」的な構造があり、`track_type=tuckin/middle/spreadout`で位置が決まる
  - つまり「任意数トラック」「任意半径帯」「任意順序」「任意角度範囲（部分円）」を表現できない
- **IO/ユーティリティが巨大で、責務が混ざっている**
  - `file_processing.py` に GenBank/GFF、色、保存、pyodide分岐などが同居
  - `utility_functions.py` が肥大化し、ラベル配置などが別ファイルと重複している箇所もある
- **ライブラリ利用がしづらい**
  - `sys.exit()` 前提のエラーハンドリングが多く、外部から呼び出すと例外で拾いづらい
  - “関数APIが安定していない”ので、他人がパイプラインに組み込む前提が弱い

---

### 提案：ターゲットアーキテクチャ（分離する境界）

今後の拡張要件を満たすために、`gbdraw` を以下のレイヤに分離します。

- **core（ドメイン）**: Diagram/Genome/Track/Style/Config など “概念と整合性”
- **io（入出力）**: GenBank/GFF/FASTA + カスタムトラックデータ読み込み、インデックス、バリデーション
- **layout（レイアウト）**: linear/circular の座標変換・トラック配置・ラベル配置・衝突回避
- **render（描画）**: SVG/（将来）Canvas/JSON scene graph などのバックエンド
- **ui（CLI/Web）**: コマンドライン引数やWebフォームは core を呼ぶだけにする
- **plugins（拡張）**: トラック追加・データソース追加を外部パッケージから可能にする

#### ディレクトリ案（例）

```text
gbdraw/
  __init__.py
  api.py                 # 安定API（Diagram生成/レンダリングの入口）
  exceptions.py          # Library向け例外（CLIはここを捕まえてsys.exit）

  core/
    model.py             # Genome/Region/Interval/Variant/Series/Link など中間表現
    diagram.py           # Diagram + TrackRegistry + Render pipeline
    config.py            # Typed config（pydantic/dataclasses） + schema versioning
    theme.py             # パレット、スタイル、既定値

  io/
    genome.py            # GenBank/GFF+FASTA 読み込み（現 load_gbks/load_gff_fasta の整理）
    tracks/
      bed.py
      wig.py
      vcf.py
      sam.py              # (server mode) pysam前提
    registry.py          # format推定、DataSource登録

  layout/
    circular.py          # 角度変換・半径帯・トラック配置
    linear.py            # x/y配置・行/レーン・高さ
    labels.py            # label placement / collision
    scales.py            # tick/scale bar/ruler

  render/
    svg.py               # svgwrite backend（現状資産を集約）
    scene.py             # (将来) interactive用のscene graph生成

  cli/
    main.py              # typer/click で引数→core.api呼び出し

  web/
    index.html           # 現状継続（serverless pyodide）
    ...                  # (将来) scene graphを扱うUI
```

この形にすると、将来 “別UI（Streamlit/Jupyter/サーバ）” を増やしても **core/layout/render/io は使い回し**できます。

---

### 提案：ライブラリAPI（安定APIの形）

「外部がパイプラインに組み込む」ためには、CLIオプションの寄せ集めより **小さく安定したAPI** が必要です。

#### 最小API（まずこれを固定する）

- `gbdraw.api.render(config: dict|Config, *, outputs=[...]) -> RenderResult`
- `gbdraw.api.load_config(path) -> Config`
- `gbdraw.api.save_config(config, path)`
- `Diagram` クラス（高度利用）

#### 例：Pythonからの利用イメージ（案）

```python
from gbdraw.api import Diagram
from gbdraw.io.genome import read_genome
from gbdraw.io.tracks.bed import read_bed
from gbdraw.tracks import AnnotationTrack, IntervalTrack
from gbdraw.layout.circular import CircularLayout

genome = read_genome("NC_000913.gb")
bed = read_bed("peaks.bed")

diagram = Diagram(
  genome=genome,
  layout=CircularLayout(start_angle=-90, end_angle=270),
)

diagram.add_track(AnnotationTrack(), radius=(0.85, 1.00))
diagram.add_track(IntervalTrack(bed), radius=(0.70, 0.82), style={"color": "#ef4444"})

diagram.save("out.svg")
```

ポイント：
- 位置指定（`radius=(r0,r1)` など）が **第一級** になる
- “トラックを足す”がAPIの中心（Circos思想）
- CLI/Webはこの API を呼ぶだけにする

---

### 提案：設定DSL（Circos風だがシンプル）

Circosのような「ルール駆動」は強いですが、学習コストが高いです。`gbdraw` では以下を狙います。

- **YAML/JSONを第一候補**（Webでも扱いやすい）
- 「トラックを積む」概念をそのまま書ける
- 位置指定は“絶対値”より“相対（0..1）”を基本にし、必要ならpx/ptも許容
- 将来の互換性のため **schema_version** を入れる

#### 例：circular + 任意トラック

```yaml
schema_version: 1
diagram:
  layout: circular
  genome:
    source: genome.gb
    topology: auto  # auto|circular|linear
  canvas:
    width: 1200
    height: 1000
    start_angle: -90
    end_angle: 270

  tracks:
    - id: annotation
      type: annotation
      radius: [0.86, 1.00]
      strand: split            # split|merged
      labels:
        mode: none             # none|auto|all

    - id: gc
      type: gc_content
      radius: [0.74, 0.84]
      window: auto
      step: auto

    - id: variants
      type: vcf
      file: sample.vcf.gz
      radius: [0.66, 0.72]
      glyph: lollipop          # lollipop|ticks|heatmap
      filters:
        min_qual: 30

    - id: depth
      type: wig
      file: depth.wig
      radius: [0.52, 0.64]
      scale:
        min: 0
        max: 200
```

#### 例：linear + マルチトラック（行の積み重ね）

```yaml
schema_version: 1
diagram:
  layout: linear
  genomes:
    - id: g1
      source: a.gb
    - id: g2
      source: b.gb

  tracks:
    - id: features
      type: annotation
      height: 80
      labels: { mode: first } # none|first|all
    - id: gc
      type: gc_content
      height: 20
    - id: links
      type: blast_links
      file: a_b.blast.outfmt6
      height: 60
```

---

### 提案：トラックモデル（任意トラックを増やすための核）

**Track = データ + ルール + 描画** をカプセル化します。

#### 基本インターフェイス案

- **DataSource**（IO）
  - `query(region, resolution) -> TrackData`
  - `iter_all()`（小さいデータ向け）
  - （大きいデータは index + query が基本）
- **Track**（意味）
  - `prepare(context)`: 事前計算（色、スケール、集約）
  - `render(renderer, layout, context)`: 描画
- **Renderer**（出力）
  - SVG / JSON scene / Canvas（将来）

#### 内部の共通データ型（例）

- **Interval**: BED/GFF由来（start, end, strand, attributes）
- **Point**: SNP/variant など（pos, value, attributes）
- **Series**: coverage や GC%（(start,end,value) or binned array）
- **Link**: BLAST/PAF/MAF/シンテニー（regionA ↔ regionB）

この中間表現があると、同じデータを circular/linear どちらにも載せられます（座標変換は layout 側の責務）。

---

### 提案：カスタムトラック（BED/Wig/VCF/SAM…）の現実的な実装戦略

データ形式ごとに「serverless（Pyodide）」で扱えるかが変わります。

#### 1) “純Pythonでいける”枠（serverless でも対応しやすい）

- **BED**: テキスト、純PythonでOK
- **WIG**: テキスト、純PythonでOK（大きい場合はダウンサンプル必須）
- **VCF（小規模）**: 純PythonパーサでOK（ただし巨大VCFは厳しい）

#### 2) “C拡張依存が強い”枠（ローカルサーバ/通常Pythonで対応）

- **BAM/SAM/CRAM**: `pysam`（htslib）前提が現実的
- **BigWig/BigBed**: `pyBigWig`（C拡張）前提が現実的
- **VCF.gz + index**: `pysam.VariantFile` / `cyvcf2` が現実的

#### 3) だから推奨：Web UIは「serverless」と「server mode」の二刀流

- **serverless（現状のPyodide）**: 小さい入力・純Python依存で動く
- **server mode（ローカル/クラスタ上Python）**: 重い依存・巨大ファイル・インデックス問い合わせ

UIは同じでも “計算実行場所” を切り替えられる設計にすると、機能拡張が詰まりにくいです。

---

### 提案：環状マルチトラック（Circos的配置）の仕様案

やりたいことは「トラックを任意の半径帯に置く」だけではなく、以下まで含まれます。

- **任意数トラック**
- **任意の半径帯（内外）**
- **任意の描画順（z-order）**
- **任意の角度範囲（全周/部分）**
- **スケール・グリッド・ラベルの独立トラック化**
- **複数配列（複数コンティグ）を“ideogram”として並べる**（将来）

#### 最小仕様（Phase 1）

- 全周1コンティグ（既存と同等）を前提に
  - `radius: [0..1]` 指定で任意トラックを積める
  - `order`（描画順）を指定できる
  - ラベルは“別トラック”として扱える（labelの手動位置指定も可能に）

#### 拡張仕様（Phase 2+）

- 複数コンティグ/複数配列を “セグメント” として配置
  - `segments: [{id, length, gap, color, label}]`
  - `region` は `segment_id + coordinate` を持てる
- chord/link を別レイヤで描く
  - BLAST/PAF/MAFのマッチをリンクとして表現

---

### 提案：インタラクティブ化（“画面上で直接いじれる”）

現状のWeb UIは「フォーム → 再生成」のワークフローです。ここに “直接編集” を足すには、出力をただのSVGにするのではなく、**編集可能なメタ情報**を持たせるのが近道です。

#### 方針A：SVGを“編集可能”にする（最短ルート）

- SVG要素を **feature単位 / track単位で `<g>` に分ける**（既存ROADMAPにも近い）
- それぞれに `data-id`, `data-track`, `data-type` 等を付与
- UI側で
  - クリックで選択・ハイライト
  - ラベルをドラッグ（位置をoverrideとして保存）
  - trackのドラッグで順序変更（→ config更新 → 再レンダ）

メリット：既存のSVGレンダリング資産を最大限使える  
デメリット：SVG編集の自由度/パフォーマンスに限界（超巨大データは重い）

#### 方針B：Scene Graph（JSON）を出す（中期の本命）

- Python側で「図形（path/rect/text）+ 属性 + バインドデータ」を **JSON** としても出力
- Web側は Canvas/WebGL/SVG のいずれでも描画可能
- 編集結果は JSON/Config に反映して再現可能

メリット：インタラクティブに強い、将来のバックエンド切替が容易  
デメリット：初期実装コストは上がる

#### 推奨：A → B の段階的導入

まずAで「ID付与 + グルーピング + ラベル手動調整」を成立させ、その後Bへ拡張すると破綻しにくいです。

---

### 提案：移行計画（破壊的変更を管理する）

大規模リファクタは “動き続ける” ことが重要です。以下の段階で進めるのが安全です。

#### Phase 0（準備）

- `gbdraw.api` を追加し、内部は現行関数を呼ぶだけでも良いので **ライブラリ入口を固定**
- 例外型 `gbdraw.exceptions` を用意（CLIは捕まえて `sys.exit(1)`、ライブラリは例外を投げる）

#### Phase 1（責務分離）

- `file_processing.py` を `io/` と `render/` に分割
- `utility_functions.py` の巨大化を止め、重複関数を整理（labels/geometry/config へ）
- 既存 CLI の挙動は維持（内部呼び出し先を差し替える）

#### Phase 2（トラック化：最小の任意トラック）

- circular/linear ともに `Track` 抽象を導入
- 既存の annotation/gc/skew/blast を “Track実装” に置き換え
- 設定DSL（schema_version付き）で track を列挙できるようにする

#### Phase 3（カスタムトラックの拡張）

- BED/WIG/VCF（小規模）を built-in として追加
- 大規模データ用に “server mode” を追加（`gbdraw serve`）

#### Phase 4（インタラクティブ）

- SVGグルーピング + data属性付与（選択/編集の前提）
- ラベル位置/track順序の保存・読み込み（configへ反映）
- （余力）scene graph 出力

---

### 提案：互換性（既存CLIを壊しにくくする）

- 既存の `gbdraw circular ...` / `gbdraw linear ...` は当面維持
- 内部で新APIへ変換する “adapter” を置く
- 新機能は `gbdraw render config.yml` のような config駆動を主導線にする

---

### 提案：開発体験（地味に効く）

- packaging: `pyproject.toml` 化（PEP517）を検討（長期的に配布が楽）
- 型: `mypy` または最小限の型整備（Track/Config周辺だけでも）
- lint/format: `ruff` + `black` の採用を検討
- テスト: “出力SVGがある程度不変”をスナップショット的にテスト（最小の回帰防止）

---

### 次に決めたい質問（ここが決まると設計が固まる）

- **ターゲットは「小規模微生物ゲノム中心」か、「巨大真核も見据える」か？**
  - 後者なら server mode とダウンサンプルは必須
- **インタラクティブはどこまで？**
  - まずは「ラベル移動/色変更/トラックON/OFF」程度か
  - それとも「トラック半径をドラッグで再配置」「クリックでフィルタ作成」までやるか
- **Web UIは Vue 継続でOK？**
  - 現状の `web/index.html` が完成度高いので継続が最短
- **Configは TOML 継続 or YAML/JSON に寄せる？**
  - Webとの相性はJSON/YAMLが良い（TOMLでも可だが編集体験は差が出る）

---

### 付録：この提案が刺さる理由（短く）

- “トラック”を第一級にすると、**新フォーマット対応 = 新トラック追加**に落とせる  
- “layout/render/ioの境界”があると、**Web/CLI/サーバ/Notebook**を同じコアで回せる  
- serverless と server mode を分けると、**Pyodideの制約**に引っ張られずに拡張できる  


