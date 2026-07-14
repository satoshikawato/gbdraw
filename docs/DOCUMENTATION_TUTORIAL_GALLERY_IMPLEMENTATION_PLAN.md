# Documentation, Tutorial, and Gallery Improvement Plan

- 作成日: 2026-07-14
- 対象: `README.md`, `docs/`, `examples/`, `gbdraw/web/gallery/`, Gallery 関連の Web UI、生成・検証ツール
- 目的: ドキュメントと例を、正確で再現可能かつ目的から選びやすい教材へ改善する
- 前提: 2026-07-14 に実施した文章、教程、Gallery、例、用語、生成物の監査結果に基づく

## 1. 結論と実施方針

文章そのものに重大な AI-writing pattern はない。優先して直す対象は文体ではなく、Gallery コマンドの再現性、初心者向け導線、教程の重複、用語の正確さ、生成物の管理である。

実施順序は次のとおりとする。

1. 誤解を招く用語と再現不能なコマンド表示を直す。
2. Gallery の各例に、対象ユーザー、学習内容、難易度、所要時間、再現方法を付ける。
3. Beginner の最小例を追加し、既存の高度な例を case study として整理する。
4. Markdown tutorial を目的別に再編し、不足している API、GFF3 + FASTA、depth、export の教程を追加する。
5. 静的 Gallery と palette 生成物を整理し、すべての公開図に再生成経路を持たせる。

既存の高度な例は最初から削除しない。分類と説明を改善してから利用価値を再評価し、重複が確認できたものだけを archive または統合する。

## 2. 成功条件

計画完了時には、次の条件をすべて満たす。

1. Gallery のすべてのコマンドが `Runnable` または `Provenance` に分類されている。
2. `Runnable` コマンドが参照する入力と補助ファイルを、Gallery の Files/Downloads から取得できる。
3. `Provenance` コマンドには、そのままでは実行できない理由と、実際に試せる session または最小レシピが表示される。
4. WSSV の例が、完全な入力 bundle を公開しない限り「session-first の高度な case study」として表示される。
5. Gallery に少なくとも Circular と Linear の Beginner 例が 1 件ずつある。
6. 各 Gallery card に一文の用途説明、難易度、主な入力形式、所要時間または計算負荷の目安がある。
7. renderer 選択肢は、実装名と一致する `Dinucleotide content` / `Dinucleotide skew` を維持する。具体的な track は `GC content` / `AT skew` のように表示する。
8. `orthogroup` は CLI の互換値として維持しつつ、ユーザー向け説明では原則 `similarity group` と表示される。
9. Tutorial 2 と Tutorial 4 の protein-search 説明が重複せず、目的別 index から到達できる。
10. Python API、現実的な GFF3 + FASTA、BAM からの depth、出版用 export の案内がある。
11. 公開 Markdown が参照するローカル図の 100% が、再生成 manifest または明示的な手動管理リストに登録されている。
12. 公開ページから未参照の生成画像を 0 件にするか、保持理由を inventory に記録する。
13. Gallery の JSON、session、tutorial media、desktop/mobile 表示、offline packaging の既存検証が通る。

## 3. 非目標と互換性方針

- CLI の `--protein_blastp_mode orthogroup` は、この計画では改名しない。
- track renderer の内部名 `dinucleotide_content` と `dinucleotide_skew` は、session と CLI table の互換性のため維持する。
- `--conservation_blast` も、この計画では削除しない。説明文では `comparison` または `similarity` を優先し、将来 alias を追加する場合は別の互換性計画で扱う。
- Gallery の session、interactive SVG、既存 hash link は壊さない。
- `gbdraw/web/gallery/examples.json` を手作業の source of truth にしない。`tools/prepare_interactive_gallery_assets.py` の定義を更新して再生成する。
- 見栄えだけを理由に新しい例を増やさない。各例には固有の学習目的を要求する。
- palette SVG を削除する前に、公開リンク、README、packaging、生成ツールからの参照を移行する。

## 4. 監査時点の基準値

以下は実装開始前に再計測する。値が変わっていた場合は、計画書の数値を更新してから作業する。

| 項目 | 監査値 | 問題 |
|---|---:|---|
| Interactive Gallery | 7 例 | Circular 3、Linear 4 |
| Beginner tutorial | 0 例 | 全件が Intermediate または Advanced |
| Linear Gallery | 4 例 | すべて protein-search 系 |
| 実行に必要なファイルが表示されない command | 4 例 | Human mtDNA、BGC、WSSV、Majanivirus |
| Gallery tutorial media | 80 files | strict validation は通過 |
| Tutorial 2 | 310 行、約 1,935 語 | 扱う課題が多すぎる |
| Palette | 55 palettes、110 SVG | 同じ入力の Circular/Linear を大量保持 |
| `examples/` | 約 677 MB | 生成可能な SVG の比率が高い |
| 公開図の再生成 manifest 未登録 | 30 files | 新しい tutorial 図が主 |
| 公開 Markdown から未参照の `examples/` media | 58 files、約 120.4 MiB | 削除・archive・公開の判断が必要 |

監査時には、公開 Markdown の相対リンク欠損は 0 件、Gallery tutorial media の strict check は 80/80 件通過している。

## 5. 用語の決定

### 5.1 ユーザー向け表示

| 現在の語 | 採用する語 | 実装上の扱い |
|---|---|---|
| Dinucleotide content | Dinucleotide content | renderer 選択肢では実装名を維持し、`nt` が `GC`/`AT` の実 track は具体名を表示する |
| Dinucleotide skew | Dinucleotide skew | renderer 選択肢では実装名を維持し、実 track は GC skew、AT skew のように表示する |
| Orthogroups | Similarity groups | CLI の literal value は補足として表示する |
| Nucleotide-similarity rings | BLAST/LOSAT comparison rings | 初出で `nucleotide-similarity rings` を併記する |
| Conservation rings | Comparison rings または Similarity rings | 生の HSP span を進化的 conservation と呼ばない |
| Collinear blocks | Collinear blocks | 標準的なので維持する |
| Publication-quality | Publication-quality | 分野で一般的なので維持する |

根拠として、BLAST ring/track の一般的な表現は [Proksee の BLAST track 説明](https://proksee.ca/posts/47)、orthogroup の定義は [OrthoFinder の公式説明](https://orthofinder.github.io/OrthoFinder/tutorials/understand-orthology/)、collinear block と anchor は [syntenet manual](https://www.bioconductor.org/packages/release/bioc/manuals/syntenet/man/syntenet.pdf) を参照する。

### 5.2 文章ルール

- 「何が描かれているか」と「そこから解釈できること」を分ける。
- BLAST/LOSAT の生 HSP を conservation、orthology、shared function の証拠として断定しない。
- 操作名は Web UI の現在の literal label と一致させる。
- `exact`, `same`, `match the Gallery session` を繰り返す代わりに、設定の目的と変化を説明する。
- 各 tutorial の冒頭に「この教程で答える問い」と「完了後にできること」を置く。
- AI-writing detector を満たすためだけの言い換えは行わない。直接的で技術的な現行 voice を維持する。

## 6. 実装段階

### Phase 0: 非回帰基準と inventory を固定する

目的は、修正前後の差を意図した変更として説明できる状態にすることである。

実施内容:

1. `tools/prepare_interactive_gallery_assets.py` の `EXAMPLES` と生成済み `gbdraw/web/gallery/examples.json` の差を確認する。
2. 7 tutorial の `difficulty`, `estimatedTime`, `requirements`, `downloads` を一覧化する。
3. 各 Gallery command からファイル引数を抽出し、次の区分を記録する。
   - Gallery から取得可能
   - repository checkout のみに存在
   - session 内に埋め込まれている
   - 生成方法のみ記載されている
   - 入手不能または provenance 不明
4. 公開 Markdown の画像参照、`examples/` の画像、`build_figure_specs()` の output path を比較する inventory check を追加する。
5. 現在の Gallery desktop/mobile screenshot と tutorial contact sheet を基準として保存する。新しい基準画像を増やす場合は、既存の screenshot maintenance 手順に従う。

完了条件:

- 7 例すべてに command classification の予定値がある。
- 公開図の manifest 未登録と未参照画像の一覧を機械的に再生成できる。
- 既存の strict media check と Gallery browser test が改修前に通る。

### Phase 1: 正確さと再現性を直す

#### 1-A. GC/AT content と skew の表示

対象:

- `gbdraw/web/js/app/circular-track-slots.js`
- `gbdraw/web/js/app/linear-track-slots.js`
- `gbdraw/web/gallery/tutorials/HmmtDNA_ATskew.json`
- `tests/web/gallery-tutorial.playwright.spec.js`
- `tests/test_web_packaging.py`
- 関連する README、tutorial、alt text

実施内容:

1. renderer 選択肢の一般表示は `Dinucleotide content` / `Dinucleotide skew` を維持する。
2. `nt=GC` と `nt=AT` の表示箇所では `GC content`, `AT content`, `GC skew`, `AT skew` を優先する。
3. 内部値 `dinucleotide_content` / `dinucleotide_skew`、session JSON、CLI track-slot syntax は変更しない。
4. tutorial の本文、表、alt text、browser test の期待値を同時に更新する。

完了条件:

- UI と公開教程では renderer 名として `Dinucleotide content` を使い、具体的な GC/AT track は個別名で表示する。
- 既存 session が同じ track renderer を復元する。

#### 1-B. Gallery command の分類

`GallerySessionExample` に必要最小限の metadata を追加する。

```text
command_kind: runnable | provenance
command_note: string
```

必要であれば、取得可能な補助ファイルを tutorial の `downloads` から Files panel に集約する。別の command schema や汎用 workflow DSL は作らない。

例ごとの方針:

| Example | 分類・対応 |
|---|---|
| Human mtDNA AT skew | qualifier TSV を配布し `Runnable` にする |
| BGC | default/specific/qualifier TSV を配布し `Runnable` にする |
| WSSV | 入力 bundle を公開しない限り `Provenance`; session-first と表示する |
| Majanivirus | 色 TSV を配布できれば `Runnable`; できなければ `Provenance` |
| Hepatoplasmataceae collinear | `Runnable`; 重複した `--losatp_threads 32` を削除する |
| Hepatoplasmataceae similarity groups | `Runnable`; user-facing label を更新する |
| Vibrio multi-record | `Runnable`; GenBank source の取得方法を Files に表示する |

全コマンドの format spelling は CLI reference で推奨する `interactive_svg` に統一する。旧 `interactive-svg` は互換 alias として受理する。

完了条件:

- `Runnable` の Copy ボタン付近に不足ファイルがない。
- `Provenance` は Copy してすぐ実行できるように見せない。
- command kind と note の欠損を `tests/test_web_packaging.py` で検出する。

#### 1-C. WSSV tutorial の位置付け

実施内容:

1. タイトルまたは summary で `Advanced session-based case study` と明記する。
2. 第一の経路を「bundled session を読み込んで結果と設定を調べる」にする。
3. 20 FASTA の手動 rebuild は、完全な prepared input bundle と checksum を公開できた場合だけ exact reproduction と呼ぶ。
4. SRA 由来入力には assembly/preparation が必要で、生 run は直接の代替ではないと維持する。
5. accession 不明の Shantou2019 は provenance gap と明示する。推測した accession は追加しない。
6. 所要時間を `input preparation`, `search`, `render/interaction` に分ける。

完了条件:

- 公開情報だけで不可能な手順を「manual exact reproduction」と呼ばない。
- session の inspect path は、外部入力を取得しなくても完了できる。

### Phase 2: Gallery を目的から選べる構成にする

#### 2-A. Card metadata と並び順

`GallerySessionExample` に、実際に card と詳細画面で使う metadata だけを追加する。

```text
description
difficulty: Beginner | Intermediate | Advanced
workflow
input_summary
estimated_time
display_order
```

`description` field はすでに generator と Gallery UI が扱えるため、まず全例を埋める。`examples.json` の file size 順 sort はやめ、`display_order` を使う。file size は Files panel に残し、card では difficulty と workflow を優先する。

各 description は、機能の列挙ではなくユーザーの目的を書く。

例:

- 悪い例: `Circular, LOSAT, comparison rings.`
- 良い例: `Compare one viral reference with many assemblies as concentric nucleotide-similarity rings.`

完了条件:

- 7 例すべての card だけを読んで、用途と難易度を区別できる。
- 最初に Beginner、次に代表的な Intermediate、最後に Advanced case study が並ぶ。

#### 2-B. 最初の Beginner 例

一度に多数の新規例を作らず、まず次の 2 例を追加する。

1. `HmmtDNA_basic_circular`
   - 入力: 小さい GenBank 1 file
   - 学習内容: upload、Circular mode、feature tracks、GC content/skew、SVG export
   - 計算: 検索なし
   - 目標時間: 5 分以内
2. `lambda_basic_linear`
   - 入力: 小さい GenBank 1 file
   - 学習内容: Linear mode、strand separation、labels、region/scale、SVG export
   - 計算: 検索なし
   - 目標時間: 5 分以内

既存の Human mtDNA AT skew tutorial は、basic circular から続く Intermediate tutorial として関連リンクを付ける。

追加後の判断:

- Beginner 2 件が機能した後に、2-genome pairwise comparison を追加する。
- 既存の Hepatoplasmataceae 2 件は、同じ入力から similarity links と collinear blocks の違いを比較する paired case study として表示する。
- Majanivirus と WSSV は密度の高い Advanced 例として後方に置く。

完了条件:

- 新規ユーザーが検索 runtime や 5～20 input files を準備せず、最初の図を作れる。
- Beginner tutorial の各 step が既存の screenshot validation に登録される。

#### 2-C. Workflow filter

例が 10 件未満の間は、複雑な検索 UI を追加しない。まず tag を次の統制語彙にそろえる。

- Layout: `Circular`, `Linear`
- Input: `GenBank`, `GFF3 + FASTA`, `Multi-record`, `Depth`
- Comparison: `BLAST`, `LOSAT`, `Similarity groups`, `Collinear`
- Level: `Beginner`, `Intermediate`, `Advanced`

10 件以上になった時点で、既存 tag を使う client-side filter を追加する。自由入力検索は、それでも選択が困難な場合だけ検討する。

### Phase 3: Markdown tutorial の情報設計を直す

#### 3-A. Index を目的別にする

`docs/TUTORIALS/TUTORIALS.md` を番号順の講座一覧ではなく、次の入口へ再編する。

1. First figure
2. Compare genomes
3. Labels, colors, and feature visibility
4. Quantitative tracks
5. Layout and record selection
6. Interactive output and reproducibility
7. Automation with Python

既存 filename は外部リンク保護のため直ちに変更しない。分割後に不要になった旧ページは、短い案内と新ページへのリンクを残す。

#### 3-B. Tutorial 2 と 4 の責務を分ける

Tutorial 2 に残す内容:

- 2 record の precomputed BLAST
- query/subject と input order
- HSP ribbon の意味と限界
- multi-record comparison への拡張
- circular BLAST/LOSAT comparison rings へのリンクまたは独立節

Tutorial 4 に集約する内容:

- CDS translation と runtime selection
- `pairwise`, `orthogroup`, `collinear` の違い
- similarity group が phylogeny-based orthogroup ではないこと
- collinear anchors、blocks、thresholds、evidence scope

移動する内容:

- record selector、crop、reverse complement は layout/input selection tutorial へ移す。
- TSV manifest の繰り返しは Tutorial 5 に一本化する。
- circular LOSATN/TLOSATX が長い場合は独立 tutorial にする。

完了条件:

- 同じ protein mode の説明と command が複数 tutorial に重複しない。
- Tutorial 2 を終えるために、protein search runtime の理解を要求しない。
- Tutorial 4 だけで 3 protein modes を比較できる。

#### 3-C. Quickstart を短く保つ

実施内容:

1. Circular の最初の図を主目的として維持する。
2. `wget` に加え、NCBI の直接リンクまたは `curl` 例を示す。OS 別手順を大量には増やさない。
3. Linear region selector は短い Next step に移し、Basic Linear Gallery tutorial へリンクする。
4. 目標時間を 10 分以内とする。

### Phase 4: 不足している高価値ドキュメントを追加する

#### 4-A. Workflow chooser

新規候補: `docs/WORKFLOW_GUIDE.md`

次の判断を短い表または decision tree で説明する。

- Circular と Linear
- GenBank と GFF3 + FASTA
- precomputed BLAST と built-in protein search
- pairwise links、similarity groups、collinear blocks
- CLI、Web app、Python API
- static SVG/PDF/PNG と interactive SVG/session

#### 4-B. Python API guide

新規候補: `docs/PYTHON_API.md`

`gbdraw.api` だけを公式 entry point として使用し、内部 module の import を教程に載せない。

最低限含める内容:

1. GenBank を読み込む。
2. `build_circular_diagram` または `build_linear_diagram` を呼ぶ。
3. options を指定する。
4. SVG または bytes/file に出力する。
5. error handling と API stability の範囲を説明する。

教程中の code block は smoke test または doctest 相当で実行し、API drift を検出する。

#### 4-C. 実用的な GFF3 + FASTA tutorial

toy の手書き GFF だけでなく、現実的な annotation output を使う。

含める内容:

- GFF3 と FASTA の record ID 一致
- CDS phase、strand、translation の扱い
- 複数 contig
- annotation に translation がない場合の protein search 制約
- よくある validation error と直し方

データは小さく、license と provenance を明示し、repository または固定 URL から取得可能にする。

#### 4-D. Depth tutorial の入力生成

既存の ready-made depth TSV に加え、BAM からの最小経路を追加する。

```bash
samtools depth -aa input.bam > sample.depth.tsv
```

説明対象:

- `-aa` と zero-coverage positions
- contig/record ID の一致
- missing positions
- log/linear scale
- window/step と raw per-base depth の違い
- 複数 record と複数 sample

大きな BAM は repository に追加しない。小さい fixture または既存 TSV の生成 provenance を示す。

#### 4-E. Publication/export guide

新規候補: `docs/EXPORT.md`

含める内容:

- SVG、PDF、PNG、EPS/PS の選択
- vector と raster の使い分け
- PNG の DPI と最終掲載サイズ
- CairoSVG など optional dependency
- font と italic markup の注意
- vector editor で修正する場合に保持すべき metadata
- interactive SVG と static export の違い
- color-vision accessibility と grayscale check

### Phase 5: Static Gallery と palette を整理する

#### 5-A. `docs/GALLERY.md`

Gallery を次の二層に分ける。

1. Capabilities gallery
   - 重複の少ない代表図
   - 一文の用途説明
   - 関連 tutorial へのリンク
2. Reproducible recipes
   - 入力取得方法
   - 補助 TSV
   - 完全な command
   - expected output

`./in_gbk/` のような repository に存在しないローカルパスを、そのまま runnable recipe として載せない。歴史的な図を残す場合は Archive/Provenance と明示する。

#### 5-B. Palette page

`examples/color_palette_examples.md` は、55 palette ごとに 2 個の大型 SVG を並べる方式をやめ、次の構成へ移行する。

1. palette 名、feature 色、簡易 preview を一覧表示する contact sheet
2. palette ごとの hex code を展開表示する compact table
3. representative な Circular/Linear SVG を数件だけ保持する
4. contrast、color-vision accessibility、print suitability の注記を追加する

contact sheet と table は `gbdraw/data/color_palettes.toml` から生成する。Markdown と palette data の二重管理は避ける。

削除条件:

- 既存 110 SVG への公開リンクが contact sheet または代表例へ移行済みである。
- `tools/reproduce_examples.py` が新しい成果物を再生成できる。
- packaging test と docs link check が通る。
- `docs/CODEBASE_REDUCTION_PLAN.md` の palette 削減作業と重複実装しない。

#### 5-C. 公開図の manifest coverage

`tools/reproduce_examples_manifest.py` を公開図の source of truth とする。

実施内容:

1. 現在未登録の tutorial 図を追加する。
2. 意図的に手動管理する画像には理由と更新手順を登録する。
3. 公開 Markdown のローカル画像参照が manifest または手動管理リストに含まれることを test する。
4. manifest にあるが公開ページから未参照の図を検出する。
5. 生成済み `.actual.svg` や wheel は追跡対象へ追加しない。

### Phase 6: 継続的な検証を固定する

新しい独立 validator を多数作らず、既存の生成・packaging・screenshot test に小さい invariant を追加する。

#### 6-A. Static checks

- Gallery entry に description、difficulty、workflow、command kind がある。
- `Runnable` command の repository-managed dependency が存在する。
- Gallery tutorial の JSON と media が一致する。
- 公開 Markdown の相対 link target が存在する。
- 公開ローカル画像が再生成 inventory に登録されている。
- user-facing prose に禁止した旧表示語が残っていない。

#### 6-B. Browser checks

- card に description と difficulty が表示される。
- Beginner が最初に並ぶ。
- Runnable/Provenance の表示が区別できる。
- Copy button は現在選択中の command をコピーする。
- desktop/mobile で card、tabs、Files、Tutorial が読める。
- existing feature/match popup と interactive SVG は非回帰である。

#### 6-C. 共通検証コマンド

各 PR では変更範囲に応じて次を実行する。

```bash
python tools/capture_gallery_tutorial_screenshots.py --all --check --strict
python tools/prepare_interactive_gallery_assets.py
python -m pytest tests/test_web_packaging.py -q
python -m pytest tests/test_reproduce_examples.py -q
python -m pytest tests/ -v -m "not slow"
ruff check gbdraw/ --select=E,F,W --ignore=E501,W503
```

Gallery UI または tutorial media を変更した PR では追加で実行する。

```bash
node -e "console.log(require.resolve('@playwright/test'))"
playwright test tests/web/gallery-tutorial.playwright.spec.js
```

Node の `@playwright/test` がない環境では、Python Playwright で対象導線を確認する。Chromium が sandbox error で失敗した場合は、AGENTS.md に従って同じ check を sandbox escalation 付きで再実行する。

## 7. PR 分割

| PR | 内容 | 主な対象 | 完了条件 |
|---|---|---|---|
| 1 | 用語と明白な command bug | track-slot labels、Hmmt tutorial、重複 thread option | UI label、JSON、browser test が一致 |
| 2 | Command classification と WSSV 再定義 | Gallery generator、Gallery UI、WSSV tutorial | 全 7 例が Runnable/Provenance に分類 |
| 3 | Gallery metadata と順序 | generator、gallery.js/css/html、packaging test | 用途・難易度・workflow が card に出る |
| 4 | Beginner Circular/Linear | sessions、SVG、tutorial JSON、media | 5 分以内の検索なし tutorial が 2 件 |
| 5 | Tutorial 2/4 と index 再編 | `docs/TUTORIALS/`、Quickstart | 重複を除き、旧 link を維持 |
| 6 | Workflow/API/GFF/depth/export | 新規 docs、small fixtures、tests | 各教程に実行可能な最小例がある |
| 7 | Static Gallery と palette | `docs/GALLERY.md`、palette page、reproducer | 代表例へ集約し、公開 link が通る |
| 8 | Figure inventory gate | reproduce manifest、tests | 公開図の coverage 100%、未参照 0 または理由あり |

PR 1～3 は同一 release 内で先に完了させる。PR 4 以降は独立して進められるが、Gallery metadata schema を別々に増やさないため PR 3 完了後に着手する。

## 8. 例を採用する基準

新しい tutorial または Gallery example は、次の質問にすべて答えられる場合だけ追加する。

1. どのユーザーの、どの作業を助けるか。
2. 既存の例にはない学習内容は何か。
3. 入力を合法かつ安定して取得できるか。
4. ゼロから再現できるか。できない場合、なぜ case study として残す価値があるか。
5. 典型的な所要時間と外部 runtime は何か。
6. 図を見て何を確認すべきか。
7. 結果から解釈してはいけないことは何か。
8. 入力、command/session、expected output、生成 version を追跡できるか。
9. automated check または明示的な手動検証手順があるか。

単に別の生物種、palette、label 配置を見せるだけで既存例と学習内容が同じ場合は、新規 tutorial にせず Gallery thumbnail、palette preview、または recipe variation として扱う。

## 9. リスクと対策

### Generated file の手修正

リスク: `examples.json` や Gallery asset を直接直すと、次回生成時に消える。

対策: generator、session、tutorial JSON の source を先に直し、生成後の diff を review する。

### 互換用語とユーザー向け用語の混同

リスク: `orthogroup` や `dinucleotide_content` を一括置換すると session/CLI compatibility を壊す。

対策: UI label、説明文、internal identifier を分け、identifier の変更はこの計画の対象外とする。

### Beginner example の増加による repository 肥大化

リスク: session、source SVG、interactive SVG、thumbnail、operation media が例ごとに増える。

対策: 小さい入力を使い、既存 fixture を再利用し、generic screenshot を増やさない。追加例による byte 増加を PR に記載する。

### 外部 source の変化

リスク: NCBI や外部 repository の URL、record version、annotation が変わる。

対策: accession version、取得日、checksum を記録する。教程を特定の未固定 `latest` URL に依存させない。

### Screenshot の過剰保守

リスク: UI 文言変更のたびに多数の media が無効になる。

対策: screenshot は例固有の状態や視覚結果に限定し、generic controls は本文または表で説明する。

### Palette 削減による既存 link 切れ

リスク: GitHub、README、外部記事から個別 SVG が参照されている可能性がある。

対策: tracked file を削除する PR では repository 内参照に加えて GitHub 上の既知 link を確認し、必要なら archive path または release asset を用意する。

## 10. 完了チェックリスト

- [ ] 用語表のユーザー向け名称が Web UI、tutorial、Gallery、alt text で一致している。
- [ ] すべての Gallery command に Runnable/Provenance が表示される。
- [ ] Runnable command の全依存ファイルを取得できる。
- [ ] WSSV が session-first case study として正確に説明される。
- [ ] Circular/Linear の Beginner tutorial が各 1 件以上ある。
- [ ] Gallery card に用途、難易度、workflow、入力、時間の情報がある。
- [ ] Tutorial index が目的別になっている。
- [ ] Tutorial 2/4 の重複が除かれている。
- [ ] Python API guide の code example が自動実行される。
- [ ] GFF3 + FASTA tutorial が実際の multi-feature input を使う。
- [ ] Depth tutorial に BAM からの作成経路がある。
- [ ] Export guide に形式、DPI、font、accessibility がある。
- [ ] Static Gallery が showcase と reproducible recipe を区別する。
- [ ] Palette page が palette data から生成される。
- [ ] 公開ローカル画像の再生成 coverage が 100% である。
- [ ] 未参照画像がないか、保持理由が記録されている。
- [ ] Markdown link、Gallery media、packaging、browser、fast test が通る。
- [ ] 生成された wheel、`dist/`、`gbdraw.egg-info/` を手作業で変更していない。
