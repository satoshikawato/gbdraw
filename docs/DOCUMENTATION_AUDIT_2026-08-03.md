# ドキュメント監査と改修優先順位（2026-08-03）

## 結論

現状のドキュメントは、リンク、見出し、CLI help、最低限のサンプル実行についてはよく保守されている。一方、公開チュートリアルと掲載図の間に、より重大な問題が残っている。

- 入力ファイルから掲載図までの再現経路がない、または途中で途切れる。
- 文書のコードと掲載図が別のレシピから作られ、見た目の一致を検証していない。
- 実行できないチュートリアル手順と、現行実装に反するバージョン記述がある。
- 内部設計用の名前が一般向け文書に出ており、利用者が覚える必要のない概念を増やしている。

したがって、最初に直すべきなのは語句だけではない。P0の再現性・正確性・掲載図品質を直し、その後に「hosted」や `Web UX profile v1` を含む情報設計を整理するのが妥当である。

## 優先度の定義

| 優先度 | 判断基準 |
|---|---|
| P0 | 公開手順が実行できない、掲載図を再現できない、現行仕様と矛盾する、または完成例として明確な品質問題がある。次の公開前に修正する。 |
| P1 | 一部の利用者が再現で止まる、案内と実体が一致しない、または重複管理による陳腐化リスクが高い。P0に続けて修正する。 |
| P2 | 読みやすさ、発見性、用語、保守性の問題。P0/P1の正本を決めた後に整理する。 |

## 優先順位一覧

| 順位 | ID | 対象 | 問題 | 推奨する改修 |
|---:|---|---|---|---|
| 1 | P0-1 | [Release notes](./RELEASE_NOTES_0.14.0b0.md) | session v39を現行としているが、正本と実装はv40。 | v40導入内容と互換範囲を履歴として記録し、`current` を固定されたリリース時点の記述へ変える。 |
| 2 | P0-2 | [Tutorial 7](./TUTORIALS/7_Linear_Layout.md) | 宣言した作業ディレクトリから、未準備ファイルと解決不能な相対パスを使う。 | 全入力の取得・コピーを冒頭に集約し、クリーンな作業ディレクトリで全コマンドを通す。 |
| 3 | P0-3 | [Tutorial 9](./TUTORIALS/9_Feature_Visibility_Shapes.md) | 主要図の再現方法がGallery sessionの読み込みに依存し、元の`.gbk`からの完全なCLI/GUI手順がない。 | 実入力、完全なCLI、同値なGUI操作を本文の主経路にし、sessionは検算用の近道へ降格する。 |
| 4 | P0-4 | [Python API](./PYTHON_API.md) | MjeNMV図のラベル選択に説明がなく、ラベルと定量トラック/tickが競合する。コードと掲載PNGも別レシピ。 | quick startと完成例を分け、ラベル方針を説明し、実際のPythonコードから掲載図を生成する。 |
| 5 | P0-5 | [Tutorial 5: GFF3 + FASTA](./TUTORIALS/5_Table_Driven_Inputs.md#3-linear---records_table-for-gff3--fasta-rows) | 87 bpの合成例で、実ファイルでも完成図でもない。 | 既存のNCBI由来lambda fixtureを使い、records table、コマンド、結果図、検証値を掲載する。 |
| 6 | P0-6 | [Tutorial 5: similarity ring](./TUTORIALS/5_Table_Driven_Inputs.md#5-circular-blast-similarity-rings-with---conservation_table) | 23,709 HSP中5件だけが残り、掲載縮尺ではリングがほぼ空に見える。 | 閾値とHSP重複を再検討し、保持件数とreference上の被覆を明記してから図を再生成する。 |
| 7 | P1-1 | [Tutorial 2](./TUTORIALS/2_Comparative_Genomics.md) | 後半のmulti-record例が、冒頭で準備していないGenBank/BLAST入力を要求する。 | `examples/linear_multi_*.tsv` と関連入力を使う一つの準備手順に統合する。 |
| 8 | P1-2 | 全Tutorial | source checkout内の`tests/test_inputs`、暗黙のカレントディレクトリ、先行節の変数やファイルに依存する。 | 各TutorialにRequirements、Inputs、Working directory、Expected outputs、Validationを置く。 |
| 9 | P1-3 | 掲載図の生成契約 | Markdownのコマンド、Pythonコード、manifest recipe、掲載画像が別々に変わり得る。 | 一つの実行可能recipeからコマンドと図を作るか、両者の正規化SVG同等性をテストする。 |
| 10 | P1-4 | [Gallery](./GALLERY.md) | 「commands付き」と案内されるが、各図に入力、コマンド、versionがない。再生成に不足する入力もある。 | showcaseと呼び直すか、各図にcanonical recipe/session、入力入手先、versionを付ける。 |
| 11 | P1-5 | [Workflow guide](./WORKFLOW_GUIDE.md)ほか | `Web UX profile v1` は内部監査用の名前で、一般利用者に保証する内容が不明。`hosted web app` も反復される。 | 見出しを削除し、必要な「Web appの既定値/入力制限」だけを書く。製品名は`web app`に統一する。 |
| 12 | P1-6 | [Quickstart](./QUICKSTART.md) | 現在のNCBI annotationを取得する手順と、2025-07-17の固定fixtureから作った掲載図を同一視している。 | 固定fixtureとchecksumを配布するか、annotation更新により図が変わることを明記する。 |
| 13 | P2-1 | 情報構造 | README、Docs、Installation、Workflowにinterface選択とprivacy説明が重複する。内部計画文書も公開文書と同じ階層にある。 | 文書ごとの役割を決め、内部文書を`docs/internal/`等へ分離する。 |
| 14 | P2-2 | Reference | CLI helpと手書き説明、Typed APIとsession compatibility、release notesに同じ契約が重複する。 | 現行仕様の正本を一つにし、他文書は利用手順とリンクだけを持つ。 |
| 15 | P2-3 | Installation/Export | InstallationにPyPI経路がない一方、Exportは`pip install "gbdraw[export]"`を前提にする。 | PyPIを正式経路として説明するか、ExportをBioconda/source別の追加手順に直す。 |
| 16 | P2-4 | 細部 | 長いRecipes/FAQに目次がなく、Tutorial 3では`feature_specific`と`feature_specifc`が混在する。 | 大分類目次を追加し、出力prefixと掲載ファイル名を一致させる。 |

## P0の根拠と完了条件

### P0-1: release notesのsession versionを正す

[Release notesの145–149行](./RELEASE_NOTES_0.14.0b0.md#L145-L149)、[175–177行](./RELEASE_NOTES_0.14.0b0.md#L175-L177)、[192–194行](./RELEASE_NOTES_0.14.0b0.md#L192-L194)、[400–417行](./RELEASE_NOTES_0.14.0b0.md#L400-L417)は、session v39をcurrent writerとしている。これに対して、[Session compatibility](./SESSION_COMPATIBILITY.md#L10-L25)はwriter v40、accepted versions 27–33/39–40、typed bridge 31–33/39–40を正本としており、実装の[`CURRENT_SESSION_VERSION`](../gbdraw/session_io.py#L40-L46)も40である。さらに[手動browser acceptance script](../tests/run_losat_cache_browser_acceptance.py#L23-L31)が39を固定しているため、この更新漏れを検出できていない。

完了条件は次のとおり。

- v40で変わった内容を0.14.0b0の履歴へ追加する。
- release notesの「Current」は、将来変わらない「0.14.0b0 writes ...」へ置き換える。
- acceptance scriptが実装定数または共通fixtureから期待値を得る。
- release notes、compatibility table、実装定数の不一致を検出するcontract testを追加する。

### P0-2: Tutorial 7を空の作業ディレクトリで通す

[Tutorial 7の23–31行](./TUTORIALS/7_Linear_Layout.md#L23-L31)は`tutorial-7-work`を作って移動するが、その後の手順は次の入力を準備していない。

- 234–243行: `annotations.tsv`
- 257–266行: 移動後には解決しない`tests/test_inputs/...`
- 274–286行: `PemoMJNVA.gb`と`PeseMJNV.gb`

文書の途中でrepository rootへ戻すより、冒頭で必要入力を一覧化し、同じworking directoryへ取得・コピーする方が分かりやすい。完了条件は、空の一時ディレクトリから順番どおりに全code blockを実行し、記載した全出力が生成されることである。

### P0-3: Tutorial 9を`.gbk`起点にする

[Tutorial 9の31–65行](./TUTORIALS/9_Feature_Visibility_Shapes.md#L31-L65)は、Circular/Linear Gallery sessionをダウンロードし、GUIでshaft ratioだけを変える手順である。[67–92行](./TUTORIALS/9_Feature_Visibility_Shapes.md#L67-L92)のCLIは最小例で、掲載されているGallery品質の図を再現しない。掲載図のmanifestも、入力コマンドではなく[`SessionVariantRecipe`](../tools/reproduce_examples_manifest.py#L2095-L2141)から生成される。

Circularには、すでに[`HMMTDNA_ATSKEW_COMMAND`](../tools/prepare_interactive_gallery_assets.py#L106-L117)という完全なcanonical commandがある。これを起点に、`HmmtDNA.gbk`、qualifier-priority TSV、feature types、palette、window/step、labels、track slots、legend、font sizesを保ったまま、次の2項目だけを明示的に変える。

```text
--arrow_head_length_ratio auto
--arrow_shaft_width_ratio 0.75
```

Linearも同じ形式にする。5件のBGC GenBank、比較方法、色/qualifier table、labels、ruler、title、legendを含む完全なCLIを載せ、shaft ratioを`0.5`にする。GUI手順は「Upload」「Mode」「Features」「Labels」「Tracks/Comparisons」「Layout」「Arrow Geometry」「Export」の順に、CLIと同じ状態になる値を列挙する。session importは「完成状態との比較」または「短縮手順」として最後に置く。

完了条件は、sessionを使わず、記載された生入力とCLIだけで掲載図が生成でき、同じ生入力をGUIへ入れた場合も主要な視覚要素が一致することである。

### P0-4: Python APIのquick startと完成例を分離する

[Circular quick start](./PYTHON_API.md#L16-L83)は、4つの環境変数、2つのcolor table、label whitelist、4つのtrack slot、多数のconfig overrideを必要とする。これはAPIの最小導入ではなく、完成図のレシピである。

現在の[`python-api-label-whitelist.tsv`](../examples/python-api-label-whitelist.tsv)は、10件の異なるproduct名を一つの正規表現で選ぶが、選択基準を説明していない。掲載PNGでは複数の内側ラベルがGC content/skew帯や`300 kbp` tickと競合し、左側のleader lineも密集している。さらに、本文はPythonでSVGを保存する一方、掲載PNGは[`CliRecipe`](../tools/reproduce_examples_manifest.py#L1376-L1423)で作る「CLI-equivalent」図である。[実行テスト](../tests/test_api_library_usage.py#L54-L85)はSVGの存在だけを確認するため、Python出力と掲載図がずれても通る。

推奨する構成は次のとおり。

1. 最初の例は、1 input、1 `draw_circular()`、1 `save()`だけにする。
2. MjeNMVの完成例は別節に置き、必要ファイルの入手方法を示す。
3. label whitelistは、系統的な名称または明示した編集目的に沿って選ぶ。選択理由のない「一部のタンパク質」は残さない。
4. 完成例では外側ラベルを基本とし、GC/tickとの衝突、leader交差、可読倍率での密集を目視確認する。
5. 掲載SVG/PNGを、文書から抽出したPython recipeそのものから生成する。

完了条件は、独立したcode blockが追加の環境変数なしで動き、掲載図の生成元がそのPythonコードであり、白背景・通常表示倍率でラベルと定量トラックが衝突しないことである。

### P0-5: Tutorial 5のGFF3 + FASTAを実データへ置き換える

[現行手順](./TUTORIALS/5_Table_Driven_Inputs.md#L90-L133)は87 bp、CDS 1件の合成データをheredocで2組作る。ファイル形式のsmoke testとしては有効だが、公開Tutorialの完成例には情報量が足りず、結果図もない。

新しいデータを増やす必要はない。[`examples/gff3_lambda`](../examples/gff3_lambda/)には、RefSeq `NC_001416.1`から抽出した[`lambda_two_contigs.gff3`](../examples/gff3_lambda/lambda_two_contigs.gff3)と[`lambda_two_contigs.fna`](../examples/gff3_lambda/lambda_two_contigs.fna)がある。[既存の検証](../tests/test_api_library_usage.py#L87-L102)で、`lambda_left`/`lambda_right`、CDS 45件、両strand、translation保持まで確認されている。

同じGFF3/FASTA pairをrecords tableの2行で参照し、`record_id`で`lambda_left`と`lambda_right`を選ぶ。出典accession、抽出範囲、ファイルへの直接リンク、checksum、期待record/CDS数も書く。新しい掲載図をmanifestへ登録し、toy heredocは必要なら「format-only minimal example」としてreferenceへ移す。

完了条件は、実ファイルの取得またはcheckout内パスから、records tableと掲載図をそのまま再生成できることである。

### P0-6: similarity ringを教材として読める密度にする

[現行コマンド](./TUTORIALS/5_Table_Driven_Inputs.md#L163-L197)は`--identity 95 --alignment_length 1000`で、23,709 HSPから5件だけを残す。SVG内のpath数とalt textは一致しているため、生成バグではない。しかし、完成図では5本が短く、feature/GC/skewの間でほぼ空のリングに見える。

同じBLAST表に対する予備集計は次のとおり。

| 条件 | raw HSP数 |
|---|---:|
| identity >= 95、length >= 1000（現行） | 5 |
| identity >= 95、length >= 750 | 17 |
| identity >= 95、length >= 500 | 98 |
| identity >= 98、length >= 500 | 60 |

`95/750`は最初の候補になるが、閾値だけを下げて終わりにはしない。重複HSPを整理し、reference上のunion coverageと分布を見て、リング幅や背景track数も調整する。本文では`conservation`を推論しないことを維持し、利用者向けの見出しと凡例は可能な範囲で`BLAST similarity`に寄せる。

完了条件は、保持HSP数と被覆を本文に明記し、通常の文書幅でリングの位置と複数のspanを目視できることである。

## P1/P2の横断改修

### Tutorialを入力から検証まで一つの契約にする

[Tutorial index](./TUTORIALS/README.md#L3-L17)は番号付き文書を「complete set of command-line guides」と呼ぶが、現状はその契約を満たしていない。

- [Tutorial 2](./TUTORIALS/2_Comparative_Genomics.md#L19-L35)は最初の2 recordだけを準備し、後のmulti-record節は別の未準備入力へ切り替わる。
- [Tutorial 4](./TUTORIALS/4_Protein_Comparisons.md#L63-L91)は`tests/test_inputs`のstyle tableを使い、checkout外では省略するよう案内するため、掲載図と同じ結果にならない。
- [Tutorial 6](./TUTORIALS/6_Depth_Quantitative_Tracks.md#L17-L30)は前半のdepth例がsource checkout専用である。
- [Tutorial 7](./TUTORIALS/7_Linear_Layout.md)には前述の未準備入力がある。
- [Tutorial 8](./TUTORIALS/8_Interactive_SVG_Sessions.md#L42-L55)はschema v3と過去互換の説明がTutorialの操作手順に入り込んでいる。詳細はcompatibility referenceへ移せる。

全Tutorialで次の短い型を共有する。

1. Requirements
2. Working directory
3. Inputs（source、accession/version、checksum、取得コマンド）
4. CLI procedure
5. GUI procedure（対象機能にGUIがある場合）
6. Expected outputs
7. Validation（record数、feature数、hit数、主要な視覚要素）

placeholderを含む短いコマンドはRecipesへ、現在/過去のschema説明はReferenceへ、smoke inputはtestsへ分ける。

### recipe、文書、掲載図を同じsourceにする

現在のテストは、Markdown local target、Tutorial numbering/navigation、CLI option集合、figure inventoryをよく守っている。一方、Markdownのコマンドとmanifest recipeの一致、Python例と掲載図の一致、図の意味上の品質は検査しない。

追加すべきcontractは次のとおり。

- 公開recipeを一つの構造化データまたは実行可能scriptとして持ち、Markdown blockとfigure manifestがそれを参照する。
- クリーンな一時ディレクトリでTutorial code blockを順次実行する。
- SVGに期待するrecord/feature/HSP/track数をassertする。
- showcase図について、labelとtick/quantitative trackの重なり、leader crossing、空の比較ringを検出または人手確認する。
- Markdown heading fragmentの検査をCIへ追加する。

### Galleryを「眺めるページ」か「再現するページ」か決める

[Docs landing](./DOCS.md#L20-L32)と[README](../README.md)はGalleryをcommands付きとして案内するが、[Gallery本文](./GALLERY.md#L100-L104)には一般的なQuickstart/Tutorialへのリンクしかない。`M16-5_fugaku`、`Pandoravirus_salinus_forest`、4-record Escherichia-Shigella等には個別recipeがない。

figure inventoryには61件のdocs figureがあるが、10 recipeは必要入力を解決できない。公開Galleryに直接出る不足入力は`Shigella_dysenteriae.gbk`、`M16-5.gb`、`Pandoravirus_salinus.gb`で、archived recipeにも`NC_007205.gb`、`NC_005042.gb`、`NC_016510.gb`、`NZ_CP010822.gb`、`NC_000921.gb`、`NC_000962.gb`が不足している。追跡済み出力が存在することと、現在再生成できることは別である。

各Gallery entryに最低限、input accession/version、recipe/sessionへのリンク、gbdraw version、expected outputを付ける。付けない場合は、Docs/READMEの説明を「curated showcase」に直す。

### `Contents can’t be shown` の発生範囲と対応

初回の並行監査では、`Contents can’t be shown`相当のサービス停止が一度発生した。そのタスクに割り当てていた文書範囲は、[Tutorial 5の90–206行](./TUTORIALS/5_Table_Driven_Inputs.md#L90-L206)と[Tutorial 9の6–197行](./TUTORIALS/9_Feature_Visibility_Shapes.md#L6-L197)、および両者の掲載資産である。サービスは原因となった単一行・画像・語句を返さなかったため、これ以上狭い原因箇所は断定できない。後続監査をCLI、ファイル、再現手順に限定したところ、同じ停止は再発しなかった。

監査時点では、[Gallery 55行](./GALLERY.md#L55)が41,525,346-byteのSVGを、[Gallery 64行](./GALLERY.md#L64)が9,049,364-byteのSVGを`<img>`として直接埋め込んでいた。両方の`<img>`要素を削除し、full-size出力へのリンクだけを残した。Markdown previewはこの2ファイルを画像として読み込まなくなったため、同じ表示が再発する場合は別の箇所を調べる。

### 一般向け用語から内部契約名を外す

[`Web UX profile v1`](./WORKFLOW_GUIDE.md#L45-L62)は、[内部architecture audit](./CROSS_SURFACE_ARCHITECTURE_DOCUMENTATION_AUDIT.md#L403-L419)のP-01を閉じるために導入された名前である。実装上のdefault差を明示する目的は正しいが、利用者が`v1`を選択、保存、確認できるわけではない。表には実装のplot-title defaultも含まれず、公開version contractとしても不完全である。

この見出しとversion名は削除し、次の二つだけを平易に残す。

- Web appの既定値: separate strands、legend位置、multi-record表示。
- Web appの入力制限: Circularは1 GBFFまたは1 GFF3 + FASTA pair、複数source fileはCLI/Pythonを使う。

`hosted web app`は`web app`へ統一する。`hosted`または`self-hosted`を残すのは、Cloudflare配信、analytics、Galleryがlocal installに同梱されないことなど、配備場所そのものが意味を持つ文だけでよい。privacy説明の正本をInstallationに置き、README、Docs、FAQ、Workflowは一文とリンクにする。

### 現行仕様と履歴を分離する

- [CLI Reference](./CLI_Reference.md)は巨大なgenerated helpの直後に同じoptionを手書きで説明する。generated appendixとcurated semanticsを分ける。
- [Typed API](./TYPED_API.md)はsession/schemaのaccepted versionsを重複記載せず、[Session compatibility](./SESSION_COMPATIBILITY.md)へリンクする。
- FAQやRecipesの`now`、`new diagrams`、`older`といった差分起点の説明は現在形に直し、移行履歴はRelease notesに置く。
- [Recipes](./RECIPES.md)は「copy-paste commands」より「command templates」が正確である。外部BLAST+、入力ファイル、working directory等の前提を冒頭に置く。

## Quickstartと一般導線

[Quickstart](./QUICKSTART.md#L17-L31)はNCBI EUtilsから現在の`NC_000913.3` annotationを取得するが、掲載図の再生成は[`NC_000913.gbk`を固定された`tests/test_inputs/MG1655.gbk`へalias](../tests/test_reproduce_examples.py#L247-L269)している。accessionのsequence versionが同じでもannotationは更新され得る。再現性を優先するなら固定fixtureの公開URLとchecksumを使い、最新annotationを使うなら掲載図との差が出ることを明記する。

文書の役割は次のように整理する。

| 文書 | 役割 |
|---|---|
| README | 製品概要とDocsへの入口 |
| Docs landing | task別のrouting |
| Installation | web app、Bioconda、source/PyPI、privacyの正本 |
| Workflow guide | input/interface/outputを決める短いdecision guide |
| Tutorial | 実入力から検証済み出力までの連続手順 |
| Recipes | 前提を明示した短いcommand template |
| Reference | option、schema、compatibilityの現行契約 |
| Release notes | 特定versionで変わった履歴 |

完了済みのimplementation plan、capture plan、architecture audit、maintenance skillは一般向け文書と同じ階層に置かず、`docs/internal/`等へ移してaudienceとownerを明記する。これは内部用語がWorkflow guideへ漏れた問題の再発防止にもなる。

## 監査範囲と確認結果

対象は、現在の未commit working treeにある`README.md`、本報告書を除く`docs/`配下のMarkdown 35件（一般向け26件、内部/保守向け9件）、掲載図inventoryである。リリース済みtagではなく、2026-08-03時点の作業ツリーの評価である。

確認した文書群は次のとおり。

| 群 | 対象 | 主な判定 |
|---|---|---|
| Entry/general | README、ABOUT、DOCS、INSTALL、QUICKSTART、WORKFLOW_GUIDE、FAQ、GALLERY、RECIPES、EXPORT、GFF3_FASTA | 導線は成立。巨大SVG埋め込みは対応済み。用語重複、Gallery/Quickstart再現性、pip案内を改訂。 |
| Reference/API | CLI_Reference、PYTHON_API、TYPED_API、SESSION_COMPATIBILITY、SVG_SEMANTIC_HOOKS | CLI helpはcurrent。Python図、互換性重複、release note矛盾を改訂。 |
| Tutorials | index、1–9 | numbering/navigationは成立。2、5、7、9を優先し、全体をself-contained化。 |
| Release | RELEASE_NOTES_0.14.0b0 | package versionは一致するが、session v39記述を修正。 |
| Internal/maintenance | CROSS_SURFACE audit、IMPLEMENTATION_PLAN 4件、WEB_GALLERY plan/register 2件、docs/skills 2件 | 内容の公開対象境界を明示し、一般向け階層から分離。 |

次の自動確認は成功した。

```text
13 passed:
- documentation contracts
- tutorial documentation contracts
- public Markdown local targets
- public figure inventory coverage
- documented Python API execution
- documented GFF3/FASTA fixture semantics

CLI reference generated-help check: passed
Missing local heading fragments: 0
Key external entry pages and DOI redirects: HTTP 200 at audit time
```

この成功は、文書の機械的整合性が良いことを示す。一方、code blockがクリーンな環境で完結するか、掲載図がそのcode blockから作られたか、図が教材として読めるかは保証しない。本監査では指摘対象の図を目視し、全61 docs figureについてrecipe/input inventoryを確認したが、全図のpixel-level比較は実施していない。

## 推奨する実施順序

1. session v40矛盾とTutorial 7の実行不能を修正する。
2. Tutorial 9、Python API、Tutorial 5の掲載図とcanonical recipeを決め直す。
3. Tutorial 2を含む全Tutorialを空のworking directoryで実行し、入力契約を統一する。
4. recipeから文書と図を作るcontract test、意味上のSVG assertion、figure QA checklistを追加する。
5. Gallery/Quickstartの固定入力と再生成経路を埋める。
6. `Web UX profile v1`、`hosted web app`、重複するprivacy/interface説明を整理する。
7. Referenceの正本化、内部文書の分離、TOC・ファイル名等のpolishを行う。

掲載図を先に手作業で差し替えるべきではない。入力、canonical recipe、文書内コマンド、完了条件を先に確定し、そのrecipeから図を再生成して目視確認する順序が、同じずれを再発させにくい。
