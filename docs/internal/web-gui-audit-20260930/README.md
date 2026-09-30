# Web GUI バグ監査 — dev 4c89bab1（2026-09-30）

対象は `origin/dev` の `4c89bab1`（PR #648 のマージ時点）。この commit の CI（Tests、Gallery publication）は成功している。
コードの変更、commit、issue 作成は行っていない。

**固有の不具合は 65 件（P1 10、P2 33、P3 22）。** 64 件は dev を実際に動かして再現し、PV-12 の 1 件だけがコードからの判断である。
既知項目（F01–F15、G01–G10、#543–#619）の退行は見つからなかった。#599 のドラッグ位置の保持も確認した。

## 方法と表記

- 7 領域を並列に監査した。コードから仮説を立て、Chromium（Playwright 1.61）で dev を動かして再現した。必要に応じて dev の CLI（`PYTHONPATH=<dev> python -m gbdraw.cli`）と比べた。
- P1 の IN-02、PV-08、CO-05 は CLI または pandas で、FE-01、TR-01 は spec の再実行で、担当とは別に再確認した。
- 重大度は次のとおり。
  - **P1**: 誤った科学的出力、編集の消失、有効な入力で Generate できない、CLI も同じく誤る。
  - **P2**: 機能が壊れている、または Generate まで古い状態が残る。回避策はある。
  - **P3**: 軽微、端のケース、表示やアクセシビリティ。
- ID の接頭辞: IN 入力、SE Session/History、GE Generate/オプション、CO 比較、FE Feature 編集、PV 凡例/Preview/出力、TR トラック/画面、X 領域横断。
- JS のパスは `gbdraw/web/js/` からの相対パスで書く。それ以外はリポジトリ直下からのパス。行番号は `4c89bab1` のもの。

## P1（10 件）

| ID | 症状 | 原因 |
|---|---|---|
| IN-01（=GE-01） | Circular で Generate 後に Species、Strain、Plot title、定義フォントなどを編集すると、Result の定義行がファイル全体から作り直される。Region の長さと GC%、Record label が失われる（例: `8,001 bp / 44.41% GC` → `16,569 bp / 44.36% GC`）。Generate していない別ファイルの accession や record ID が書き込まれ、grid では定義が入れ替わる。書き換わった Result がそのまま出力と Session 保存に使われる | `app/results.js:234-368` が draft の `files.c_gb` を `app/python-helpers.js:1693-1827`（regenerate_definition_svgs）へ渡す。この処理は record 選択、region、逆相補、ラベル、並び順を無視し、`record_index` で照合して `:354-361` で Result に書き戻す。起点は `app/watchers.js:574-585` |
| IN-02 | Prokka 形式（`ACCESSION`/`VERSION` 行が空）の GenBank で record ID が `KEYWORDS` になる。単一 record の Circular は RECORD_SELECTION/NO_MATCH で Generate できず、複数 record では出力名が `KEYWORDS.svg` になる。CLI は `contig_1.svg` | `app/record-discovery.js:137-138` の `/^VERSION\s+(\S+)/m` の `\s+` が改行をまたぐ。ACCESSION も同様に `VERSION` を拾う。`0ce21a32`（2026-08-10）から存在し、main にも同じ正規表現がある（main での挙動は未確認） |
| IN-03 | GFF3+FASTA で、GFF に feature のない配列が FASTA にあると、Circular（All records）と Linear の自動展開が NO_MATCH で失敗する。CLI は描画できる。回避策は、注釈のある record を選ぶか、その行を削除すること | JS の高速経路（`app/record-discovery.js:213-227`）は FASTA ヘッダをすべて列挙する。一方 Python（`gbdraw/io/genome.py:105-135`）は feature のない record を捨てる |
| PV-08 | 既定設定では、やや長い `/organism`（例 "Salmonella enterica subsp. enterica serovar Typhimurium"、"Mycobacterium tuberculosis variant bovis BCG"）を持つ record を Generate できない。TRACK_LAYOUT/CANNOT_FIT（`ticks`）になる。CLI の既定では成功するが、Web の既定と同じ `--multi_record_canvas --track_type tuckin --gc --skew --separate_strands` では失敗する。エラーは track を指し、原因である定義文字列に触れない | 推定: grid 配置が定義文字列のために中央の領域を確保し、ticks スロットを圧迫している。Web は record が 1 つでも grouping grid を使う |
| FE-01 | ラベル文字とラベル表示の編集が黙って消える。batch で Result を切り替えたとき、feature を非表示にして Generate したとき、Circular の record 選択を往復したときに、override が `{}` になる。History には記録されず、その後の Generate と Session 保存で編集が失われる | `app/feature-editor/label-actions.js:689-702` の syncLabelEditor は、表示中 SVG の feature-ID の組が変わると `clearOverrides()` を呼ぶ（buildContextKey `:342`、hasFeatureScopedOverrideInSvg `:353`） |
| TR-01 | Multi-Record Canvas（Web の既定で ON）では、凡例が custom track slot を反映しない。スロットの色、legend label、追加した skew 行が凡例に出ず、リングは指定色なのに凡例は既定色のまま。CLI の `--multi_record_canvas` も同じ | `gbdraw/api/diagram.py:3114-3128` は `prepare_legend_table(show_gc/show_skew)` だけを使い、`_sync_legend_table_for_circular_slots`（`gbdraw/diagrams/circular/assemble.py:1280-1400`）を呼ばない。呼んでいるのは単一 record の経路（`:2766`）だけ |
| CO-02 | LOSATP の実行後に Feature visibility を変えて Generate すると、古い protein 結果が再利用され、非表示にした CDS へのリンクが残る（19 ribbon / 18 group。新規に実行すると 18 / 17） | `app/run-analysis.js:468-541` の canReuseResolvedProteinArtifacts が visibility table（featureVisibilityCacheKey、`:2352-2357`）を比較しない |
| CO-03 | Feature popup の「Rotate record using this feature」を Orient forward ON で使うと、LOSATN の ribbon が相同性のない領域に描かれる。新たに Generate しても同じ | 向きは `recordDisplayDrafts[].reverseComplementOverride`（`app/record-display-options.js:205-233`）に保持される。しかし LOSAT の表示座標への変換は `seq.region_reverse` しか見ない（`app/run-analysis.js:3177-3208`）。`services/session-request.js:5009-5017` も比較を投影し直さない（OIC-021 違反） |
| CO-04 | 1 ファイルにまとめた複数 record は 1 つの batch として検索される。そのため同じ record でも、ファイルの分け方でリンクが変わる。Max target seqs=1 では別ファイルなら 19、1 ファイルなら 0（CLI は 19 / 19）。上限なしでも E-value が変わるので、閾値付近の判定が変わりうる | `app/linear-sources.js:151-212` の prepareLosatSourceBatches。自己ヒットが max-target-seqs を埋め、DB サイズも合算される |
| CO-05 | 13 列以上の outfmt 6（例 `-outfmt "6 std qlen slen"`）が警告なしに誤読される。14 列では ribbon が 0 本。13 列では列が 1 つずれ、identity 200% などの ribbon が描かれる。CLI も同じ | `pd.read_csv(names=<12 列>)` で余った先頭の列が index になる（`gbdraw/io/comparisons.py:73-78`、`gbdraw/session_request_codec.py:3330-3335`、`gbdraw/api/record_planning.py:1426-1431`、`gbdraw/analysis/conservation.py:175-180`）。`conservation.py:147-149` の対処は DataFrame 入力だけに効く（設計段階で訂正） |

## P2（33 件）

### 入力と Session

| ID | 症状 | 原因 |
|---|---|---|
| IN-04 | Linear で region を指定したファイルを複数 record のファイルに置き換えると、古い region と定義が残り、record が展開されない。Generate は汎用エラーになる | `app/app-setup.js:4103-4132`。record が 2 つ以上のとき（`:4121`）だけ region を解除し、定義は常にコピーする。region があると展開を飛ばす（`:1049`） |
| IN-05 | Circular の Multi-Record Canvas の並び順が、Linear へ切り替えて戻るとリセットされる | `app/run-analysis.js:1735-1738, 1748, 1777-1780`（入力の有無の判定にモードが含まれる）、`app/watchers.js:586-603` |
| IN-06 | 複数 record を含む Linear の GenBank で、全 record の定義が 1 番目の record の organism になる | `app/app-setup.js:1064-1081`、`app/linear-sources.js:247-255, 271-276` |
| IN-08 | 最初の Generate の前に、いずれかのモードの入力が不完全だと、Save Session が UNKNOWN で失敗する | `services/config.js:4057-4058, 4078-4102`、`services/session-request.js:641` |
| SE-05 | v39 など旧形式の Session は読み込めるが、Save が UNKNOWN で失敗する。本来の案内は "Generate again before using Save Session…" | `services/config.js:224, 231, 4018-4019`（旧形式の読み込みでは feature catalog が null のまま） |
| SE-06 | CLI の `--session_output` で作った Linear+BLAST の Session が読み取り専用扱いになり、比較付きで Generate できない。Inherit を選んでも失敗する | CLI は `cli-seq-N` を書くが Web は `record-N` を使う（`gbdraw/cli_utils/session.py:1321`）。無効な protein 比較の項目も拒否される（`services/imported-comparison-intent.js:191, 437-447`） |
| SE-07 | CLI の Session を読み込むと凡例の位置が失われ、最初の Generate で凡例が動く（upper_left → left） | `services/session-request.js:4410-4412` の項目順と、`services/config.js:1902, 1845-1884`、`state.js:237` |
| SE-08 | Gallery の tobacco-chloroplast Session では Qualifier Priority の編集が無視される（preserved の `labels.filtering.raw` が優先される）。読み込みにも約 10 秒かかる（preview の前に Pyodide を起動するため） | `services/config.js:1759-1769, 4272-4275`、`services/session-request.js:1272-1280` |

### History（Undo/Redo）

| ID | 症状 | 原因 |
|---|---|---|
| SE-01 | checkpoint 型の履歴（凡例項目が増える色変更、Reset Settings）を Undo/Redo した後に Linear → Circular へ切り替えると例外が出て、Features が次の Generate まで 0 件になる | `services/history-snapshot.js:1322-1324, 1523-1528`、`services/config.js:1045, 1169-1171`。複製された catalog が admit されず、`services/feature-catalog.js:964-972` で例外になる |
| SE-02 | checkbox や radio をラベル文字のクリックで変更すると、Undo できない | `app/history-inputs.js:90-91, 131-137` |
| SE-03 | テキスト欄にフォーカスがある状態のクリック（checkbox、モード切替）が、History から落ちるか、その欄の編集と 1 つにまとめられる | `services/history.js:368`、`app/history-inputs.js:51-88, 93-116, 161` |
| FE-04 | Exact product や protein ID の範囲で非表示にした feature が、Redo や無関係な Undo のあと preview に再び表示される | `app/feature-editor/visibility-actions.js:434, 491`（手動ルールを無視）、`app/app-setup.js:2659` |
| GE-06 | Generate 中に Ctrl+Z を押すと、確定済みの request と Result が古いものに戻り、実行中の Generate は黙って捨てられる | `state.js:793-802`（processing を見ない）、`app/history-shortcuts.js:13-30` |

### Generate とオプション

| ID | 症状 | 原因 |
|---|---|---|
| GE-02 | stroke の即時反映が、空の値・古い値・不正な値を Result に書き込み、Auto に戻らない。その後 Generate が失敗すると「Result は保持された」と表示されるが、Result はすでに変わっている | `app/svg-styles.js:70-78, 407-525, 650-667`、`index.html:4082, 4352, 4493` |
| GE-03 | Run Info の Source recipe では図を再現できない。Circular で定義フォントを 30 にすると、CLI では間隔が変わって CANNOT_FIT になる。Linear では ruler label のフォントが異なる。`--session` による再現は一致する | `gbdraw/circular.py:951-954` と `services/session-request.js:1164-1165` の食い違い、`app/run-info.js:879-880, 965-973`、`gbdraw/linear.py:1165-1166` |
| X-01 | 明確な検証メッセージが UNKNOWN の「Retry / save a Session」に置き換わる（IN-07、GE-05、TR-05、SE-09）。例: "choose a Record before setting a region…"、"Feature Width must be Auto or a positive finite number."、"Invalid session file."。depth の値が数値でないと "Replace or reselect" が 2 回繰り返される。identity 150 などは項目名のない VALIDATION_UNCLASSIFIED になる | `services/error-normalization.js:131, 264-310, 387-425`（NATIVE_VALIDATIONS が足りない）、`gbdraw/web_support/error_adapter.py:43-62, 118, 122` |
| X-02 | 不正な数値が黙って既定値や Auto に置き換わる（GE-04、GE-08、TR-04、CO-09）。E-value `1e-50x` は既定の 0.01 で描かれ、Run Info には `--evalue NaN` と出る。identity −5 → 0、bitscore −1 → 50、長さ 12.5 → 0。GC と Depth の window/step の 0・負の値・小数は Auto になるが、入力欄は入力した値のまま。OIPC-C01/C03 違反 | `app/run-analysis.js:885-896, 2371-2384`、`services/session-request.js:561-564, 2463`、`app/run-info.js:1419` |

### 比較

| ID | 症状 | 原因 |
|---|---|---|
| CO-01 | cross-origin isolation のない環境では、既定の Threaded LOSAT が UNKNOWN で失敗し、serial に切り替わらない。gbdraw.app と `gbdraw gui` は COOP/COEP を送るので、影響は分離されていない環境に限られる | 既定値が `'threaded'`（`services/session-active-config-contract.js:51`、#526 で auto から変更）。`services/losat.js:792-793` が例外を投げる |
| CO-06 | アップロードした表は並び順で record に割り当てられ、ID が照合されない。query と subject を逆にした表では誤った ribbon が描かれ、metadata も矛盾する | `gbdraw/render/groups/linear/pairwise_match.py:436-445` |
| CO-07 | 逆相補にした record の Save Raw LOSAT TSV が元の座標で出力され、CLI や再アップロードでは誤った領域に描かれる | `app/run-analysis.js:1616-1641` |
| CO-08 | Similarity group の名前と説明が不安定な `og_*` ID に紐づくため、グループが組み直されると別のグループへ移るか消える | `app/run-analysis.js:1452-1461`、`app/orthogroups.js:711-725` |

### Feature 編集

| ID | 症状 | 原因 |
|---|---|---|
| FE-02（+PV-09） | batch では、範囲を指定した色・非表示・凡例の編集が表示中の Result にしか反映されず、他の Result の出力は古いまま。次の Generate では凡例の削除だけが全 Result に適用される | `app/svg-styles.js:370-411`、`app/feature-editor/visibility-actions.js:210`、`app/app-setup.js:2266-2313`、`app/candidate-render.js:292-296` |
| FE-03 | batch で Result 2 を表示していても Features drawer は record 1 の一覧を出し、Edit を押しても何も起きない | `state.js:711-725`、`app/feature-editor/svg-actions.js:404-413` |
| FE-05 | Circular の複数 record（single/batch）で Region Annotations の「Selected features」を使うと record が固定され、Generate が失敗する | `app/annotations/target-actions.js:91-106`、`app/annotations.js:71-83`、`app/annotations/validation.js:53-55`、`app/annotations/record-catalog.js:145` |
| FE-06 | Feature 検索の既定の All では、IUPAC の文字だけでできた遺伝子名（CYTB、GATC、dnaA など）が配列に一致し、ほぼすべての feature がヒットする（CYTB は 34/37） | `app/feature-search/search-core.js:139-159, 391-397, 446-466`。Interactive SVG の検索も同じ実装（コードからの判断） |
| FE-07 | `/translation` のない CDS（GFF3 入力はすべて該当）で Copy aa FASTA を使うと、GTG/TTG の開始コドンが M ではなく V/L のまま翻訳される | `gbdraw/web_support/feature_metadata.py:290`。LOSATP の protein 入力も同じ可能性がある（コードからの判断） |

### 凡例・Preview・出力

| ID | 症状 | 原因 |
|---|---|---|
| PV-01 | ドラッグした凡例やタイトルは、Position を変えて Generate すると canvas の完全に外へ出る。出力にも含まれない | `app/legend-layout/decoration-continuity.js:72-107`、`app/legend-layout/composition-actions.js:914-944` |
| PV-02 | feature のない凡例項目（GC content、GC skew ±）の名前を変えても、Generate で元に戻る | `app/feature-editor/color-actions.js:661, 666`、`app/candidate-render.js:228-242` |
| PV-03 | 凡例の並べ替え（Sort、Move up/down）が Generate で失われる。並び順は state、request、session のどこにも保持されない | `app/legend/sort-actions.js:71-128` |
| PV-04 | feature のある凡例項目を既存の名前に変えると、衝突ダイアログを出さずに UNKNOWN エラーになる | `app/feature-editor/color-actions.js:780-797`、`app/legend/entry-actions.js:735` |
| PV-05 | Web の PDF は px を pt として扱うため、CLI（CairoSVG）の PDF より 33% 大きくなり、PNG の DPI とも合わない | `services/export.js:382-396` |

### トラック

| ID | 症状 | 原因 |
|---|---|---|
| TR-02 | 削除・無効化・移動した Circular の Depth 行が、Hide GC などの無関係な切り替えで元に戻る。Linear でも "Add Depth TSV series" で再び有効になる | `app/app-setup.js:2185-2206`、`app/circular-track-slots.js:1712-1798`、`app/depth-track-state.js:255-268, 367-383`、`app/linear-track-slots.js:1318-1398` |
| TR-03 | Depth ファイルをアップローダーの Remove で消すと、管理対象の depth 行が残り、Generate が失敗する | `app/app-setup.js:1803-1819, 1964-1970, 2199-2201` |

## P3（22 件）

| ID | 症状 | 原因 |
|---|---|---|
| FE-08 | 1 つの feature のラベル編集が、Generate 後に大文字小文字だけが違う feature（orfA / ORFA）にも適用される | `app/feature-selector.js:219-221` は大文字小文字を区別するが、Python は区別しない（`gbdraw/labels/filtering.py:292`） |
| FE-09 | Linear で record ID が重複していると、「This feature only」の色編集が効かない | `app/feature-editor/rule-actions.js:462-470`、`app/rule-matching.js:35` |
| FE-10 | キャンセルした Reset fill ダイアログの既定色が、次に別の type を Reset したときに使われる | `app/feature-editor/color-actions.js:1299-1317` |
| FE-11 | 原点をまたぐ feature や分割された feature の位置と長さが誤って表示される（`1..4000`、4,000 bp）。drawer と Location 検索は 0 始まりの座標を使う | `app/feature-editor/svg-actions.js:130-136, 175-181`、`index.html:6606`、`app/feature-search/search-core.js:356` |
| FE-12 | 色に `none` を指定した Specific table の行を、CLI は受け付けるが Web は拒否する | `app/file-imports.js:53-57` |
| GE-07 | Linear で凡例を None にして Generate した後、Position を選ぶと未捕捉の例外が出て、案内が表示されない | `app/legend-layout/reposition-actions.js:116-119`、`app/watchers.js:234-252` |
| GE-09 | JS 側の準備中に Cancel しただけでも、起動済みの Worker が終了し、次の実行が cold start になる | `services/diagram-generation.js:603-616` |
| PV-06 | PDF のテキスト層で、曲線ラベルと目盛のスペースが失われ、検索やコピーができない | `services/export.js:163-209` |
| PV-07 | Generate するたびに canvas の padding が 0 に戻る | `app/app-setup.js:2292-2302` |
| PV-10 | Linear で凡例をその場で移動した配置と、Generate 後の配置が異なる。Generate 前に出力すると移動時の配置で出力される | `app/legend-layout/reposition-actions.js:108-142` |
| PV-11 | Escape で Editor を閉じると、フォーカスが body に落ちる | `app/ui.js:350-358` |
| PV-12 | Legend editor は色をその場で編集できると表示するが、色を変える操作がない（コードからの判断） | `index.html:6512, 6547` |
| CO-10 | Match popup と FASTA ヘッダは crop 後の表示座標を示し、feature popup は元の座標を示す | `app/pairwise-match-popup.js:284-326, 1410-1413`、`app/match-sequences.js:754-760` |
| SE-04 | select にフォーカスがあると、Ctrl+Z などのショートカットが反応しない | `app/history-shortcuts.js:7, 18` |
| SE-10 | Reset Settings で、Linear の record ごとの Definition/Subtitle と alignment plan が残る（Circular では初期化される）。仕様の判断が必要 | `services/reset.js:146-166` |
| TR-06 | stack が無効な間に Hide GC を解除すると、stack を再び有効にしても gc_content 行が無効のまま | `app/circular-track-slots.js:2186-2200` |
| TR-07 | slot の legend label に `,` や ` #` があると、Run Info の CLI コマンドが壊れる | `app/circular-track-slots.js:1149`、`app/linear-track-slots.js:438-439`、`gbdraw/tracks/parsing.py:118-146` |
| TR-08 | 無効な行に、別の行の "(auto)" の寸法が表示される | `app/track-slot-display.js:53-79`、`app/circular-track-slots.js:2437-2445` |
| TR-09 | Show Coordinate Scale を OFF にしていても、最初の stack 有効化とプリセットの Reset で ticks 行が追加される | `app/circular-track-slots.js:1851-1858` |
| TR-10 | アクセシビリティ: 名前のない入力が 8 種類ある（Window、Step、GC Content Mode など）。help tip は 177 個中 175 個が hover 専用で、キーボードやタップでは開けない | `index.html:3086, 3295, 3414, 3425, 3469, 4526, 4532, 4540, 4589, 7429-7430` |
| TR-11 | 1280 px と 1920 px の幅で、Circular stack の行の slot id 欄と renderer 欄が 31 px しかなく読めない | `index.html:3338-3367` |
| TR-12 | Gallery のチュートリアルのリンクが切れている（`7_Linear_Layout.md`、gbdraw.app でも 404） | `gbdraw/web/gallery/tutorials/vibrio-harveyi-group-collinear.json:578` |

## 共通する原因

修正計画のために、原因が共通する項目をまとめる。

1. **複数 Result（batch/grid）で、編集が表示中の Result にしか適用されない、または消える**: FE-01、FE-02、FE-03
2. **即時反映が Result を Generate とは違う内容で書き換える**: IN-01、GE-02、PV-10
3. **Generate で即時編集が引き継がれない**: PV-01、PV-02、PV-03、PV-07、CO-08
4. **JS の高速な record 検出と Python の読み込みが一致しない**: IN-02、IN-03
5. **比較の座標系と再利用の条件**: CO-02、CO-03、CO-07、CO-10（関連: CO-04、CO-06）
6. **エラー変換が足りず、案内が UNKNOWN になる**: X-01、SE-05、SE-06、PV-04、CO-01、IN-04
7. **History のトランザクション境界**: SE-01、SE-02、SE-03、SE-04、GE-06、FE-04
8. **Web と CLI の互換（recipe、session、表）**: GE-03、SE-06、SE-07、TR-07、CO-07、FE-12
9. **Python 側の不具合（CLI も同じく誤る）**: CO-05、TR-01、PV-08、FE-07

## 推奨する着手順

1. CLI も誤る P1（CO-05、TR-01）。修正箇所が小さく、CLI と Web の両方が直る。
2. 編集の消失と誤った出力（FE-01、IN-01）。
3. Generate できない問題（IN-02、IN-03、PV-08）。IN-02 は正規表現を 1 行の中に限定すれば直る。
4. 比較の正しさ（CO-02、CO-03、CO-04）。
5. X-01（エラー変換）。多くの P2/P3 の使い勝手がまとめて改善する。

## 問題がなかった範囲

- **XSS**: GenBank の各 qualifier、ファイル名、細工した Session JSON に入れたペイロードは、preview、popup、Interactive SVG（file:// で開いたもの）のいずれでも実行されなかった。
- **Session**: Gallery の全 10 Session で Save → Load → Save が一致した。Undo/Redo と Generate の混在、モードの往復、読み込み失敗時の状態保持も正しい。
- **CLI での再現**: Circular 32 件、Linear 25 件の Source recipe を dev の CLI で再生し、GE-03 以外は一致した。
- **比較**: outfmt 6/7、CRLF、空ファイル、LOSAT cache の再利用、threaded と serial の結果の一致、Cancel はいずれも正しい。
- **配列のコピー**: ± 鎖、原点をまたぐ join、codon_start=2 はいずれも正しい（FE-07 を除く）。
- **その他**: PNG の寸法（72/300/600 DPI）は正しい。390〜1920 px で横スクロールは出ない。読み込み時に console error は出ない。Gallery のリンクは 159 件中 158 件が正常。

## 未検証

- macOS の Cmd キー、実機でのタップ、非常に大きいゲノム
- Circular conservation、GFF3+FASTA での LOSAT、TLOSATX の遺伝暗号キャッシュ
- batch での Label TSV 取り込み（コードからの判断のみ）、record 自体のドラッグ
- IN-02 の main（リリース版）での挙動

## GUI の外で見つかったもの

- CLI の `--block_stroke_width -1` は traceback を出して終了する。
- 負の label font size や stroke 幅を GUI と CLI の両方が受け付け、そのまま SVG に書き込む。
- CLI の `--window 0` や `--step 0` では GC track が空になる。
- Dinucleotide に `G` を指定すると CLI が IndexError で終了する。

## 証拠と再現方法

- **領域ごとの詳細**（再現手順、観測値、証拠ファイル名）: `/home/kawato/gbdraw-baselines/web-gui-audit-20260930/evidence/<area>/FINDINGS.md`。area は input、session、generate、comparison、feature、preview、tracks のいずれか。
- **再現 spec**: `/home/kawato/gbdraw-baselines/web-gui-audit-20260930/specs/audit-<area>/`。リポジトリには置かない。修正の PR で、assert を持つテストとして `tests/web/` に取り込む。
  - dev worktree の `tests/web/audit-<area>/` に置き、`GBDRAW_WEB_TEST_PORT=<port> npx playwright test tests/web/audit-<area>/<spec> --workers=1` で実行する。
  - spec は assert せず、観測値を JSON に書き出す。
  - fixture と出力先は scratchpad の絶対パス（`/tmp/claude-1000/...`）を参照しているため、再実行するときは書き換える。
- **比較用の大きな SVG（約 2 GB）**: 一時領域の scratchpad にしか残していないため、消えることがある。

## 設計段階での訂正（2026-09-30）

修正案を作る段階でコードを読み直し、次の点を訂正した。詳細と、新たに見つかった不具合（N-01〜N-16）は [01_REMEDIATION_PROPOSAL.md](01_REMEDIATION_PROPOSAL.md) の第 3 節にある。上の表の記載は監査時点のまま残す。

- **PV-08:** 原因は定義文字列用の予約ではない。Multi-Record Canvas が depth 入力なしでも空の depth スロットを確保していることが原因（`gbdraw/api/diagram.py:2264, 2897, 2994`、`021b4c49` からの退行）。この状態は dev の CLI で確認した。その結果、Web 既定の Circular 出力はすべて単一 record 経路と幾何が違う（N-01）。
- **CO-05:** conservation の reader もファイル入力では誤読する（上の表で直した）。
- **IN-02 と IN-03:** 同じ正規表現が `app/run-analysis.js:772-774` と `app/match-sequences.js:820-822` にもある。IN-03 は Worker の helper（`app/python-helpers.js:1645-1660`）も原因の一つ。
- **IN-06:** 行番号は `app/linear-sources.js:107-145`。
- **既存の決定との関係:** PV-01 は PD-OI-052 が受容した残余リスクで、欠陥ではない。PV-03 と PV-07 は、同じ決定が継承を保証しないと明記している（機能要望）。
- **TR-09 と TR-10(b):** TR-09 の「初回の stack 有効化」は UI 契約どおり。TR-10(b) は意図された設計。
- **CO-01 と CO-06:** serial へ fallback しないことと、位置で割り当てることは契約どおり。欠陥は、それぞれ UNKNOWN の表示と metadata の矛盾。
- **X-02:** window/step の 0 と負値は CLI も受理する。GE-08 は TR-04 と同じ P2 に揃える。
