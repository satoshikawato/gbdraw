<!-- Raw design report of workstream W8 (sub-agent output, 2026-09-30, base origin/dev 4c89bab1). Evidence for 01_REMEDIATION_PROPOSAL.md, not Product or architecture authority. -->

# W8 報告：65 件が dev 4c89bab1 の CI を通過した理由と、クラスごとの再発防止策

注記：並列で起動した 3 つの調査エージェントの結果を待たずに締め切った。
- 各クラスの「最も近い既存テスト」は、主要なファイルを自分で読んだ範囲に基づく。
- C1、C3、C5、C6、C7 の一部は網羅的に確かめていない。該当箇所には「要確認」と書いた。
- 推論には次の根拠を使った。不具合は決定的に再現し、CI は緑だった。したがって「その入力と経路を実行し、正しい結果と比べるテスト」は存在しない。
- 実行時間は計測していない見積りである。

## 0. 要旨
- **主因は 3 つ。**
  1. テストの多くは、Result と mounted SVG と export が互いに一致すること（**内部の一貫性**）を検査している。Generate し直した結果や CLI の結果と比べる**正解の基準**を持たない。
  2. テストが単一 Result、単一 record、既定モードの経路に偏っている。
  3. 規範（OIPC-C01/C03/C04/C06/C07、OIC-021 AC-15、`gbdraw/web/CLAUDE.md:69-71` の single-thread fallback）はあるが、実行できる汎用テストがない。
- **提案は 10 個のガード群。** 新しいフレームワークや workflow、規約の層は作らず、既存の仕組みを拡張する。
  - 使う既存資産：`generateAndWaitForResult`、`svg-visual-semantics.mjs`、`tests/utils/svg_compare.py`、`tests/fixtures/*.json` の JS/Python 共有ベクター、privileged capability の allowlist、pytest の strict xfail、Playwright の `test.fail`。
  - 65 件のうち**約 50 件はクラス単位のガードで検出できる**。残り約 15 件は個別の回帰テストで扱う。
- **PR smoke は上限に達している**（19/19、`tests/ci/playwright-inventory.test.mjs:23-26`）。
  - 速い検査は Node/pytest で PR 段階に置く。
  - ブラウザの行列テストは dev SHA staging（functional-full）に置く。
  - 57 probe の完全な sweep は release 段階か、main 昇格前の監査で回す。
- **architecture rule の登録枠も上限に達している**（`MAXIMUM_RULE_COUNT = 4`、`tools/web-architecture-evaluation.mjs:2`。4 件登録済み）。静的ガードは既存の privileged capability の allowlist を狭める方向で作る。
- **修正の進め方**は、決定パック 4 つと検査基盤の PR を先に出し、その後 P1→P2→P3 の順に進める。合計約 50〜55 PR、dev staging と同期したマージ列で流す。

## 1. ギャップ分析

### 1.1 系統的な原因
| # | 原因 | 根拠 |
|---|---|---|
| S1 | **内部の一貫性だけで、正解の基準がない** | `tests/web/helpers/visual-state.cjs:11-40` の `capture`/`assertCoherent` は、選択中の Result、mounted、export の 3 つが一致するかだけを見る。`tests/web/mode-transition-result.playwright.spec.js:16-56` は、編集後の Result をモード往復後の自分自身と比べる。IN-01、GE-02、PV-10 のように 3 つとも同じ誤った内容になると通過する。 |
| S2 | **比較の範囲が狭い** | `tests/web/helpers/mode-transition.cjs:101-110` の `semantics` は `[data-gbdraw-feature-id]` の d/fill/stroke/display しか比べない。定義行、凡例、目盛、タイトルは対象外なので、IN-01 は原理的に検出できない。 |
| S3 | **Web の既定 topology（grid）が Python テストの既定にない** | Web は record が 1 つでも Multi-Record Canvas（grid）を使う。TR-01 と PV-08 は `--multi_record_canvas` のときだけ起きる。CLI も同じく誤るのに CI が緑なので、grid で slot 凡例や長い `/organism` を検査するテストは存在しない。Product Impact map 自身も残存リスクとして「one representative ... not exhaustively」と書いている（`tools/web-product-impact-map.json:129,139,148,439,452`）。 |
| S4 | **個別の回帰テストが一般化されていない** | 前回監査（2026-09-17）の回帰は `tests/web/gui-audit-regressions.playwright.spec.js` に 1 件 1 テストで積まれた。たとえば `:118` はグループ名がモード変更を生き残ることを検査するが、グループの組み直し（CO-08）は対象外。同じクラスの不具合が別のきっかけで再発した。 |
| S5 | **高速化で入れた JS 再実装に、意味の等価性の証拠がない** | `gbdraw/web/js/app/record-discovery.js` の高速経路は 0ce21a32「Fix and speed up multi-record rendering」で入った。共有ベクター `tests/fixtures/record_metadata_inference_cases.json:2-7` は「the two implementations cannot drift apart」と明記するが、対象は定義文字列だけで、record ID と record の列挙を含まない（→ IN-02、IN-03）。 |
| S6 | **エラー分類が文字列照合で、生成側と変換側の対応を検査しない** | `gbdraw/web_support/error_adapter.py:42-62`（`_CONSTRAINTS`）と `:205-210`（該当しなければ `VALIDATION_UNCLASSIFIED`）。`gbdraw/mode_profiles.py:43` の "identity must be a finite number in [0, 100]." はどの条項にも一致しない。監査員は手書きの一覧（`evidence/input/evidence/error-normalization-check.mjs`）で漏れを見つけた。 |
| S7 | **規範が実行できない** | OIPC-C01（不正値は拒否し、黙って補正しない、`docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md:285-291`）、C03（`:299-302`）、C04（cache identity、`:304-310`）、C06（空値や非 active モードは削除を意味しない、`:320-324`）、C07（`:326-329`）、OIC-021 AC-15（比較図形も同じ変換を使う）。どれも汎用の性質テストがなく、各 OIC 項目で選んだシナリオだけが検査されている。 |
| S8 | **検査の投資先が統治の仕組みと計算回数に偏っている** | Web の Node テスト約 670 ケースのうち約 239 ケースは統治機構そのものの検査（`architecture-contracts.test.mjs` 143 ケース / 5,278 行、ratchet、promotion 系）。CW-01..06 と `session-regeneration-contract.cjs` は計算回数とライフサイクルを測るが、出力の正しさは測らない。 |
| S9 | **速い層は枠が埋まり、重い層はマージ後にしか走らない** | PR smoke は 19/19。規約文書の記述は古いまま（`docs/internal/SELECTIVE_CI.md:76` は "thirteen"、`:124` は "8 or more than 12"）。functional（395 ケース、CI では retries 2、`playwright.config.js:11`）は dev push のときだけ走る。 |
| S10 | **不変条件が「現在の Result」しか定めていない** | `gbdraw/web/CLAUDE.md:145-156` は「update it and the current Result synchronously」とだけ書く。batch や grid の他の Result の扱いが規範にないため、クラス 1 は製品判断が必要になる。 |
| S11 | **テストの入口が 1 か所にない** | Web の「Changing a setting」チェックリスト（`gbdraw/web/CLAUDE.md:235-252`、6 番目が `:244`）に次の 4 点がない。モードの行列（single/grid/batch × Circular/Linear）、不正値の拒否、CLI での再現、即時編集と Generate の等価性。 |

副次的な発見：
- `tests/web/helpers/session-regeneration-contract.cjs:10` の `number()` は、存在しない metric を 0 として扱う。そのため `exact(key, 0)`（`:40-43`）は metric がなくても通る。これは CW の原則「A missing metric fails」（`docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md:390-446`）に反する可能性がある（要確認）。
- `docs/internal/ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md:592` は「at most three rules」と書くが、実装の上限は 4。

### 1.2 クラスごとの分析
| クラス | 最も近い既存テストと検査内容 | 漏れた点 |
|---|---|---|
| 1 複数 Result（FE-01/02/03） | mapped contract の `tests/web/contracts/active-result-edit-transaction.playwright.spec.js` は、tRNA の fill を History、Session、Generate、export まで追う。単一 Result である（map の残存リスク `tools/web-product-impact-map.json:439`）。`tests/web/right-drawer.playwright.spec.js:702`（smoke）も単一 Result。functional のタイトルに batch/grid/multi-record を含むものは 395 件中 45 件あるが、編集 → Result 切替 → Generate → Save を組み合わせるものは確認できなかった（全件の精査は未了）。`syncLabelEditor`/`clearOverrides` を直接呼ぶ Node テストは、見た範囲（`label-override-empty`、`feature-override-identity`）にはない。 | Result topology の次元がない。表示を同期する経路が canonical な状態を書き換えないこと（観測だけで変化しないこと）を検査していない。 |
| 2 即時編集が Result を書き換える（IN-01、GE-02、PV-10） | `visual-state.cjs` の一貫性検査と `mode-transition-result` の自己比較（S1、S2）。即時編集と Generate の一致を検査するのは、色 rule についての `tests/web/python-rule-parity.playwright.spec.js:18,138` だけ。 | 「即時編集した Result ≡ 同じ draft で Generate した結果」を比べていない。crop、逆相補、record ラベル、grid のような前提条件を持つ fixture もない。構造面では、Result を直接書き換える箇所が 14 モジュールに約 20 か所ある（`grep 'results.value = nextResults'`）。canonical な `gbdraw/web/js/app/preview-runtime.js:738-767`（`flushActiveResult`）を通らない複製がある（例: `app/svg-styles.js:70-78`、`app/results.js:354-361`）。allowlist の "Mounted SVG/Result replacement" は 18 owner を許可している（`tools/web-change-policy.json`）。 |
| 3 Generate が即時編集を引き継がない（PV-01/02/03/07、CO-08） | 即時編集のテストは、編集直後の mounted と Result を見る。Generate 後に編集が残ることを見るのは feature の色（上記 contract）と Python rule だけ。`gui-audit-regressions:118` はグループ名をモード変更についてだけ見る。 | 編集 → Generate → 保持の往復が、凡例、decoration、padding、group 名に一般化されていない。 |
| 4 JS と Python の record 検出の不一致（IN-02、IN-03） | 共有ベクター `record_metadata_inference_cases.json` を `tests/test_record_metadata.py:22-29` と `tests/web/record-metadata-inference.test.mjs:15-18` が読むが、対象は定義文字列だけ。 | record ID と record 集合のベクターがない。Prokka 形式（ACCESSION/VERSION が空）や、FASTA に注釈のない配列を含む GFF の fixture が `tests/test_inputs` にない（ファイル名で確認）。 |
| 5 比較の座標系と再利用の条件（CO-02/03/07/10、関連 CO-04/06） | `tests/ci/playwright-inventory.test.mjs:45-46` の "protein raw cache survives cancellation and derived options preserve search identity" は**再利用されること**を検査する。OIC-019 が無効化の次元として挙げるのは raw 設定、入力、Clear Cache、Session/History だけ（`OPTION_INTEGRITY_PRODUCT_CONTRACT.md` の catalog）。OIC-021 AC-14/15 は規範として要求されている。 | 次元ごとの感度（visibility を変えたら再利用しない）を検査する負のテストがない。Orient forward ON で LOSATN の ribbon 位置を見るテストは見つけられなかった（要確認）。ファイルのまとめ方による差（CO-04）、raw TSV の往復（CO-07）、ID による照合（CO-06）の性質テストもない。 |
| 6 エラー変換の漏れ（X-01 ほか） | 対応表が文字列照合（S6）。`generateAndWaitForResult`（`tests/web/helpers/app-lifecycle.cjs:332-362`、15 の spec が使う）は、エラーのとき summary が空でないことしか見ない。 | すべての検証メッセージが mapped かどうかの網羅テストがない。UNKNOWN、pageerror（GE-07）、Run Info の `NaN`（X-02）も、共通 helper の自動検査で検出されていない。 |
| 7 History の境界（SE-01〜04、GE-06、FE-04） | `history-inputs.*`、`history-generated-authority` などは存在する。監査で再現に使った入力経路（ラベルのクリック、text にフォーカスしたままのクリック、select にフォーカスした状態のショートカット、Generate 中の Ctrl+Z）の次元は、見た範囲では確認できなかった（要確認）。 | 入力方法 × 処理中という次元がない。Undo^n Redo^n の往復で状態が元に戻るかも検査していない。 |
| 8 Web と CLI の互換（GE-03、SE-06/07、TR-07、CO-07、FE-12） | `tests/test_run_info_exact_replay.py` は `--session` の exact replay（Issue #469）だけを実行する。監査でも一致した。`tests/web/run-info.test.mjs:807-943` は recipe の文字列と `unavailableReason` を検査するが、生成した Source recipe を CLI で実行しない。 | 「Source recipe を CLI で実行した結果 ≡ GUI の結果」「CLI の session を読み込んで Generate した結果」の検査がない。表の読み込み（color `none`）も JS/Python の共有ベクターがない。 |
| 9 Python 側の不具合（CO-05、TR-01、PV-08、FE-07） | outfmt 6 の reader が 4 つある。どれも `names=COMPARISON_COLUMNS`：`gbdraw/io/comparisons.py:73-78`、`gbdraw/session_request_codec.py:3330-3335`、`gbdraw/api/record_planning.py:1426-1431`、`gbdraw/analysis/conservation.py:175-180`。pandas で試したところ、14 列では先頭 2 列が MultiIndex になり、列が 2 つずれた。conservation の `:147-149` は DataFrame を受け取る経路の対処で、ファイルを読む `:175` でも同じずれが起きる可能性がある（要確認）。TR-01 は単一 record 経路（`gbdraw/diagrams/circular/assemble.py:2766`）と grid 経路（`gbdraw/api/diagram.py:3114-3128`）が別になっている。 | 13〜14 列、`--multi_record_canvas` での slot 凡例、Web 既定 profile と長い `/organism` の組合せ、`/translation` がなく GTG で始まる CDS。いずれも実行するテストが存在しない（CI が緑であることから推論）。 |

## 2. 再発防止策（クラスごとの最小限のガード）

実行段階の略称：
- PR = PR から dev の速い検査（Node/pytest）
- STG = dev SHA staging（functional-full / browser / core）
- REL = release 段階または main 昇格前の監査

時間は未計測の見積り。

| ID | ガード | owner（所有ファイル） | 段階 | コスト | 検出できた不具合 |
|---|---|---|---|---|---|
| G-A | **即時編集 ≡ Generate（メタモルフィック）**。即時編集の表（Species/Strain/定義フォント/Plot title、stroke の設定と解除、凡例のドラッグと Position、凡例名の変更（feature のない項目を含む）、並べ替え、padding、OG 名）× Circular/Linear。即時編集後の Result を `svgVisualSemantics` で取り、Generate 後と `compareVisualSemantics` で比べる。さらに Save → Load → Generate 後も編集が残ることを見る。fixture は crop、逆相補、record ラベル、grid を含むものにする | 新規 `tests/web/contracts/live-edit-regeneration-equivalence.playwright.spec.js`。比較は既存の `tests/web/helpers/svg-visual-semantics.mjs:24-72` | STG。PR smoke には代表 1 件を入れ替えで置く | 約 2〜3 分 | IN-01、GE-02、PV-10、PV-01、PV-02、PV-03、PV-07、CO-08（metadata 版）、SE-07 |
| G-B | **Result topology の行列**。{Circular single/grid/batch、Linear multi-record} × {ラベル文字、ラベル表示、範囲を指定した色、非表示、凡例}。Result 2 で編集 → 1 と 2 を往復 → Generate → Save/Load。override が残ること、drawer が表示中の record を出すこと、範囲指定が各 Result の export に反映されること（DP-1 に従う）を見る | 新規 `tests/web/multi-result-edit-matrix.playwright.spec.js`。mapped contract とは別ファイルにして、mapped evidence の制約（`docs/internal/WEB_CHANGE_POLICY.md:393-397`）を避ける | STG | 約 2〜3 分 | FE-01、FE-02、FE-03、FE-05、PV-09、IN-01(c) |
| G-C | **無断変更の検出**。編集ではない操作（Result の選択、変更なしの Generate、モード往復、Undo+Redo の対、無関係な切替）の前後で、利用者が持つ状態（`buildConfigData`/`buildUiStateData`/`buildEditorStateData`。catalog は除く。監査員の `harness.snapshot` と同じ）を比べ、差がないか、History に記録されていることを求める | `tests/web/helpers/app-lifecycle.cjs` に helper を追加し、表形式の spec を 1 本 | STG（Node で書ける部分は PR） | 約 1 分 | FE-01、IN-05、TR-02、TR-06、TR-09、PV-07、SE-07、FE-04 |
| G-D | **JS/Python 共有ベクターの拡張**（既存の `tests/fixtures/*.json` の方式を使う）。record 検出は小さな合成ファイル（Prokka 形式、注釈のない配列を含む GFF+FASTA、複数 record、CRLF/BOM、ID の重複）とし、期待値は pytest が `gbdraw.io.genome` の loader から検証する（Python が正、`gbdraw/web/CLAUDE.md:44-45`）。Node が `record-discovery.js` と照合する。同じ方式を color table の `none`、ラベル照合の大文字小文字、outfmt の 12/13/14 列と outfmt 7 の header（4 つの reader すべて）にも使う | `tests/fixtures/`、`tests/test_genome_loading.py`、`tests/web/record-metadata-inference.test.mjs`（または兄弟ファイル）、`tests/test_comparisons.py` | PR | 5 秒未満 | IN-02、IN-03、IN-06、FE-12、FE-08、CO-05 |
| G-E | **Source recipe ↔ CLI の一致**。監査員の `parity-circular/linear` と `replay.py` を 1 本にまとめる。gallery parity と同じく、spec から `python -m gbdraw.cli` を起動し、`tests/utils/svg_compare.compare_svgs` で比べる。代表 10 probe（オプション族ごとに 1 つ。def_font_size、ruler label フォント、`,` や ` #` を含む slot ラベルを含む）。CLI の session については、`tests/web/session-cli-compatibility.playwright.spec.js` に Linear+BLAST の `--session_output` → 比較付き Generate、および凡例位置を追加する | 新規 `tests/web/contracts/source-recipe-cli-parity.playwright.spec.js`。57 probe の完全版は tool として残す | STG（完全版は REL） | 約 3〜4 分（完全版は約 20 分） | GE-03、TR-07、SE-06、SE-07、CO-07（再アップロード版） |
| G-F | **比較のメタモルフィック**。(1) cache key の感度表：visibility table、region、reverse、`reverseComplementOverride`、gencode、record 選択を変えたら再利用しない。表示だけの変更では再利用する（CW-06 の無効化テストを兼ねる）。(2) rotate+orient ≡ region_reverse の投影（`gbdraw/web/js/app/run-analysis.js:3177-3208`）。(3) 同じ record を 1 ファイルにまとめても別ファイルにしても同じ batch 計画とリンクになる（`app/linear-sources.js:151-212`）。(4) raw TSV を保存して再アップロードしても同じ ribbon になる。(5) query/subject を入れ替えた表や行を並べ替えた表でも、同じ結果になるか拒否される（pytest、`gbdraw/render/groups/linear/pairwise_match.py:436-445`）。(6) match popup、FASTA header、feature popup の座標系が一致する | `tests/web/run-analysis-derived-cache.test.mjs`、`tests/web/linear-sources.test.mjs`、`tests/web/match-sequences.test.mjs`、`tests/test_linear_multi_record_comparisons.py`。ブラウザで確かめる 2 ケースは functional | PR（Node/pytest）、STG（ブラウザ） | PR 数秒、STG 約 2 分 | CO-02、CO-03、CO-04、CO-06、CO-07、CO-10 |
| G-G | **エラーと数値の網羅**。(1) pytest：AST で `raise ValidationError(...)` を集め、placeholder を埋めたメッセージを `serialize_web_error` に通し、UNKNOWN/UNCLASSIFIED にならないことを求める。例外は縮小しかできない allowlist にする。(2) Node：検証を担う JS ファイルの `throw new Error` を同様に `normalizeUserFacingError` に通す（監査員の `error-normalization-check.mjs` を種にする）。(3) **共通 helper の自動検査**：`generateAndWaitForResult` と Save の helper で、error code ≠ UNKNOWN（`allowUnknown` で明示的に外せる）、pageerror がない、Run Info に `NaN`/`undefined` がないことを確かめる。15 の spec が自動で検査を受ける。(4) 数値の表（OIPC-C01/C03）：request に投影される数値欄すべてに NaN、`1e-50x`、-5、0、12.5 を与え、拒否されるか、Python が拒否するようリテラルのまま渡されることを求める。既定値への黙った置き換えは禁止 | `tests/test_web_error_adapter.py`、`tests/web/error-normalization.test.mjs`、`tests/web/helpers/app-lifecycle.cjs:332-362`、新規 `tests/web/option-input-integrity.test.mjs`（mapped の `session-request.test.mjs` は避ける） | PR（自動検査は STG の既存ケースに乗る） | 各 2 秒未満 | X-01（IN-07、GE-05、TR-05、SE-09）、SE-05、PV-04、IN-04（メッセージ部分）、IN-08、GE-07、X-02（GE-04、GE-08、TR-04、CO-09） |
| G-H | **History の境界の行列**。入力方法（直接クリック、ラベルのクリック、キーボード、text にフォーカスしたままのクリック、select にフォーカスした状態のショートカット）× control の種類、および処理中の次元（Generate 中の Undo は拒否し、確定済みの request と Result を変えない、OIPC-C07）。Undo^n Redo^n の往復の後にモードを切り替えても壊れないこと（G-C の検出を再利用） | `tests/web/history-inputs.playwright.spec.js`、`history-generated-authority.playwright.spec.js` | STG | 約 1〜2 分 | SE-01、SE-02、SE-03、SE-04、GE-06、FE-04 |
| G-I | **Python コアの性質テスト**。(1) record 1 つの grid ≡ single：凡例項目と slot の集合が一致し、single が成功するなら grid も成功する。(2) Web 既定 profile（`--multi_record_canvas --track_type tuckin --gc --skew --separate_strands`）× 長い `/organism` の合成 fixture。(3) `/translation` を除いた GenBank の aa ≡ 元の `/translation`（GTG/TTG）。(4) 数値 option の制約（負の stroke、window 0、Dinucleotide `G`）が、traceback ではなく名前付きの ValidationError になる | `tests/test_circular_multi_canvas.py`、`tests/test_circular_track_slots.py`、`tests/test_web_feature_metadata.py`、`tests/test_mode_profiles.py` | PR（core-pr） | 合計 10 秒未満 | TR-01、PV-08、FE-07、および「GUI の外で見つかったもの」 |
| G-J | **既存の仕組みだけで作る静的ガード**。(1) privileged capability の "Mounted SVG/Result replacement" を 2 つに分ける。"Result content commit"（`results.value =`、`state.results.value =`、`flushActiveResult(`）と "SVG serialization"（`serializeCleanSvg(`）。修正 PR ごとに commit 側の allowlist を狭め、最終的には `preview-runtime.js`、`run-analysis.js`、`services/config.js`、`services/history-snapshot.js` 程度にする。rule registry の枠（満杯）は使わない。(2) CO-05 の統合後に、outfmt の reader が 1 つだけであることを pytest で検査する（前例 `tests/test_ci_import_boundaries.py`）。(3) `tests/test_documentation_reference_contracts.py:241-248` の検査対象に `gbdraw/web/gallery/tutorials/*.json` の href を加える。(4) a11y：すべての input/select に accessible name があり、help tip に keyboard でフォーカスできる | `tools/web-architecture-detectors.mjs` と `tools/web-change-policy.json`（分離の手順は §4）、pytest 2 本、Playwright 1 本 | PR（静的検査）、STG（a11y） | 数秒 | クラス 1/2 の再発を構造的に防ぐ、CO-05 の再発、TR-12、TR-10 |

不採用（YAGNI）：
- `clearOverrides()` を view-sync 経路から呼ぶことを禁じる静的 detector は採らない。3 か所の呼び出しがすべて `gbdraw/web/js/app/feature-editor/label-actions.js`（`:699`、`:990`、`:1077`）にあり、ファイル単位の detector では区別できない。G-C の振る舞いテストで代える。
- 新しい rule kind を追加して registry の枠を広げることも採らない。スキーマの計画が必要になり、費用に見合わない。

個別の回帰テストで扱うもの（約 15 件）：FE-06、FE-09、FE-10、FE-11、GE-09、PV-05、PV-06、PV-11、PV-12、SE-08、SE-10（仕様の判断が必要）、TR-03、TR-08、TR-11、CO-01。
- CO-01 は COOP/COEP なしで既定の threaded を実行し、serial に切り替わることを確かめる spec 1 本で扱う。
- 前提として、Playwright の webServer は `python3 -m http.server`（`playwright.config.js:21-26`）で COOP/COEP を送らない。既存の LOSAT テストが並列方式をどう指定しているかは要確認。

検出できる件数：P1 は 10/10、P2 は約 28/33、P3 は約 12/22。合計約 50/65。

## 3. 監査手順の定着

**`tests/web/` に取り込むもの。** いずれも assert を持つ形に書き直し、`/tmp` の絶対パスを除き、fixture を小さな合成ファイルとして `tests/fixtures/` に置く。

| 監査員の資産 | 取り込み先 |
|---|---|
| `audit-generate/parity-*.spec.js`、`replay.py` | G-E（比較器は `tests/utils/svg_compare.py` に統一し、`svgcmp.py` は捨てる） |
| `audit-generate/live-*.spec.js`、`audit-input/in01*.spec.js`、`audit-preview/legend-persistence.spec.js`、`layout-persistence.spec.js` | G-A |
| `audit-feature/batch*.spec.js`、`annotation-selected.spec.js` | G-B |
| `evidence/session/harness.py` の `snapshot`/`diff` | G-C（Node helper に移植） |
| `audit-tracks/depth-slot-*`、`slot-history*` | G-C |
| `evidence/input/evidence/error-normalization-check.mjs`、`fastpath.mjs` | G-G、G-D |
| `evidence/session/e01*`、`e02*`、`e16*`、`history-adapter.spec.js` | G-H |
| `audit-comparison/losatp-reuse`、`rotate-orient`、`source-batch`、`raw-export`、`swapped`、`popup` | G-F |
| `audit-tracks/a11y-probe.spec.js` | G-J(4) |
| 合成 fixture（prokka_like_single.gbk、gffpair.*、multi.gb、case.gb、dup.gb、none_color.tsv、`make_fixture.py`） | `tests/fixtures/`（数 KB） |

**tool として残すもの（CI には入れない）：**
- 57 probe の完全な sweep
- XSS のペイロード sweep（`audit-preview/xss.spec.js`。sanitizer 自体は `svg-sanitization.test.mjs` が担う）
- 画面幅の sweep（responsive/gallery-responsive）
- Gallery 全 session の往復
- 読み込み時間の測定
- `cdp-*`/`debug-*` の調査用スクリプト

これらは `tools/audit/`（新設の 1 ディレクトリ）に置き、実行方法を README 1 つにまとめる。

**予算：**
- PR smoke には追加しない（19/19）。
- 代表ケース（G-A の 1 件）を入れる場合は、既存の 1 件と入れ替える。上限を 21 に上げる場合は、`tests/ci/playwright-inventory.test.mjs:23-26` と `docs/internal/SELECTIVE_CI.md:76,124` の修正を 1 つの PR で行う。
- PR 段階の追加は Node/pytest だけとし、合計 15 秒未満を目安にする。
- STG の追加は約 12〜15 分で、4 shard に分散させる。shard 1 はすでに 30 分を超えた記録がある（`docs/internal/SELECTIVE_CI.md` の「Smoke inventory and regression retention」節）。新しい spec は `test.describe.configure` とファイル分割で偏らせず、JSON report で所要時間を確認する。

**定期監査（最小限の定義）。** 新しい gate や workflow は作らず、昇格 PR 本文のチェック項目にする。
- **契機**: 各 dev→main 昇格の前。または前回監査から runtime PR が約 30 件溜まったとき。
- **入力**: 前回昇格からの ci-impact 計画に現れた capabilities の和集合。
- **手順**: (1) release 段階の dispatch（既存）、(2) `tools/audit/` の sweep、(3) 変更された領域の探索的監査（時間を区切る）。
- **出力**: 今回の README と同じ形式の `docs/internal/web-gui-audit-<date>/`。確認した不具合のクラスごとに「ガードを追加する」か「ガードを置かない判断」を記録する。
- **終了条件**: P1 はすべて修正済みか、Owner が明示的に waiver を出したもの。

## 4. 修正の進め方（ワークフロー）

**前提となる制約：**
- Ordinary profile は 8 ファイル / churn 800 / 純増 100、Architecture profile は 12 / 1,500 / 400（`docs/internal/WEB_CHANGE_POLICY.md:151-165`）。数えるのは production scope だけで、テストは含まない。
- checker の実装と authority を同じ PR で変えてはならない（同 `:196-243`）。
- mapped contract は evidence を先に出す（同 `:393-397`）。
- 許可範囲の拡大は authority PR を先に出す。縮小は runtime と同時でよい（同 `:407-431`）。
- 未解決の Product 判断があるものは実装を止める（`docs/internal/PRODUCT_IMPACT_RATCHET.md:296-317`）。
- 狭い PR 経路を使うには、base の exact SHA に staging の証拠が必要（`docs/internal/SELECTIVE_CI.md`「Evidence required for selection」）。
- dev への push は concurrency で古い staging を cancel する（`.github/workflows/test.yml` の concurrency）。

**Phase 0（証拠・判断・基盤、約 8 PR ＋ 決定パック 4 つ）**

0-1. policy-documentation PR 1 本。smoke 件数と rule 上限の記述を実装に合わせる。

0-2. Decision Pack を 4 つにまとめ、Owner に一括で提示する。Lane B、`docs/internal/PRODUCT_DECISION_PACKET_TEMPLATE.md` を使う。

| パック | 対象 | 判断すること |
|---|---|---|
| DP-1 複数 Result の編集範囲 | FE-02/PV-09、FE-03 | 範囲を指定した編集を全 Result に即時反映するか |
| DP-2 即時編集 ≡ Generate を契約にする | IN-01、GE-02、PV-10、PV-01〜03、PV-07、CO-08 | 等価性を OIC に登録するか。PV-03 の並び順を session に持つか（スキーマの変更） |
| DP-3 エラーと既定値 | SE-10 | Reset の範囲。X-02 は OIPC-C01/C03 で、CO-01 は `gbdraw/web/CLAUDE.md:69-71` で既に決まっているので `IMPLEMENT_EXISTING_AUTHORITY` |
| DP-4 Web と CLI の互換 | SE-06、PV-05、CO-04/CO-06、FE-12 | SE-06 は CLI session の ID の互換 reader が要るか。永続形式の方針に従い、main の履歴または release tag での存在証拠が要る。PV-05 は PDF の単位 |

- 科学的な正しさの問題（CO-05、TR-01、FE-07、CO-02、CO-03）は、決定的な規則があるので `IMPLEMENT_EXISTING_AUTHORITY` とし、判断を待たない。
- 判断を受けたら、パックごとに authority だけの OIC PR を 1 本出す（合計 4 本、Review REQUIRED）。

0-3. 検査基盤の PR（tests-only、約 4 本）：
- G-C と G-G(3) の helper
- G-D のベクター
- G-A、G-B、G-H の spec
- G-G(1)、G-G(2)、G-I の pytest と Node テスト

既知の不具合は次の形で「期待どおり失敗する」と印を付ける。
- Playwright：`test.fail(true, 'IN-01')`
- pytest：`xfail(strict=True, reason='CO-05')`
- ベクター：`knownDefect: 'IN-02'` 欄。この欄がある間は乖離が続いていることを assert する

**修正 PR はその印を外すことが完了の条件になる。** 直ったのに印が残っていれば、テストが失敗して知らせる。

0-4. G-J(1) の順序：
1. authority だけの PR：`tools/web-change-policy.json` に新しい key を現在の owner 全体で追加する。
2. checker だけの PR：detector の分離と `architecture-contracts.test.mjs:516` の特性記録を更新する。

注意：key を先に追加した状態を checker が許容するかは要確認。許容しない場合は、同じ policy の規則内で順序を入れ替える。

**Phase 1（P1 10 件、約 9 PR）**

| 順 | 内容 | 補足 |
|---|---|---|
| 1 | CO-05 | reader を 4 つから 1 つの owner に統合し、G-J(2) を付ける |
| 2 | TR-01 | |
| 3 | FE-07 | P2 だが Python の同じ系統なので先に行う |
| 4 | IN-02+IN-03 | 同じ `record-discovery.js` なので 1 PR |
| 5 | FE-01 | |
| 6 | IN-01 | committed request から定義を作り直す。architecture-change label の見込み。DP-2 の後 |
| 7 | CO-02 | |
| 8 | CO-03 | |
| 9 | CO-04 | |
| 10 | PV-08 | 原因が推定なので `EVIDENCE_REQUIRED` として特定してから行う |

- 参照 SVG が変わる PR は、`--update-reference-outputs` と Review REQUIRED で扱う。

**Phase 2（P2、約 25 PR。クラスと衝突するファイルでまとめる）**

| まとまり | PR 数 |
|---|---|
| エラー変換（X-01 系、SE-05、PV-04、IN-04 のメッセージ） | 3 |
| 数値の黙った補正（X-02 系） | 2 |
| History（SE-01、SE-02/03、GE-06、FE-04） | 4 |
| 引き継ぎ（GE-02、PV-01/02/03、CO-08。DP-2 の後） | 5 |
| 複数 Result（FE-02/03/05。DP-1 の後） | 3 |
| 入力（IN-04 の状態、IN-05、IN-06、IN-08） | 3 |
| Session/CLI（SE-06/07/08、GE-03） | 3 |
| 比較（CO-01、CO-06、CO-07） | 2 |
| FE-06、PV-05、TR-02/03 | 各 1 |

**Phase 3（P3 22 件、約 8 PR）** owner ファイル単位でまとめる。例：TR-10/11 の a11y と CSS、TR-12 の 1 行、FE-09/10/11 の editor など。

**Phase 4（仕上げ、約 3 PR）**
- "Result content commit" の allowlist を最終形まで狭める（縮小の経路を使う）。
- G-A と G-B を Product Impact map の contract に昇格させる。authority だけの PR で、`docs/internal/PRODUCT_IMPACT_RATCHET.md:617-628` の条件を満たした場合に限る。
- smoke の代表ケースを見直す。

**並べ方：**
- ホットなファイルは直列にする：`gbdraw/web/js/app/run-analysis.js`（CO-02/03/07/08、X-02、IN-05）、`app/app-setup.js`、`services/config.js`、`services/session-request.js`。
- ファイルが重ならないものは並列にしてよい。
- **マージは staging と同期させる**：
  - 同じ staging の周期で runtime のマージは 2〜3 本までにする。
  - dev の先端で Dev staging / gate が緑になってから、次の同じクラスの修正をマージする。
  - staging が赤になったらマージ列を止め、その周期でマージした範囲の中で原因を切り分ける。
- **自動マージ**：記憶にある 2026-09-29 の許可に従い、`gh pr merge --merge --match-head-commit` で dev にだけマージしてよい。main は対象外。ただし次の PR は自動マージから外すことを推奨する。
  - 理由：Review REQUIRED は exit 0 で、Gate を失敗させない（`docs/internal/WEB_CHANGE_POLICY.md:86-149`）。
  - 対象：Product Contract の同時変更、参照出力の変更、size 超過、architecture-change label、policy 文書。
  - これらは Owner のレビューを待つ。
- **決定パック**：standing instruction（推奨案で進める）を使う場合は、receipt に「Owner-delegated（2026-09-29）」と明記する。AGENTS.md の「製品の選択肢を自律的に選ばない」との緊張があるため、4 パックを一度提示して確認する運用を推奨する。

## 5. 横断的な提案（10）とリスク

**提案：**
1. 「正解の基準のない一貫性検査」を合格の根拠にしない。新しい編集機能のテストには、fresh Generate か CLI との比較を最低 1 つ含める。
2. `gbdraw/web/CLAUDE.md:235-252` の Changing a setting に 4 項目を加える：モードの行列、不正値の拒否、CLI での再現、即時編集と Generate の等価性。
3. 共通 helper の自動検査（G-G(3)）を最初に入れる。既存の 15 の spec がそのまま検査網になり、費用がもっとも小さい。
4. JS による Python の再実装（高速経路）は、共有ベクター（Python が正）がない限り受け入れない。
5. エラーは文字列照合から構造化された code へ移す方針とし（修正側の workstream）、移行が終わるまでは G-G(1)/(2) の allowlist を縮小だけの ratchet にする。
6. Result への書き込みを `flushActiveResult` 1 か所に集める（G-J(1)）。クラス 1 と 2 の方針を 1 か所で実装できるようになる。
7. 既知の不具合は `test.fail`、strict xfail、`knownDefect` で先に記録し、修正 PR が印を外す。
8. 「CW の計測は、欠けている metric を 0 と読まない」を helper で強制する（`session-regeneration-contract.cjs:10`、要確認）。
9. 規約文書の記述と実装のずれ（smoke 件数、rule 上限）を検出する小さな pytest を置く。
10. 定期監査は昇格前のチェック項目として置き、新しい gate や workflow は作らない。

**リスク：**
- **ブラウザテストの flaky**：
  - retries 2 が不安定さを隠す。JSON report の retries を監視し、新しい spec は時間待ちをせず、操作の完了を待つ。
  - `test.fail` と retries が重なると、たまたま通ったときに誤って失敗することがある。行列系の spec は `retries: 0` にし、1 ワーカーで順に実行する。
- **CI 時間の増加**：STG が約 15 分増える見込み（未計測）。shard が 45 分を超えないよう、導入時に計測し、超える場合は完全版を REL に移す。
- **テストが実装に合わせ込まれる**：
  - G-A と G-B は、DOM の構造ではなく意味比較（`svgVisualSemantics`）と利用者が持つ状態で判定する。
  - 無視する属性の一覧は gallery parity と共通にし、増やすときはレビューを必須にする。
- **等価性の厳しすぎ**：PV-10 などは、製品として「即時編集は近似でよい」と判断される余地がある。DP-2 より前に G-A を hard にしない。
- **ガバナンスの順序による停滞**：G-J(1) は authority → checker → runtime の 3 段になる。Phase 0 で先に入れないと、Phase 1 の修正が縮小の経路を使えない。
- **静的走査の誤検出**：G-G(1)/(2) の AST/正規表現の走査は、f-string やテンプレートの展開で誤検出しうる。allowlist は理由付きで、縮小しかできない形にする。
