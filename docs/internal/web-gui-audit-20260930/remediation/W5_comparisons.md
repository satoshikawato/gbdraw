<!-- Raw design report of workstream W5 (sub-agent output, 2026-09-30, base origin/dev 4c89bab1). Evidence for 01_REMEDIATION_PROPOSAL.md, not Product or architecture authority. -->

# W5 比較 — 修正提案（dev 4c89bab1、読み取りのみ）

## 前提と自分で集めた証拠
- パスはリポジトリ直下から、行番号は 4c89bab1 のもの。OIPC は `docs/internal/OPTION_INTEGRITY_PRODUCT_CONTRACT.md` を指す。
- 実験はすべて scratchpad の `design/w5/` で行い、dev は変更していない。
  - **E1**: dev の CLI で `SA_SB.tsv`（SB 1101..1300）を描いた。
    - `--reverse_complement 0 --reverse_complement 1` では、逆相補にした SB の record-local 位置 1101..1300（x 734..867、source では 1900..1701）に描かれた。
    - `--region SB:1001-2000` では、crop 後の record 端（x 666.7）を越えた x 734..867 に描かれ、警告も出ない。
    - つまり、表は crop と逆相補を適用した後の表示座標として読まれ、範囲も検証されていない（`gbdraw/linear_comparison.py:25-26` が `:31-33` の範囲検査より前に return する）。
  - **E2**: native LOSAT（`/mnt/c/Users/genom/GitHub/LOSAT/LOSAT/target/release/LOSAT`）で BGC0000708/709 の protein を比べた。
    - 対ごとに A×B を検索すると、A→B は 32 行（max_target_seqs 1）、232 行（既定）。
    - query 側だけをまとめた AB×B は、A→B の行が両条件とも完全に一致した。
    - subject 側もまとめた AB×AB（現在の Web の方式）は、0 行（limit 1、自己ヒット）と 148 行になり、E-value は約 1.8 倍（1.8→3.2、5.2→9.2）になった。
    - blastp と blastn に `-dbsize` / `-searchsp` はない。
    - 時間（native、平均）
      - 2 record: 1 job 886 ms、subject-record 2 job 1641 ms、対ごと 4 job 1312 ms。
      - 3 record: 1135 / 2523 / 2839 ms。
  - **E3**: `reverseComplementOverride` と `projectCommittedRecordTransform` は main（4556e04e）にない（dev のみ）。一方、source batch、`searchContext`、threaded の既定値、`canReuseResolvedProteinArtifacts` は main にもある。

## 共通の規則（提案）

### R-FRAME: 座標系
座標系を 3 つに分ける。

- **F（探索座標系）**: 選択と crop の後の配列で 1 始まり。元の鎖の向きで、逆相補と回転の前。LOSAT の出力そのもの。
- **V（表示座標系）**: F に実効の逆相補を適用したもの（x→L+1−x）。回転は含まない。回転は renderer の `project_match_endpoints`（`gbdraw/linear_comparison.py:15`）が適用する。
- **Src（source 座標系）**: 入力ファイル上の座標。feature popup はこれを使う。F→Src の写像は `RecordDisplayTransform`（`gbdraw/layout/record_coordinates.py`）が持つ。

現状は次のとおり。

| 対象 | 使っている座標系 | 場所 |
|---|---|---|
| raw cache と Save Raw | F | |
| 生成した LOSATN の request 行 | V | JS から helper `convert_losat_nucleotide_to_display_tsv`（`gbdraw/web/js/app/python-helpers.js:439-471`）へ渡す |
| upload、CLI `-b` | V（E1） | |
| protein 行 | V と view hash | planner が向きを投影し直す（`gbdraw/api/record_planning.py:216-300`） |
| popup と FASTA ヘッダ | V | |

提案する規則は次のとおり。

1. raw と Save Raw は F のまま。docs に「preserve raw search rows」とある（`docs/REFERENCE/comparison-programs-thresholds-and-results.md:150`）。
2. 保存する比較 evidence、upload、CLI の表はすべて F にする（推奨。ただし Q-FRAME の決定が必要。CO-07 を参照）。F→V の向きの投影は Python の planner だけが持つ。
3. 人が読む座標（popup、FASTA ヘッダ）は Src にする。
4. 決定が出るまでは、V を作る処理を 1 か所にし、向きは 1 つの resolver から読む（CO-03）。

### R-CACHE: cache identity
- **raw**: 向きのある record 対ごとに 1 つ。(program, outfmt, 正規化した args, query の配列/protein set の hash, subject の hash) で決まる。DB の範囲は subject の 1 record だけ（CO-04 A の場合）。向き、回転、ラベル、ファイル名、ファイルの分け方、行配置、スレッド数は含めない。
- **derived**: 既存の `buildLosatDerivedPayloadCachePayload`（`gbdraw/web/js/app/run-analysis.js:4180-4204`）が唯一の owner。
- **resolved の再利用判定**: 抽出しなくても分かる入力をすべて比べる。mutation イベントによる invalidate（`gbdraw/web/js/app/app-setup.js:508-514`）は最適化にとどめ、正しさはこれに頼らない。

---

## CO-01 Threaded の既定値で失敗する
1. **分類**: IMPLEMENT_EXISTING_AUTHORITY。
   - 根拠は PD-OI-018（OIPC:811-815「defaults Execution to `threaded`… existing execution modes and their support checks remain」）と、`docs/REFERENCE/web-app.md:19-22`（「Select **Serial**… when that capability is unavailable」）。
   - 結論: threaded は明示でも既定でも厳格なままにする。そのうえで、原因と対処の分かる型付きエラーを出す。
   - 自動で serial に落とす案は threaded の意味を変えるので、PD-OI-018 の改訂が必要になる。今は推奨しない。
2. **根本原因**: 監査のとおり。追加で 2 点ある。
   - LOSAT の実行中も `failureStage` が `'request-validation'` のまま変わらない（`run-analysis.js:1903`。変わるのは `:4675` の `'render'` だけ）。
   - 「serial への fallback がない」ことは、上の authority では不具合ではない。
3. **修正案**
   - `gbdraw/web/js/services/losat.js:101-109` の 3 つの判定を、同期関数 `losatThreadingPrecondition()` に切り出す。`getLosatThreadingSupport` と dispatch の両方がこれを使う。
   - `:792-793` の `throw new Error(...)` を、code `LOSAT_THREADING_UNAVAILABLE` を持つエラーに置き換える。定義と文言は W3（X-01）が担当する。
   - この判定は、未キャッシュの job がある dispatch の時点で行う。`run-analysis.js:3265` で先に判定すると、全件キャッシュ済みの再生成まで失敗させてしまうので不可。
   - LOSAT の実行に入る前に stage を LOSAT 用に切り替える。STAGES（`gbdraw/web/js/services/error-normalization.js:6-7`）に項目がないので、W3 と調整する。
   - `gbdraw/web/js/app/losat-settings.js:121-129` では、同じ判定を使って「Threaded（この環境では不可）」と表示する。
   - fallback、環境によって変わる既定値、新しい state は追加しない。
4. **上位の対応**: 不要。
5. **テスト**
   - node: `crossOriginIsolated=false` のとき判定が unavailable を返すこと。
   - Playwright: 既定の test server は分離されていない（`playwright.config.js:19`）。この上で、既定のまま LOSATN を Generate すると `LOSAT_THREADING_UNAVAILABLE` が出て、job 数が 0 で、前の Result が残ること。全件キャッシュ済みの再生成は成功すること。Serial では成功すること。
6. **規模・PR**: 本番ファイル 2〜3（W3 の定義を含む）、約 40 行。Ordinary、Review CLEAR。
7. **依存・リスク**: W3（X-01）。テストはすでに serial を指定している（77563b6b）。

## CO-02 LOSATP の古い結果の再利用（P1）
1. **分類**: IMPLEMENT_EXISTING_AUTHORITY。根拠は OIPC-C04（OIPC:304-310）、PD-OI-022（「Cache reuse still validates the actual inputs」）、`docs/REFERENCE/comparison-programs-thresholds-and-results.md:155-158`（selected features が変われば invalidate される）。
2. **根本原因**: 監査のとおり。再利用の仕組みは 2 つある。
   - 段階的な cache（protein 抽出のキーは visibility を含む。`run-analysis.js:3322-3332`。derived は `:4180-4204`）は正しい。
   - 近道の `canReuseResolvedProteinArtifacts`（`:468-541`、`:3069-3100` で使用）は、手書きの identity で段階的な cache を迂回する。visibility がこの identity に含まれていない。
   - `files.linearCanonicalComparisons` は成功のたびに作り直される（`gbdraw/web/js/services/config.js:3091`）。消えるのはイベントによる invalidate のときだけで、visibility の編集はそれを呼ばない。
   - この近道は、raw cache のない CLI や Gallery の Session で必要なので（PD-OI-008）、削除しない。
3. **修正案**
   - `active` に `featureVisibility` を加える（値は `:2353` で計算済みのもの）。
   - これを committed Session の `diagramOptions.featureVisibilityTableFile` のテキストと比べる。比較の前に、両方を `serializeFeatureVisibilityRules(parseFeatureVisibilityRules(text).rules)` で正規化する。
   - テキストの取得には `resourceTextFromRef`（`gbdraw/web/js/services/session-request.js:3283`）を使う。`session-request.js` が request 同値性の owner なので、小さな helper を export する。render request ではなく committed Session（resources を持つ）を読む。
   - 同じ述語で、source resource の同一性も比べるとより堅い（`matchesSessionResourceDescriptor` がある）。
   - 新しい cache 層や invalidate イベントは追加しない。
   - 不変条件: 抽出前に分かる入力がすべて committed と一致するときだけ再利用する。
4. **上位の対応**: 不要（R-CACHE を適用）。
5. **テスト**
   - `tests/web/run-analysis-derived-cache.test.mjs` の probe 方式で、visibility が違えば false、正規化後に同じなら true になること。
   - Playwright（監査の losatp-reuse を元にする）: Similarity、Pairwise、Collinear のそれぞれで、「visibility を変えて Generate」の結果が「cache を消した新規実行」と一致すること（例: 18 ribbon / 17 group）。hidden CDS を指す ribbon がないこと。Save→Load→Generate でも同じになること。
6. **規模・PR**: 本番ファイル 2、約 25 行。Ordinary。科学的出力が正しくなる変更なので Review REQUIRED。
7. **依存・リスク**: CLI の visibility ファイルは、正規化しないと不要な再実行が起きる。X-02 の閾値の正規化（W3）と関係する。

## CO-03 回転と Orient forward（P1）
1. **分類**: IMPLEMENT_EXISTING_AUTHORITY。根拠は PD-OI-032 の item 8（OIPC:1377-1381）と AC-14/15（:2389-2390）。upload の表を回転したときの扱いだけは Q-FRAME（CO-07）の決定に依存する。
2. **根本原因**
   - Linear の向きの保存先が 2 つある。
     - カードの `seq.region_reverse`（`gbdraw/web/index.html:1896`）。
     - rotate が書く draft の `reverseComplementOverride`（`gbdraw/web/js/app/record-display-options.js:205-233, 480`）。
   - request builder は両方を見る（`session-request.js:2318-2358`、`effectiveRecordReverseComplement` の `record-display-options.js:138-146`）。
   - 他の利用箇所は `region_reverse` しか読まない: `run-analysis.js:3177`（LOSAT の座標変換）、`:492`（再利用判定）、`gbdraw/web/js/app/match-sequences.js:868`、`gbdraw/web/js/app/feature-metadata-extraction.js:226`。
   - 誤るのは 2 つの経路である。
     - (a) 新規 Generate が reverse=false で F→V に変換する。
     - (b) Apply が `projectCommittedRecordTransform`（`session-request.js:4945-5031`）で向きだけを変え、nucleotide の V の行を変換し直さずに直接描画する。
   - alignment はすでに `region_reverse` に書き、override を null にしている（`record-display-options.js:264-277`）。rotate だけが例外である。
   - LOSATP は view hash による投影し直しで、どちらの経路でも正しいはず（監査の「CODE-ONLY」は「たぶん影響なし、要確認」と読むべき）。
   - 同じ原因による未確認の副作用: override が残っていると、カードの checkbox が効かず、表示も実際の向きと食い違う。
3. **修正案**
   - (a) Linear の向きの唯一の owner を `seq.region_reverse` にする（PD-OI-027 の record-owned orientation）。
     - `writeResolvedTransform` は、Linear 行で linearSeq と 1:1 に対応するとき、`region_reverse` に書いて override を null にする（`commitAlignmentOrientations` と同じ処理）。
     - capture と restore は `region_reverse` を含める（`captureAlignmentOrientationIntent` と同じ形。AC-10/12/17 のため）。
     - override は、Circular と、cardinality が `all` で複数行になる Linear（`session-request.js:2330-2340`）にだけ残す。後者のために `resolveLinearRecordReverse` を 1 つ export し、`run-analysis.js:3177, 492` と `match-sequences.js:868` はこれを使う。直接読んでいる箇所は削除する。
     - dev だけの形式なので（E3）、互換 reader は不要。dev 所有の fixture と Gallery を新しい形式で作り直す。
   - (b) Q-FRAME が決まるまでの暫定策
     - `projectCommittedRecordTransform` は、向きが変わるとき、target を端点に持つ生成済みの nucleotide 比較（`addResolvedComparisonResource` の `comparison-resolved-*`、`session-request.js:1826-1860`）について、target 側の座標を L+1−x に書き換えた text resource を candidate に入れる。crop は `:4966-4968` で拒否されるので、L は record 全長。
     - upload の比較を持つ record では、Orient forward を理由付きで不可にする（PD-OI-032 item 7、OIPC:1371-1375）。向きを変えない回転は許す。
     - Q-FRAME が C に決まれば、(b) も upload の制限も削除する。JS 側の F→V 変換（`run-analysis.js:3172-3208, 3478-3485, 4341-4377`）も削除できる。
   - 追加しないもの: 新しい向きの保存先、Apply 経路での LOSAT 再実行、Worker protocol。
   - 不変条件: request、LOSAT の投影、再利用判定、配列の materialize は、同じ 1 つの実効向きを読む。
4. **上位の対応**: R-FRAME を適用する。ownership 表に「Linear の向き: `seq.region_reverse`」と追記する。
5. **テスト**
   - Playwright の metamorphic テスト（R2c+R3c、LOSATN）。次の 4 経路で ribbon の端点がすべて同じで、期待値は R3 の表示座標 301..1300。
     - P1: カードの逆相補と start 1300。
     - P2: popup の 5′ + Orient forward で Apply。
     - P3: P2 の後に新規 Generate。
     - P4: Save→Load→Generate。
   - P2〜P4 では LOSAT job が 0 であること。P2 の後にカードの checkbox がチェック状態になること。Undo で向きと Result が一緒に戻ること。
   - unit: `tests/web/record-display-options.test.mjs`。`projectCommittedRecordTransform` を 2 回適用すると元に戻ること。変換は生成済みの行だけに適用されること。
   - LOSATP の回転と Orient forward を実際に動かして確認する。
6. **規模・PR**: 本番ファイル 4〜5、約 150 行、net 約 +70。Ordinary。Review REQUIRED（科学的出力、OIC-021 の証拠）。
7. **依存・リスク**: Q-FRAME の決定。`feature-metadata-extraction.js` の担当 workstream との調整。PD-OI-029 の Reset で行う後からの手動編集の判定が壊れないこと。

## CO-04 ファイルのまとめ方でリンクが変わる（P1）
1. **分類**: PRODUCT_DECISION_REQUIRED。PD-OI-018 item 7（OIPC:773-789）は source ごとの batch と「four directed source jobs」「must not become 64」を明記している。E-value と Max target seqs がまとめ方に依存するのは、この決定から直接出てくる。
2. **根本原因**: 監査のとおり（E2）。
   - subject 側に record を足すと、自己ヒットが Max target seqs を埋め、DB サイズも変わる。
   - query 側をまとめても結果は変わらない。
   - job 数の見積もり（`losat-settings.js:48-90`）が batch のキーを独自に計算しており、`prepareLosatSourceBatches` と二重になっている（PD-OI-020 の「表示と実行の一致」に反する危険）。
3. **選択肢**
   - **A（推奨）**: DB を subject の 1 record に限る。
     - A1: query 側は、その subject への要求がある同じ source の record だけをまとめる。
     - A2: 単純に対ごとに実行する。E2 では 2523 ms 対 2839 ms と差が小さい。
     - どちらでも Web は CLI（`gbdraw/analysis/protein_colinearity.py:6630-6760` は対ごと）と一致する。
     - 修正: `gbdraw/web/js/app/linear-sources.js:151-212` の batch キーを (querySource, subjectRecord, args) にする。`excludeSelfComparisons` と `searchContext` を削除する（`run-analysis.js:239-285, 3842, 3852, 3869, 3933, 4097` と、Python 側の cache key 生成）。純粋関数 `planLosatSourceJobs` を作り、見積もりにも使う（`losat-settings.js` の重複を削除）。
     - 不変条件: raw evidence は (program, args, query の内容, subject の内容) だけで決まる。
   - **B**: batch は維持する。要求されていない自己検索を全モードで除外し（0 対 19 の問題はこれで直る）、E-value が source DB に依存することを docs と Run Info に明記する。
   - **C**: limit の有無で方式を切り替える。挙動が 2 通りになるので非推奨。
   - A を推奨する理由: PD-OI-018 の rationale自体（「File count is not record count」）が、ファイルのまとめ方に依存しないことを求めている。
4. **上位の対応**: PD-OI-018 を revision 4 に改訂する。docs に「E-value は record 対ごと」と書く。
5. **テスト**
   - `tests/web/linear-sources.test.mjs:196-243` を更新する。subject が常に 1 つであること、要求されていない (i,i) がないこと、まとめ方の違う 2 通りで per-pair の identity が同じであること。
   - Playwright: multi.gbk と別ファイルで、各対の raw TSV が一致し、ribbon も一致すること（limit 1 / 上限なし、全モード）。
   - wasm（serial/threaded）で `playwright.perf.config.js` による計測を行い、PD に添付する。
6. **規模・PR**: 本番ファイル 3〜4、約 200 行、net ≤0。authority を先に通す PR、または reviewed co-change。
7. **依存・リスク**: 性能（native で 2〜2.5 倍）。main の Session には `searchContext` 付きの cache が残るので、reader は拒否せず無視する。PD-OI-020/022/051 を守る。

## CO-06 upload の表の ID を照合しない
1. **分類**: PRODUCT_DECISION_REQUIRED。docs（`docs/REFERENCE/comparison-programs-thresholds-and-results.md:90-93`「They must map to the intended displayed records」）は利用者の義務を書いているだけで、強制するかどうかは決めていない。一方、metadata の矛盾は IMPLEMENT で直す。
2. **根本原因**
   - 監査のとおり。端点は edge の index から決まり（`session-request.js:2071-2081`）、行は照合されない。
   - `gbdraw/render/groups/linear/pairwise_match.py:436-445` は ID に表の値を優先し、index は端点から取るので、両者が矛盾する。
3. **選択肢**
   - A: 厳格。照合できない行はすべて拒否する。
   - **B（推奨）**: 矛盾はエラー、不明は警告にする。
     - 入れ替わった端点を指す行、または他の record を指す行はエラーにする。
     - どの record にも一致しない ID（contig_A など）は、位置で割り当てて警告を出す。
   - C: metadata だけ直す。
   - B の修正: 照合は Python の planner、`project_source_bound_comparisons` の隣（`gbdraw/api/request_render.py:1582` から呼ばれる）に置く。こうすると CLI、typed API、Web がすべて同じ経路を通る。
     - alias は `record.id`、version 接尾辞を除いた id、`record.name`。
     - 対象は feature binding を持たない表だけ。
     - `data-*-record-id` は常に端点の record から取る。
     - auto-swap と、行を他の edge へ振り分ける処理は追加しない。
4. **上位の対応**: `docs/REFERENCE/input-formats-and-tsv-schemas.md:28-31` と CLI reference に照合の規則を書く。
5. **テスト**
   - pytest: 入れ替わった表がエラーになること、version 接尾辞を許すこと、不明な ID で警告が出て metadata が一致すること。CLI と codec の両経路で確認する。
   - browser: 具体的なエラーが出て前の Result が残ること。contig_A/B でも FASTA 操作が有効になること。
6. **規模・PR**: Python 2 ファイル、約 80 行。Review REQUIRED（CLI の挙動と SVG metadata が変わる。参照 SVG を確認）。
7. **依存・リスク**: W6 の CO-05（同じ reader）の後に行う。

## CO-07 Save Raw と逆相補
1. **分類**: PRODUCT_DECISION_REQUIRED（Q-FRAME）。export 自体は docs の契約（raw の行）に従っている。欠けているのは、逆相補にした record に対して upload や CLI の表をどの座標系で読むかの定義である。`docs/REFERENCE/web-app.md:335-336` は、逆相補で「comparison endpoint mapping」が変わる、つまり表が record に追従すると読める。
2. **根本原因の再確認**: 監査より範囲が広い。
   - 利用者自身の順方向の BLAST 表に逆相補を適用しても誤る（CLI の E1、Web は監査の再 upload で確認済み）。
   - crop したとき、範囲外の行が record の外に描かれる（E1、検証なし）。
3. **選択肢**
   - A: 表は V のまま。Save Raw は、逆相補にした端点について helper で V にしてから出力する。
   - B: 表は V、raw は F のまま。座標系を画面と docs に明示する。再現には Run Info のファイル（V、`run-analysis.js:4362`）を使う。
   - **C（推奨）**: すべての比較表を F にする。
     - planner（`record_planning.py:216-300`）が、feature binding を持たない表を `collection.transforms` の source_step で反転し、1..L を検証する。
     - JS 側の F→V 変換と helper（`python-helpers.js:439-471`）を削除する。CO-03 (b) も不要になる。
     - 互換: main でも保存された request の nucleotide 行は V なので、Session reader で逆相補の端点だけ V→F に一度変換する。upload の bytes も書き換わるので、これも明記する。main 由来の fixture を用意する。
     - CLI `-b` + `--reverse_complement` は後方非互換になるので、release note を書く。
   - D: source 絶対座標。crop した領域から作った表が壊れるので非推奨。
   - どの選択肢でも、範囲外の行は検証する（OIPC-C03）。
4. **上位の対応**: R-FRAME を新しい PD-OI として記録し、docs 4 か所（input-formats、comparison-programs、CLI reference、web-app.md:335）を更新する。
5. **テスト**
   - 向き 4 通りで、Save Raw → 再 upload / CLI で ribbon が一致すること。
   - C では、upload して逆相補を切り替えた結果が、LOSAT の結果と一致すること。
   - 範囲外の行で警告またはエラーになること。
   - migration の fixture。
6. **規模・PR**: C は Architecture（約 350 行、net ≈0）、A は 1 ファイル約 20 行、B は docs のみ。Review REQUIRED（科学的出力、互換）。
7. **依存・リスク**: CO-03、CO-10。Gallery を再生成する必要がある。

## CO-08 Similarity group の名前
1. **分類**: 別の group へ移る点は IMPLEMENT。残らない override を黙って消す点は OIPC-C05/C06（OIPC:312-326）に反する。消えた group の override をどう扱うかは PRODUCT_DECISION_REQUIRED。
2. **根本原因**: 監査のとおり。`og_N` は最小 member の順番で振られる（`protein_colinearity.py:5287-5294`、`:4876`）。prune（`run-analysis.js:1452-1461`、`:4934`）は、同じ ID が残っていれば別の group でも名前を残す。
3. **選択肢**
   - **A（推奨）**: member 集合が完全に一致する group へ付け替える。member の handle `h_…` は決定的なので、これで照合できる。一致しないものは dormant として webEdits に保存し（`orthogroupDormantOverrides`、無ければ []）、同じ group が再び現れたら戻す。一覧と Clear を画面に出す。
   - B: 完全一致で付け替え、一致しないものは Generate の transaction 内で通知付きで削除する（Undo で戻る）。
   - C: 重なりの割合で引き継ぐ。誤った付与になるので非推奨。
   - 修正: `pruneOrthogroupOverrides` を、`orthogroups.js` の `rekeyOrthogroupOverrides(previousGroups, candidateGroups)` に置き換える。現在の override を読む側は変更しない。
4. **上位の対応**: Session に項目を 1 つ追加する（current writer。absent は空なので migration は不要）。
5. **テスト**
   - node: og_16 の {a,b} を改名し、組み直して {a,b} が og_17 になったら、名前が og_17 に付くこと。
   - 消えた group の説明が dormant に入り、閾値を戻すと復元されること。
   - Playwright og-rename。
6. **規模・PR**: 本番ファイル 3〜4、約 100 行。Review REQUIRED（Session 項目）。
7. **依存・リスク**: history-snapshot に dormant を含める。

## CO-10 popup の座標
1. **分類**: FASTA ヘッダは IMPLEMENT（出所の正確さ。docs は「genomic spans」と書いている）。popup の表示座標は PRODUCT_DECISION_REQUIRED。
   - A: Src だけ表示する。
   - **B（推奨）**: Src を主に表示し、値が違うときだけ「表の座標」も示す。
   - C: 現状のまま、ラベルだけ付ける。
2. **根本原因**
   - `data-qstart` などは V（`pairwise_match.py:450-454`）で、それが popup（`gbdraw/web/js/app/pairwise-match-popup.js:295-321, 1411-1413`）とヘッダ（`match-sequences.js:754-760`）に出る。
   - 配列そのものは正しい（`match-sequences.js:860-887`）。
   - standalone SVG 側に同じ処理の複製がある（`gbdraw/web/js/services/standalone-interactivity-assets.js:5231, 5278`）。
   - SVG に Src への写像の情報はない（`gbdraw/render/groups/linear/seq_record.py:389-390` は id と index だけ）。
3. **修正案**
   - renderer が record group に `data-gbdraw-record-source-base` / `-source-step` を出す。値は `RecordDisplayTransform` から取る。
   - 写像の関数は 1 つにし、ヘッダと popup で使う。standalone 側も同じ式を使う。
   - C の後なら Src = base + F − 1 で済む。
4. **上位の対応**: `docs/SVG_SEMANTIC_HOOKS.md` を更新する。
5. **テスト**: `tests/web/match-sequences.test.mjs` で、crop SB 1001..2000 のときヘッダが coords=1001..2000 になること。逆相補した record の場合。standalone との一致。Python で属性値の一致を確認する。
6. **規模・PR**: Python 1 + web 3 ファイル、約 70 行。Review REQUIRED（参照 SVG が変わる）。
7. **依存・リスク**: Q-FRAME の後に行う。

## FE-07（W6 の担当）が比較に与える影響
影響はある（コードからの判断）。
- LOSATP は、`/translation` のない CDS（GFF3 はすべて）について `_translate_cds_feature`（`protein_colinearity.py:2761-2774`、呼び出しは `:2885-2890`）を使う。これは Web と CLI で共通である。
- GTG/TTG の開始コドンで先頭が V/L になるので、bitscore、E-value、identity がわずかに変わる。閾値の近くでは ribbon や group の membership が変わりうる。
- W6 の修正は次を満たす必要がある。
  - 翻訳処理を 1 つにし、`gbdraw/web_support/feature_metadata.py:~290` と `_translate_cds_feature` の両方で使う。
  - M にするのは、5′ 端が完全（`<` がなく、codon_start が 1）で、その table の start codon の場合に限る。
  - cache は content-hash なので、変更に合わせて自動で外れる。schema の bump は不要。
  - GFF3 の LOSATP の参照出力が変わりうる（Review）。
  - CLI と Web の両方で、先頭が M になるテストを追加する。

## 横断的な提案
1. Q-FRAME の Decision Pack を最初に出す（C を推奨）。CO-03 (b)、CO-07、CO-10、upload + 逆相補の問題がこれで一度に決まる。
2. R-CACHE: identity はデータとして比べ、正しさをイベントに頼らない。LOSAT の job 計画は 1 つの純粋関数にし、実行と見積もりの両方で使う。
3. Linear の向きの保存先を 1 つにする。`region_reverse` を直接読む箇所を owner の外で禁止する fitness test を追加する。
4. すべての比較行について、端点が record 内にあることを planner で検証する（`gbdraw/linear_comparison.py:25-33`）。
5. metamorphic テスト群を作る: ファイルのまとめ方、向きを変える経路、visibility の再利用と新規実行、Web と CLI の raw 行の一致。
6. LOSAT 専用の失敗 stage を用意する（W3 と調整）。
7. standalone popup の重複は別の課題として記録する。今回の PR では統合しない。
8. 順序: CO-01、CO-02 → CO-03 (a) → 各 PD（Q-FRAME、CO-04、CO-06、CO-08、CO-10）→ 実装。
9. 要確認: wasm の LOSAT で、Max target seqs を空欄（上限なし）にしたときに既定の 500 が隠れた上限として働くか（native の既定値は 500）。PD-OI-001 の「No hidden cap」との関係。

## 分類の誤り、不具合ではないもの、新たな発見
- **CO-01**: fallback がないのは契約どおり。実際の問題は UNKNOWN の表示（W3）と stage の誤り。
- **CO-03 の LOSATP**: たぶん影響なし（view hash で投影し直される）。
- **CO-04**: PD-OI-018 の帰結であり、PD が必要。
- **CO-06**: 位置による割り当ては設計どおり。不具合は metadata の矛盾で、照合をするかどうかは PD。
- **CO-07**: export は契約どおり。本質は座標系が決まっていないことで、P1 相当の広い問題（下の N2）。
- **CO-10**: 問題は FASTA ヘッダで、popup の座標は製品判断。
- **新規**
  - N1: crop したとき範囲外の行が record の外に描かれる（CLI で確認）。
  - N2: upload + 逆相補で ribbon を誤る（CLI で確認）。
  - N3: override がカードの checkbox を上書きする（コードからの判断）。
  - N4: job 数の見積もりが二重に実装されている。
  - N5: LOSAT の失敗が request-validation として報告される。
- **未検証**: LOSATP の回転 + Orient forward の実動作、wasm の max_target_seqs の既定値、Circular conservation の座標系。
