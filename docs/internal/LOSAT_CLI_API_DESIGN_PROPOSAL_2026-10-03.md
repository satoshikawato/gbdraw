# CLI / Python API から LOSAT を実行する設計案

- 日付: 2026-10-03
- 状態: 提案（Owner の判断待ち。コードは変更していない）
- 調査の基準: `origin/dev` `be1bbd46`（2026-10-03 fetch）、`origin/main` `4556e04e`、release tag `0.13.0`
- 出典の書式: `path:line` は `origin/dev` の path。`main:path:line` は `origin/main`、`0.13.0:path:line` は tag `0.13.0`
- 依頼（Owner、2026-10-03）: TLOSATX、LOSATN、similarity ring も CLI / API から LOSAT で直接実行できるようにする。`--protein_blastp_mode` は名前が分かりにくい。CLI / API は大きく変えてよい。実装は次の session で行う。

---

## 0. Owner の決定（2026-10-03）

この節は本文より優先する。ここに無い D1〜D3、D5、D6、D8、D10、D13〜D17 は推奨のまま採る（Owner-delegated、§5）。

| # | 決定 | 本文への影響 |
|---|---|---|
| D4 | **拒否する。** 旧 CLI flag は fresh run で置き換え先を示して拒否し、legacy session（version 27–30）の argv は書き換える | 推奨どおり |
| D7 | **Web に合わせる。** CLI / API の LOSATP も source-file batching の E-value scope にする | PR-5 を実施する。複数 record の file では CLI の LOSATP の E-value が変わる。Gallery の差分を review する |
| D9 | **固定しない。** track している `gbdraw/bin` の binary を v0.1.0 に置き換えない | PR-1 から binary の置き換えを外す。runtime の解決順は今のまま。結果の同等性は版の一致に依存するので、使った runtime の識別（path、版、program）を Session と Run Info に記録する（D10）。CLI Reference の「Linux x86_64 用に同梱」という誤記は B16 として別に直した |
| D11 | **固定しない。** Web の WASM を v0.1.0 から作り直さない | PR-6 から WASM の pin を外す。`dc-megablast` は Web に残す。CLI で native の runtime が task に対応しないときは、`COMPARISON_INPUT` の診断付きの error にする |
| D18 | **`--conservation_sequence` を新設する（Owner 承認 2026-10-03）。** `--conservation_fasta` は D4 の規則で退役させ、Python の `conservation_fasta_files` は `conservation_sequence_files` に改名する | 「D12 の追加の設計」の option の名前の項を確定とする |
| D12 | **ring は LOSATN と TLOSATX。比較する genome の入力は FASTA に加えて GenBank / DDBJ の flat file も受け付け、parse して配列を取り出す。CLI、Python API、Web のすべてで同じ。** | 下の「D12 の追加の設計」を見る |

### D12 の追加の設計

- **読み手は 1 つ。** 比較する genome の file から配列を取り出す処理は Python の `gbdraw.io` に 1 つだけ置き、CLI、Python API、Web（Worker の Python helper）が同じものを使う。形式は内容で判定する（`>` で始まれば FASTA、`LOCUS` で始まれば GenBank / DDBJ の flat file）。DDBJ の flat file は GenBank と同じ形式なので、同じ parser（Biopython）で読む。JavaScript に 2 つ目の配列の reader を作らない（R9、G-J(2) と同じ考え方）。
- **1 file = 1 genome。** multi-record の file は全 record を 1 つの genome（LOSAT の query の集合）として扱う。E-value の scope も file 単位で、D7 と同じ規則にする。
- **配列の無い record は拒否する。** `ORIGIN` が空、または `CONTIG` だけの GenBank は `INPUT_UNREADABLE`（reason は新設、例 `SEQUENCE_MISSING`）で止める。
- **TLOSATX の gencode。** GenBank / DDBJ から読んだ場合も既定値は 1（D17）。`/transl_table` からは推定しない。`--conservation_losat_gencode` で指定する。
- **ring の label の既定値。** GenBank / DDBJ では record の definition か organism、FASTA では file 名（今の precomputed の規則に合わせる。実装時に今の既定値を確かめる）。
- **Web。** ring の行に GenBank / DDBJ の file を置けるようにする（今の FASTA の入力欄と同じ場所）。Session には今の FASTA と同じく resource として保存する。
- **option の名前（D18、Owner 承認済み）。** 入力が FASTA に限られなくなるので、`--conservation_fasta` の名前は内容に合わない。**推奨: `--conservation_sequence FILE [FILE ...]`（FASTA、GenBank、DDBJ を受け付ける）を新設し、`--conservation_fasta` は D4 の規則で退役させる**（fresh run では置き換え先を示して拒否し、legacy session の argv は書き換える）。precomputed の ring（`--conservation_blast` と組み合わせる場合）も同じ名前にする。Python API の `conservation_fasta_files` も `conservation_sequence_files` に改名する（D5 と同じく alias なし）。D16 の「`--conservation_*` は今回変えない」はこの 1 つについて上書きする。
- **影響する PR。** PR-4 に含める（CLI、Python API、共通の reader、test）。Web の入力欄は PR-6 に含める。先に失敗させる test: GenBank と DDBJ の file を比較する genome にした LOSATN の ring が、同じ配列の FASTA から作った ring と byte 一致すること。配列の無い GenBank を拒否すること。

---

## 1. 目的と結論（推奨案の要約）

### 目的

Web が実行できる 3 つの LOSAT program（LOSATN = `blastn`、TLOSATX = `tblastx`、LOSATP = `blastp`）と Circular similarity ring を、CLI と Python API からも同じ意味で実行できるようにする。同じ入力と同じ設定なら Web と同じ hit を得る。名前は Web の UI label に揃える。

### 推奨（案 A）

1. **概念を 3 つに分ける。** (a) 検索 program は `--losat {losatn,tlosatx,losatp}`。Web の **LOSAT Mode** ボタンと同じ名前にする。(b) 検索する record pair は既定で Web の **All adjacent** と同じ。明示 pair は既存の `--comparisons_table` に `source` 列を足して表す。Circular は displayed record と各 comparison FASTA の組。(c) 表示は `--losatp_mode {similarity_groups,collinear,pairwise}`。Web の **LOSATP mode** メニューと同じで、LOSATP 専用。
2. **precomputed table は同じ comparison 概念の別 source とする。** `-b` と `--comparisons_table`、`--conservation_blast` と `--conservation_table` は変えない。LOSAT の結果は planner の中で同じ table（`LinearComparison`、conservation DataFrame）に解決する。renderer、track slot、threshold の扱いはそのまま使う。
3. **`--protein_blastp_mode` などの古い名前は「Retired inputs」として扱う。** fresh run では置き換え先を示して拒否する。0.12/0.13 の legacy session（version 27–30）の argv は既存の canonicalizer で書き換える。persisted wire 名（`generatedProteinComparison`、`losatpBin` などの settings key、`mode: "orthogroup"`、Web state 名）は変えない。新しい compatibility path は作らない。
4. **LOSAT runtime の owner を 1 つにする。** `gbdraw/comparisons/losat_runtime.py` を新設し、`analysis/protein_colinearity.py` にある runtime 解決、argv 構築、subprocess 実行を移して program 共通にする。NCBI BLAST+ fallback も program ごとに扱う。job plan は Web の `planLosatSourceJobs` と同じ規則を Python に実装し、shared vector で両者を固定する。E-value の database scope も Web と揃える。
5. **Session は Web と同じ形で書く。** LOSATN/TLOSATX の結果は `nucleotideBlast` resource と `losatCache`（schema 2、Web と同じ key）として保存する。replay に LOSAT は要らない。
6. **科学的な同等性は実測済み。** pin 済みの native LOSAT 0.1.0 は、Web 由来の tutorial table を再現する。LOSATN は byte 一致、TLOSATX は row 集合が一致し、描画した SVG は byte 一致した。mtDNA の TLOSATX 3 表は byte 一致した（§2.5）。

実装は 6 本の PR に分ける（§4）。Owner の判断が要る点は 17 個ある（§5）。

---

## 2. 現状（事実と出典）

### 2.1 Web

**Program と UI label**

- 実行できる program は `blastn`、`tblastx`、`blastp` の 3 つ（`gbdraw/web/js/services/losat.js:17`）。それ以外は例外になる（`losat.js:155-156`）。
- state は `losatProgram = ref('blastn')`（`gbdraw/web/js/state.js:92`）。
- UI label の対応は `LOSATN: 'blastn'`、`LOSATP: 'blastp'`、`TLOSATX: 'tblastx'`（`gbdraw/web/js/app/comparison-ui.js:20-24`）。
- LOSATP mode の対応は `PAIRWISE: 'pairwise'`、`SIMILARITY_GROUPS: 'orthogroup'`、`COLLINEAR_BLOCKS: 'collinear'`（`comparison-ui.js:26-30`）。
- 初期値は `losat.blastp.mode: 'orthogroup'`、`blastn.task: 'megablast'`、`executionMode: 'threaded'`（`gbdraw/web/js/services/session-active-config-contract.js:50-55`）。
- LOSATN task の選択肢は `megablast` / `blastn` / `dc-megablast`（`gbdraw/web/index.html:2219-2224`）。

**Linear の plan**

- state は `linearComparisonPlan = {mode, defaultSource, edges}`。
- `mode` は `none` / `adjacent` / `selected`、`defaultSource` は `losat` / `upload` で、既定は `none` と `losat`（`gbdraw/web/js/app/linear-comparisons.js:1-10`）。
- `adjacent` は隣接する row の record をすべて組にする（`adjacentRowPairs(..., true)`、`linear-comparisons.js:155-173, 268-280`）。
- `selected` は明示した edge を使う。edge は隣接 row の間に限られ、同じ row の中や隣接しない row の間は issue になる（`linear-comparisons.js:182-231`）。
- program は 1 回の Generate につき 1 つだけ。
- Similarity groups と Collinear は `adjacent + losat` で、record が 2 つ以上のときだけ動く（`linear-comparisons.js:323-327`）。selected と LOSATP の非 pairwise mode を組み合わせると `selected-losat-requires-pairwise` になる（`linear-comparisons.js:302-310`）。

**Job の展開** — `buildLosatJobSpecs`（`linear-comparisons.js:634-690`）

- orthogroup: 全 record pair を両方向で検索し、self も含める。
- collinear: inference が ON のときだけ self を含める。scope が `all` なら全 pair、`adjacent` なら各 edge を両方向で検索する。

**実行の単位** — `planLosatSourceJobs`（`gbdraw/web/js/app/linear-sources.js:170-213`）

- source file を 1 つの genome として扱う。
- 別々の file の間では、query source から subject file 全体を検索する。
- 同じ file の中では、query record を除いた file を検索する。
- self は要求されたときだけ実行する。
- 複数 record からなる batch では、`searchContext` を cache key に入れる（`linear-sources.js:252-254`）。
- Web の CLAUDE.md は、この関数を「the single LOSAT job plan」として、Generate と Settings の job 数見積もりの両方に使うと定めている（`gbdraw/web/CLAUDE.md:259-261`）。

**LOSAT の引数** — `buildLosatArgs`（`gbdraw/web/js/app/run-analysis.js:3485-3499`）

- blastn: `--task <task>`。
- tblastx: `--query-gencode <query record losat_gencode> --db-gencode <subject record losat_gencode>`。
- blastp: Pairwise mode のときだけ `--max-hsps-per-subject 1`。candidate limit が null でなければ `--max-target-seqs N`。
- e-value や filter は渡さない。threshold は検索の後に掛ける（`gbdraw/web/js/mode-profiles.js:114-136`）。
- `losat_gencode` は record ごとの入力値で、既定は 1（`state.js:136`）。`/transl_table` からの推定はしない（`gbdraw/web/js` に `transl_table` は 0 件）。

**Thread**

- blastp 以外は 1 job あたり 1 thread に固定し、job の並列で速くする（`losat.js:232`、`run-analysis.js:1553-1554`）。

**Circular ring**

- Web は Circular ring にも LOSAT を実行する。program は blastn と tblastx だけ（`run-analysis.js:892-895`、`index.html:2748-2749`）。
- query は各 comparison FASTA。subject は displayed circular の入力 file で、全 record を連結する（`run-analysis.js:2643-2669`）。
- 結果の reference 側は `'subject'` に固定する（`run-analysis.js:2797`）。
- tblastx の gencode は ring ごとの `losat_gencode` を `--query-gencode` に、`subject_gencode`（UI label は "Reference gencode"）を `--db-gencode` に渡す（`run-analysis.js:2630-2641`）。ring 行の UI label は "Subject gencode" だが、値は query 側に渡る（`index.html:2843-2845`）。
- 初期値は `{source:'losat', losat_program:'blastn', subject_gencode:1, reference:'auto', ...}`（`session-active-config-contract.js:57-58`）。

**Threshold**

- 既定値は mode profile が持ち、Web に生成される。Linear は `evalue=1e-2, bitscore=50, identity=0, alignment_length=0`、Circular は `1e-5, 50, 70, 0`（`gbdraw/mode_profiles.py:130-152`）。

**Cache**

- nucleotide の raw key は `{cacheSchema:2, program, outfmt, args, queryCanonicalHash, subjectCanonicalHash, flow?, searchContext?}` の SHA-256（`run-analysis.js:235-279`）。
- protein の raw key は Python の `build_protein_losat_cache_key` が作る（`gbdraw/analysis/protein_colinearity.py:1447-1462`）。
- derived cache は Web 専用。上限 16 entry で、mode ごとに持つ（`run-analysis.js:204, 329-427`）。
- orthogroup と collinear は同じ args なので raw entry を共有する。pairwise は `--max-hsps-per-subject 1` が付くので共有しない。
- どの key にも LOSAT / WASM の version は入っていない。

**Python への受け渡し**

- request の形は `renderRequest.schema = 8`（`gbdraw/web/js/services/session-request.js:161`）。
- comparison の kind は `nucleotideBlast`、`precomputedProteinComparison`、`orthogroupResult`、`collinearityResult`、`generatedProteinComparison` の 5 種（`session-request.js:1893-2159`、`gbdraw/session_request_codec.py:3361-3496`）。
- LOSATN/TLOSATX の結果は、search frame の table として `nucleotideBlast` resource に入れて送る（`run-analysis.js:4126-4152`）。

**Source recipe / Run Info**

- LOSAT の設定を CLI の flag に変換しない。生成した TSV を `-b` または `--comparisons_table` として出力する（`gbdraw/web/js/app/run-info.js:1353-1416`）。
- LOSATP の recipe は作れない。`nucleotideBlast` 以外の kind は `Source recipe unavailable` になる（`run-info.js:538-544`）。
- `LOSAT_DATABASE_SCOPE_NOTE` は「The CLI searches each record pair separately」と明記している（`run-info.js:1690-1692`）。

**WASM**

- `gbdraw/web/wasm/losat/{losat.wasm,losat-threaded.wasm}` は version で pin されていない。build 手順は README だけにある（`gbdraw/web/wasm/losat/README.md:6-20`）。
- 最後の更新は `6178f53d`（2026-05-25）。threaded 版には旧 CLI v1 の help 文字列が入っている。

### 2.2 Python core

**Typed options** — `gbdraw/api/options.py`

- `LinearDiagramOptions`（`options.py:980-1018`）:
  - `protein_blastp_mode: Literal["none","pairwise","orthogroup","collinear"]="none"`
  - `protein_comparison_pairs`
  - `losatp_bin="losat"`、`ncbi_blastp_bin=None`、`losatp_threads=None`
  - `protein_blastp_max_hits=5`、`protein_blastp_candidate_limit=None`
  - `orthogroup_member_max_hits=None`
  - `collinear_infer_orthogroups=True`
  - `collinearity_*`
  - `comparison_table_file`
- `CircularDiagramOptions`（`options.py:864-881`）: `conservation_blast_files`、`conservation_fasta_files`、`conservation_dataframes`、`conservation_reference="auto"`、labels、colors、ring width/gap、`conservation_table_file`。

**入門 API** — `gbdraw/interface.py`

- `LinearComparisonOptions`（`interface.py:333-383`）は、すでに分かりやすい名前を使っている（`protein_mode`、`losat_executable`、`blastp_executable`、`threads`、`max_hits`、`candidate_limit`、`pairs`）。これらを typed の名前に写す（`interface.py:904-926`）。
- `ComparisonRingOptions` / `ComparisonRingTrackOptions(source, label, color, comparison_sequence_source)`（`interface.py:305-322`）。`ConservationOptions` はその alias（`interface.py:326-327`）。

**Request と similarity alignment**

- `LinearDiagramRequest.similarity_alignment: SimilarityAlignmentPlan | SimilarityAlignmentReference | None`（`gbdraw/api/requests.py:707-725`）。Reference は `protein_blastp_mode == "orthogroup"` を要求する（`requests.py:748-756`）。
- planner は analysis を 1 回だけ実行して Reference を plan に解決し、request を `protein_blastp_mode="none"` と precomputed artifact に書き換える（`gbdraw/api/request_render.py:1572-1573, 1629-1825`）。この「intent を planner が解決済み artifact に置き換える」前例を、nucleotide 検索にも使う。

**LOSAT runtime** — すべて `analysis/protein_colinearity.py` の中にあり、`blastp` に固定されている

- runtime 解決 `_resolve_protein_blastp_runtime`（`protein_colinearity.py:3286-3341`）の順序:
  1. 明示した `losatp_bin`
  2. 明示した `ncbi_blastp_bin`
  3. conda の `$PREFIX/bin/losat`
  4. managed（`gbdraw/losat_setup.py:135-146`）
  5. bundled（source checkout だけ。`protein_colinearity.py:3146-3201`）
  6. PATH の `losat`
  7. PATH の `blastp`
- argv: LOSAT は `[exe,"blastp","-query",q,"-subject",s,"-outfmt","6",(-max_hsps),(-max_target_seqs),(-num_threads)]`、NCBI は subcommand を除いた同じ引数（`protein_colinearity.py:5885-5937`）。
- cache の args は Web の v1 表記（`--max-hsps-per-subject`、`--max-target-seqs`）を保つ（`protein_colinearity.py:6064-6074`）。
- `LosatpCacheManager` は schema 4 の raw cache（`protein_colinearity.py:2230-2581`）。
- fresh な CLI run は空の artifact から始めるので、毎回検索する（`gbdraw/linear.py:1515-1519`）。session replay は `losatCache` を再利用する（`gbdraw/api/session_compat.py:265-276`）。
- Python から blastn / tblastx を実行する path はない。nucleotide の cache entry は通過させるだけ（`gbdraw/session_io.py:1262-1267`）。

**`setup-losat`** — `gbdraw/losat_setup.py`

- lock は `gbdraw/data/losat-release.json:2-38`。version `0.1.0`、candidate `6bfb1b09…`、4 target（Linux x64 glibc≥2.34、Windows x64、macOS arm64/x64。`losat_setup.py:27-48`）。
- HTTPS で download し、size と sha256 を確かめて atomic に配置する（`losat_setup.py:215-254`）。
- cache の identity と `--version` を確認する（`losat_setup.py:110-132`）。
- CLI に option はない（`losat_setup.py:257-260`）。

**Bundled binary**

- `gbdraw/bin/linux-x86_64/losat` は git で track されている（sha256 `3d357ae9…`、2,789,728 B）。lock の binary（`fb01e0f6…`）とは別の dev build。
- wheel と sdist からは除外されている（`gbdraw/_build_support.py:119`、`MANIFEST.in:18`）。
- 0.12.0 と 0.13.0 では wheel に入っていた（`0.13.0:gbdraw/_build_support.py:30-31`）。

**Circular ring**

- 読み込みと threshold: `gbdraw/analysis/conservation.py:143-241`。
- reference 側の判定（query / subject / auto）: `conservation.py:281-313`。reference 側の record ID と一致する行だけを使う（`conservation.py:346-347`）。
- slot は `sequence_conservation` で、`conservation_<n>` として最初の features slot の内側に入る。ring は source の順に並ぶ（`gbdraw/api/diagram.py:519-568`）。

**Error**

- `GbdrawError(..., diagnostic={code, reason?, field?, ...})`（`gbdraw/exceptions.py:8-25`）。
- Web adapter の code の whitelist は `INPUT_INVALID INPUT_UNREADABLE DEPTH_INVALID TABLE_INVALID COMPARISON_INPUT TRACK_LAYOUT`（`gbdraw/web_support/error_adapter.py:45-53`）。
- LOSAT runtime と setup の error は diagnostic を持たない、ただの `ValidationError`（`protein_colinearity.py:3247-3283, 5964-6042`）。

### 2.3 CLI（現行の option）

**Linear** — `gbdraw/linear.py:259-381`

- `--comparisons_table`、`-b/--blast`
- `--losatp_bin`、`--ncbi_blastp_bin`、`--losatp_threads`
- `--protein_blastp_mode {none,pairwise,orthogroup,collinear}`
- `--protein_blastp_max_hits`（既定 5）、`--protein_blastp_candidate_limit`（既定 none）
- `--protein_blastp_output`
- `--align_orthogroup_feature`
- `--collinear_*`（うち 3 つは help を隠している）
- threshold は `linear.py:460-479`。

**Linear の検証** — `linear.py:766-788`

- table と `-b` は併用できない。
- mode と `-b` は併用できない。
- align を使うには orthogroup mode が要る。

**Circular** — `gbdraw/circular.py:254-296`

- `--conservation_blast`、`--conservation_table`、`--conservation_fasta`、`--conservation_reference {query,subject,auto}`、labels、colors、ring width/gap。
- threshold は `gbdraw/cli_utils/common.py:314-335`。

**Table の形式**

- `--comparisons_table` の列は `blast`、`query`、`subject`（`gbdraw/io/cli_tables.py:213, 566-631`）。
- `--conservation_table` の列は `blast`、`label`、`color`、`comparison_fasta`（`cli_tables.py:212, 274-334`）。

**Top-level help**

- `--protein_blastp_mode` と `--losatp_threads` を載せている（`gbdraw/cli.py:105, 116-117`）。

**Session option**

- `--session`、`--save_session`、`--session_output`（`gbdraw/cli_utils/session.py:88-112`）。
- CLI の session は raw `losatCache` を書くが、derived cache は空にする（`cli_utils/session.py:359-457`）。

### 2.4 永続形式と履歴（互換義務の根拠）

| 対象 | origin/main | release tag | dev | 互換義務 |
|---|---|---|---|---|
| Session version | 42（`main:gbdraw/session_io.py:40`） | 0.13.0 = 30（`0.13.0:gbdraw/session_io.py:25`） | 44（`session_io.py:40`） | 27–33 と 39–42 には reader が要る。44 は dev 専用なので in-place で書き直せる |
| Canonical request schema | 7（`main:gbdraw/session_request_codec.py:93`） | なし（0.13.0 には codec がない） | 8（`session_request_codec.py:109`） | 7 以下には reader が要る。8 は dev 専用 |
| `generatedProteinComparison` と settings（`losatpBin`、`ncbiBlastpBin`、`losatpThreads`、`proteinBlastpMaxHits`、`proteinBlastpCandidateLimit` など） | あり（`main:session_request_codec.py:3073, 3209`。#287、2026-07-15） | なし | `session_request_codec.py:392-410, 3304-3329` | main にあるので wire 名を変えると reader が要る |
| CLI `--protein_blastp_mode`、`--losatp_bin`、`--protein_blastp_max_hits`、`--protein_blastp_candidate_limit`、`--align_orthogroup_feature` | あり | 0.11.0 から（#192 `aca9651b`） | あり | legacy session（27–30）の `cliInvocation.args` に入っている（`0.13.0:gbdraw/session_io.py:986`）。argv の書き換えが要る |
| `--losatp_threads` | あり | 0.12.0 から | あり | 同上 |
| `--ncbi_blastp_bin` | あり | 0.12.0 から | あり | 同上 |
| `--protein_blastp_output` | あり（#319） | なし | あり | canonical session（31 以上）は `renderRequest` から replay するので、argv reader は要らない |
| `setup-losat`、`--comparisons_table`、`--conservation_table`、`linearComparisonPlan` | あり | なし | あり | 変えない |
| `SimilarityAlignmentReference` | なし（dev 専用、#580） | なし | あり | 自由に変えられる |
| hyphen 表記の alias（`--protein-blastp-mode` など） | なし | 0.13.0 にあった | なし | 既存の canonicalizer が書き換える（`session_io.py:2964-2995`） |

前例として、`docs/SESSION_COMPATIBILITY.md:236-255`「Retired inputs」がある。fresh な CLI / Python request は古い名前を拒否し、対応する古い session は replay の前に書き換える。`--collinear_max_gene_gap`（0.13.0 で release 済み）がこの扱いを受けている。

### 2.5 実測（本調査で実行。scratchpad の pinned binary を使い、repo は変更していない）

| 比較 | 結果 |
|---|---|
| lock の Linux binary（v0.1.0 の download。sha256 は lock と一致）の subcommand | `blastn`、`blastp`、`tblastx`。option は NCBI 式の単一 dash（`-query`、`-task`、`-query_gencode`、`-num_threads`、`-outfmt`） |
| `LOSAT blastn`（既定 megablast）、lambda→DE3 と Web 由来の `gbdraw/web/tutorial-data/lambda-de3-comparison/lambda-de3.losatn.tsv` | **byte 一致**（6 行） |
| `LOSAT tblastx`（既定値）と `lambda-de3.tlosatx.tsv`（397 行） | **row 集合が一致**。同点の行だけ順序が違う |
| 上の 2 つの TLOSATX 表を `-b` で描画した SVG | **byte 一致**（247,436 B）。行の順序は描画に影響しない |
| `tblastx` に gencode を渡す（danio 2/2、drosophila 5/2、C. elegans 5/2）。比較先は `metazoan-mitochondria-comparison/*.tlosatx.tsv` | **3 表とも byte 一致** |
| `-num_threads 1` と `4`、`8`（tblastx） | 出力は byte 一致 |
| NCBI BLAST+ 2.17.0 の `blastn` / `tblastx -subject`（lambda/DE3） | LOSAT と row 集合が一致。一般には保証されない（`docs/CLI_Reference.md:1391-1393`） |
| `LOSAT blastn -task dc-megablast` | **v0.1.0 は拒否する**（`possible values: megablast, blastn`）。Web は選択肢に出しており、旧 dev build は受け付ける |

`docs/internal/DOCUMENTATION_RENOVATION_PLAN_2026-08-03.md:399, 414` によると、lambda/DE3 の表は「pinned browser WASM, serial, one thread」で作られた。mtDNA の表は native LOSAT で作られた。

### 2.6 既知の不整合（この設計で解消する）

- `docs/CLI_Reference.md:1386-1391` は「The current package bundles LOSAT for Linux x86_64」と書いているが、packaging は除外している（`_build_support.py:119`）。`docs/INSTALL.md:88` と `docs/REFERENCE/command-line.md:44` が正しい。
- `docs/CLI_Reference.md:543` と `docs/REFERENCE/comparison-programs-thresholds-and-results.md:25-29` は「CLI does not run LOSATN or TLOSATX / rings」と書いている。実装後は誤りになる。
- E-value の scope が Web と CLI で違う（`comparison-programs-thresholds-and-results.md:91-100`、`run-info.js:1690-1692`）。
- LOSATP の既定値が違う。Collinear で Web は `candidateLimit=5`、`orthogroupMemberMaxHits=5`、`collinearInferOrthogroups=false`（`session-active-config-contract.js:45-55`）。CLI / API は unbounded と `True`（`options.py:1011-1016`）。違いは文書化されている（`comparison-programs-thresholds-and-results.md:83-89`）。
- `--align_orthogroup_feature` 以外にも、`--collinear_infer_orthogroups` と member hits の CLI flag がない。そのため Web の LOSATP 設定を CLI で完全には表せない。

---

## 3. 設計

### 3.1 option の model

1 つの comparison を次の要素で表す。この形を Python の typed boundary（`gbdraw/api/options.py`）に data として持たせ、CLI、入門 API、Web の request はそこへ写す。

```text
ComparisonPlan
  evidence source : table（file / DataFrame）| losat(program)
  program         : losatn(blastn) | tlosatx(tblastx) | losatp(blastp)      # 1 diagram につき 1 つ（Web と同じ）
  pairs           : Linear  = adjacent（隣接 row の全 record pair）| explicit edges（隣接 row の間に限る）
                    Circular = displayed record（subject / database）× 各 comparison FASTA（query）
  search args     : losatn_task, gencode（record ごと）, losatp の raw 上限
  display         : Linear nucleotide = pairwise matches（ribbon / curve）
                    Linear losatp     = similarity_groups | collinear | pairwise
                    Circular          = ring（既存の sequence_conservation slot）
  evidence scope  : losatp similarity_groups = 全 pair と self
                    losatp collinear         = adjacent | all（既存の --collinear_search_scope）
  thresholds      : 既存の evalue / bitscore / identity / alignment_length（mode profile の既定値。検索後に適用）
  runtime         : LOSAT（explicit / conda / managed / bundled / PATH）→ NCBI BLAST+ <program>
```

Linear で「all pairs」を ribbon として描くことはできない。edge は隣接 row の間でなければならない（`gbdraw/linear_comparison.py:123-152`、Web の `linear-comparisons.js:182-231`）。そのため nucleotide program の「all」は、multi-record row での adjacent（隣接 row の全 cross pair）か、明示 table で表す。「全 pair を evidence にする」は LOSATP の derived mode だけが持つ概念とする。

### 3.2 CLI の spelling: 案 A / B / C

**案 A（推奨）: Web の label に揃えた直交 option と、table の `source` 列**

| option | mode | 意味 / 既定値 / 制約 |
|---|---|---|
| `--losat {losatn,tlosatx,losatp}` | linear, circular | LOSAT Mode。Circular は `losatn` と `tlosatx` だけ |
| `--losatp_mode {similarity_groups,collinear,pairwise}` | linear | 既定は `similarity_groups`（Web の既定値）。`--losat losatp` のときだけ使える |
| `--losatn_task {megablast,blastn}` | linear, circular | 既定は `megablast`。`losatn` のときだけ使える |
| `--losat_gencode CODE [CODE ...]` | linear, circular | `tlosatx` のときだけ使える。Linear は全 record に 1 つ、または record ごとに 1 つ。`--records_table` の `losat_gencode` 列も使える。Circular は displayed（reference）record の gencode。既定は 1 |
| `--conservation_losat_gencode CODE [CODE ...]` | circular | ring ごとの gencode。`--conservation_fasta` と同じ順。`--conservation_table` の `losat_gencode` 列も使える |
| `--losat_bin PATH` / `--ncbi_blast_bin PATH` | 両方 | 両方を同時には指定できない。NCBI 側は選んだ program の executable（`blastn` / `tblastx` / `blastp`） |
| `--losat_threads N` | 両方 | 各 job の `-num_threads`。job は順に実行する |
| `--losatp_max_hits N` | linear | pairwise 表示の上限。既定 5 |
| `--losatp_max_target_seqs N\|none` | linear | raw の `-max_target_seqs`。既定 none |
| `--losatp_member_max_hits N\|none` | linear | 新設（Web の Member hits per protein）。既定 none |
| `--similarity_alignment_feature ID` | linear | 旧 `--align_orthogroup_feature`。`similarity_groups` のときだけ使える |
| `--losat_output_dir DIR` | 両方 | raw evidence を出力する（§3.7） |
| `--comparisons_table` の `source` 列 | linear | `table`（既定）または `losat`。`losat` の行では `blast` を空にする |

**例（案 A）**

```bash
# 1. Linear LOSATN、adjacent ribbon
gbdraw linear --gbk NC_001416.gb NC_042057.1.gb --losat losatn -o lambda-de3 -f svg

# 2. Linear TLOSATX「all pairs」。2 row × 2 record で a→c, a→d, b→c, b→d を検索する
gbdraw linear --gbk a.gb b.gb c.gb d.gb \
  --multi_record_position '#1@1' --multi_record_position '#2@1' \
  --multi_record_position '#3@2' --multi_record_position '#4@2' \
  --losat tlosatx --losat_gencode 11 -o abcd
#    明示 pair（LOSAT と precomputed の混在も可）:
#    pairs.tsv:  source  blast         query  subject
#                losat                 #1     #2
#                table   b_c.tsv       #2     #3
gbdraw linear --gbk a.gb b.gb c.gb --comparisons_table pairs.tsv --losat tlosatx -o abc

# 3. Linear LOSATP、Similarity groups を protein ID で整列
gbdraw linear --records_table bgc_records.tsv --losat losatp \
  --losatp_mode similarity_groups --similarity_alignment_feature CAG38695.1 \
  --show_labels first --pairwise_match_style curve -o BGC0000708-BGC0000713 -f interactive_svg

# 4. Linear LOSATP、Collinear blocks
gbdraw linear --gbk g1.gb g2.gb g3.gb --losat losatp --losatp_mode collinear \
  --collinear_search_scope all --collinear_min_anchors 2 --losat_threads 8 -o collinear

# 5. Circular similarity ring。LOSATN で 2 genome と比較する
#    comparison FASTA が query、displayed record が subject（database）
gbdraw circular --gbk ref.gb --losat losatn \
  --conservation_fasta s1.fna s2.fna --conservation_labels S1 S2 -o ref-rings
#    TLOSATX（T-CLI-09 と同じ入力）:
gbdraw circular --gbk HmmtDNA.gbk --losat tlosatx --losat_gencode 2 \
  --conservation_fasta NC_002333.2.fna NC_024511.2.fna NC_001328.1.fna \
  --conservation_losat_gencode 2 5 5 --identity 40 -o mt-rings

# 6. precomputed（現状のまま）
gbdraw linear --gbk a.gb b.gb -b a_b.blast.tsv
gbdraw linear --gbk a.gb b.gb c.gb --comparisons_table comparisons.tsv
gbdraw circular --gbk ref.gb --conservation_blast s1.tsv s2.tsv --conservation_fasta s1.fna s2.fna
```

**案 B: 複合 spec `--compare PROGRAM[:DISPLAY]`**

```bash
gbdraw linear --gbk a.gb b.gb --compare losatn
gbdraw linear --gbk a.gb b.gb c.gb d.gb --multi_record_position ... --compare tlosatx --compare_gencode 11
gbdraw linear --records_table bgc.tsv --compare losatp:similarity_groups --compare_align CAG38695.1
gbdraw linear --gbk g1.gb g2.gb g3.gb --compare losatp:collinear --collinear_search_scope all
gbdraw circular --gbk ref.gb --compare losatn --compare_with s1.fna s2.fna
gbdraw linear --gbk a.gb b.gb -b a_b.tsv            # precomputed は現状のまま
```

- 利点: 指定が 1 つで済む。
- 欠点:
  - 小さな言語の parser が要り、argparse の choices や help、生成した CLI_Reference に出ない。
  - `:` 以降が LOSATP 専用なので、検証は結局 2 段になる。
  - Web の 2 つの control（LOSAT Mode と LOSATP mode）に 1 対 1 で対応しない。

**案 C: table が中心。table の各行に program を書ける**

```bash
gbdraw linear --gbk a.gb b.gb --comparisons adjacent:losatn
gbdraw linear --gbk ... --comparisons_table pairs.tsv     # source 列 = losatn|tlosatx|losatp|<path>
gbdraw linear --records_table bgc.tsv --comparisons adjacent:losatp --losatp_mode similarity_groups
gbdraw circular --gbk ref.gb --conservation_table rings.tsv   # source 列 = losatn|tlosatx|<path>、comparison_fasta 列
```

- 利点: 混在する edge を自然に書ける。
- 欠点:
  - 行ごとに program が変わる。Web は 1 Generate につき 1 program なので、これと矛盾する。
  - LOSATP の derived mode は edge ごとの概念ではない。
  - 最も簡単な使い方にも table か mini-syntax が要る。

**推奨: 案 A（案 C から `source` 列だけを取り入れる）**

- Web の 2 つの control と state（`losatProgram`、`losat.blastp.mode`）に 1 対 1 で対応する。
- argparse の choices でそのまま検証でき、生成した help に全部の値が出る。
- 混在する edge は既存の `--comparisons_table` に列を 1 つ足すだけで済む（Web の selected plan と同じ意味）。
- 新しい概念は増えない。

### 3.3 Circular similarity ring と LOSAT

**Query と reference**

- Web と同じにする（`run-analysis.js:2643-2669, 2797`）。
- query は comparison FASTA。FASTA の全 record を使う。
- subject（database）は displayed circular の入力 source の全 record。
- reference 側は `subject` に固定する。
- これで E-value の database がどの ring でも reference genome になり、ring の間で比較できる。
- 既存の mtDNA の表も同じ向きで作られている（`tools/build_metazoan_mitochondria_comparison_fixture.py:186-208`）。Owner の例にある「subject genome」は、gbdraw の用語では comparison genome（LOSAT の query）にあたる。docs ではこの用語を使う。

**Program**

- `losatn` と `tlosatx` だけ。`losatp` は `COMPARISON_INPUT` / `RING_LOSAT_PROGRAM` で拒否する。Web も対応していない。

**Gencode**

- reference は `--losat_gencode`（Web の `subject_gencode` → `--db-gencode`）。
- ring は `--conservation_losat_gencode`（Web の ring `losat_gencode` → `--query-gencode`）。
- 既定はどちらも 1。

**入力**

- `--losat` を指定したときは、既存の `--conservation_fasta` が ring の query 一覧になる（1 本以上必須）。
- `--conservation_blast` は併用できない。
- `--conservation_reference` は `auto` か `subject` だけ受け付け、解決後の値は `subject`。
- labels と colors の数は `--conservation_fasta` に揃える。
- `--conservation_table` は、`--losat` があるときは `blast` 列を禁止し、`comparison_fasta` 列を必須にする。`losat_gencode` 列は任意。
- comparison FASTA は interactive SVG の span export にもそのまま使える。

**Ring の順序と描画**

- ring は `--conservation_fasta` の順に並ぶ。`--conservation_blast` と同じ規則（`track_index` 1..n、`gbdraw/api/diagram.py:519-568`）。
- planner は LOSAT の結果を `conservation_dataframes` 相当に解決する。threshold（Circular の既定 1e-5/50/70/0）、`sequence_conservation` slot、`--circular_track_slot` との結合（`source_index`）、renderer は変えない。
- multi-record canvas では、既存どおり最初の record に描く（`diagram.py:2834-2870`）。行は `sseqid` と reference record の照合で絞る（`analysis/conservation.py:346-347`）。

### 3.4 Python API

**Typed API**（`gbdraw.api`、owner は `gbdraw/api/options.py`）

```python
LosatProgram = Literal["losatn", "tlosatx", "losatp"]
LosatpMode = Literal["similarity_groups", "collinear", "pairwise"]

@dataclass(frozen=True)
class LosatRuntimeOptions:
    losat_executable: str | None = None        # None = 自動解決
    ncbi_blast_executable: str | None = None   # 選んだ program の NCBI executable
    threads: int | None = None

@dataclass(frozen=True)
class LosatSearchOptions:
    program: LosatProgram
    losatn_task: Literal["megablast", "blastn"] = "megablast"
    record_gencodes: tuple[int, ...] = ()       # 空なら全 record 1。Circular では (reference,)
    pairs: tuple[tuple[int, int], ...] | None = None   # None = adjacent rows（Linear のみ）
    losatp_mode: LosatpMode | None = None       # program == "losatp" のとき必須
    losatp_max_hits: int = 5
    losatp_max_target_seqs: int | None = None
    losatp_member_max_hits: int | None = None
    runtime: LosatRuntimeOptions = LosatRuntimeOptions()

LinearDiagramOptions.losat_search: LosatSearchOptions | None = None
CircularDiagramOptions.losat_search: LosatSearchOptions | None = None
CircularDiagramOptions.conservation_losat_gencodes: Sequence[int] | None = None
```

- 削除する typed field:
  - `protein_blastp_mode`、`protein_comparison_pairs`
  - `losatp_bin`、`ncbi_blastp_bin`、`losatp_threads`
  - `protein_blastp_max_hits`、`protein_blastp_candidate_limit`、`orthogroup_member_max_hits`
- 残す field: `collinearity_*`、`collinear_infer_orthogroups`、`orthogroup_membership_mode`、precomputed 系（`blast_files`、`linear_comparisons`、`protein_comparisons`、`orthogroups`、`collinearity_blocks`、`comparison_table_file`）。いずれも LOSATP の derived 設定か precomputed input で、名前の問題がない。
- 検証は 1 か所にまとめる（R7、`gbdraw/web/CLAUDE.md:286-300`）。
  - program と option の整合は `LosatSearchOptions.__post_init__` で見る。
  - mode ごとの制約は `LinearDiagramOptions` / `CircularDiagramOptions` の検証で見る。Circular は `losatp`、`pairs`、`losatp_*` を拒否し、Linear は precomputed 入力との衝突を見る。
  - CLI はこれらを写すだけで、検証は重複させない。
- `SimilarityAlignmentReference(feature_id)` は `losat_search.program == "losatp"` かつ `losatp_mode == "similarity_groups"` を要求する。今は `protein_blastp_mode == "orthogroup"` を要求している（`requests.py:748-756`）。dev 専用なので自由に変えられる。
- planner は nucleotide の intent を `LinearComparison` と cache entry に解決してから render と encode を行う。`SimilarityAlignmentReference` と同じく、未解決の intent は request に encode しない（`session_request_codec.py:619-630` と同じ扱い）。

**入門 API**（`gbdraw`）

```python
LinearComparisonOptions(
    losat: LosatProgram | None = None,
    losatp_mode: LosatpMode = "similarity_groups",
    losatn_task: str = "megablast",
    gencodes: int | Sequence[int] = 1,
    pairs: Sequence[tuple[int, int]] | None = None,      # 全 program の明示 edge に拡張
    max_hits: int = 5, max_target_seqs: int | None = None, member_max_hits: int | None = None,
    losat_executable: str | None = None, ncbi_blast_executable: str | None = None, threads: int | None = None,
    # 変えない field: blast_files, comparisons, protein_comparisons, orthogroups, match_style,
    #   collinearity_*, orthogroup_membership, max_paralog_links, similarity_alignment
)
ComparisonRingOptions(tracks=..., reference="auto", ring_width=None, ring_gap=None,
    losat: Literal["losatn", "tlosatx"] | None = None, losatn_task="megablast",
    reference_gencode: int = 1, losat_executable=None, ncbi_blast_executable=None, threads=None)
ComparisonRingTrackOptions(source: TableSource | None = None, label=None, color=None,
    comparison_sequence_source=None, losat_gencode: int = 1)   # losat があれば source=None、comparison_sequence_source は必須
```

例:

```python
draw_linear(records, options=LinearOptions(comparisons=LinearComparisonOptions(losat="losatn")))
draw_linear(records, options=LinearOptions(comparisons=LinearComparisonOptions(
    losat="losatp", losatp_mode="similarity_groups",
    similarity_alignment=SimilarityAlignmentReference("CAG38695.1"))))
draw_circular(ref, options=CircularOptions(comparison_rings=ComparisonRingOptions(losat="tlosatx",
    reference_gencode=2, tracks=[ComparisonRingTrackOptions(comparison_sequence_source="NC_002333.2.fna",
    losat_gencode=2, label="Danio")])))
```

### 3.5 名前の変更と互換

**CLI**

- fresh run では古い flag を拒否し、置き換え先を示す。
- legacy session（27–30）の `cliInvocation.args` は、既存の `_canonicalize_legacy_session_cli_args`（`session_io.py:2964-2995`）で書き換える。
- 古い flag と新しい flag の対応表は、拒否と書き換えの両方が参照する 1 つの table にする。
- この table は 1 token を複数 token に書き換えられるように拡張する。

| 旧 | 新 | 0.13.0 の session reader |
|---|---|---|
| `--protein_blastp_mode none` | （削除。comparison なし） | 要る |
| `--protein_blastp_mode pairwise` | `--losat losatp --losatp_mode pairwise` | 要る |
| `--protein_blastp_mode orthogroup` | `--losat losatp --losatp_mode similarity_groups` | 要る |
| `--protein_blastp_mode collinear` | `--losat losatp --losatp_mode collinear` | 要る |
| `--losatp_bin X` / `--losatp-bin X` | `--losat_bin X` | 要る |
| `--ncbi_blastp_bin X` | `--ncbi_blast_bin X` | 要る |
| `--losatp_threads N` | `--losat_threads N` | 要る |
| `--protein_blastp_max_hits N` | `--losatp_max_hits N` | 要る |
| `--protein_blastp_candidate_limit N` | `--losatp_max_target_seqs N` | 要る |
| `--align_orthogroup_feature ID` | `--similarity_alignment_feature ID` | 要る。0.13.0 の意味（group に含まれる feature）は exact ID と同じ解決になる（現行の legacy 処理、`session_request_codec.py:3569-3580`） |
| `--protein_blastp_output F` | `--losat_output_dir DIR` | 不要。release に入っておらず、canonical session は `renderRequest` から replay する |

**Python**

- typed field と入門 API field は改名し、古い kwarg は Python 標準の `TypeError` で落とす。alias を受け付ける shim は作らない。
- 対応は `docs/SESSION_COMPATIBILITY.md` の「Retired inputs」と release notes に載せる。
- 入門 API の対応:
  - `protein_mode` → `losat` + `losatp_mode`
  - `blastp_executable` → `ncbi_blast_executable`
  - `candidate_limit` → `max_target_seqs`
  - `orthogroup_member_max_hits` → `member_max_hits`
  - `losat_executable` の既定値は `"losat"` から `None` に変える。

**Persisted wire 名は変えない**

- `renderRequest.comparisons[].kind = "generatedProteinComparison"`、`mode ∈ {none,pairwise,orthogroup,collinear}`。
- settings key（`losatpBin`、`ncbiBlastpBin`、`losatpThreads`、`proteinBlastpMaxHits`、`proteinBlastpCandidateLimit`、`orthogroupMemberMaxHits` ほか）。
- Web state（`losatProgram`、`losat.blastp.mode`、`circularConservation.losat_program`、`losat_gencode`）と cache key の args 表記（v1）。
- codec は今 `_camel(field_name)` で key を作っている（`session_request_codec.py:3313-3319`）。これを明示的な wire 名 table に置き換える。typed 名と wire 名を切り離すだけで、reader も migrator も増えない。

**dev 専用の形式**

- session 44、request 8、`SimilarityAlignmentReference`、新しい `cliInvocation.args` がこれにあたる。
- 変更は in-place で書き直し、中間形の reader は作らない（root `CLAUDE.md:138-148`）。

### 3.6 Runtime

**Owner**

- `gbdraw/comparisons/losat_runtime.py` を新設する。この path は CI 上で `losat-integration` に分類される（`tools/ci-impact-policy.mjs:167-169`）。
- ここへ移すもの:
  - runtime 解決と platform 判定（`protein_colinearity.py:3146-3341`）
  - argv 構築と subprocess 実行（`protein_colinearity.py:5885-6042`）
  - raw cache manager（`protein_colinearity.py:2230-2581`）を一般化したもの
- 移した path は元の場所から削除する。
- protein の identity と analysis は `protein_colinearity.py` に残す。
- install と cache の owner は `losat_setup.py` のまま変えない。

**Program ごとの対応**

| | LOSAT argv（v0.1.0） | cache key の args（Web v1） | NCBI fallback |
|---|---|---|---|
| losatn | `blastn -query q -subject s -outfmt 6 -task T [-num_threads N]` | `["--task", T]` | `blastn -task T ...` |
| tlosatx | `tblastx ... -query_gencode Q -db_gencode D` | `["--query-gencode", Q, "--db-gencode", D]` | `tblastx -query_gencode Q -db_gencode D ...` |
| losatp | 現行のまま（`-max_hsps`、`-max_target_seqs`） | 現行のまま | `blastp ...` |

**解決順序**（全 program 共通）

1. `--losat_bin`
2. `--ncbi_blast_bin`
3. conda の `losat`
4. managed（`setup-losat`）
5. bundled（source checkout だけ）
6. PATH の `losat`
7. PATH の NCBI `<program>`
8. どれもなければ `LOSAT_RUNTIME` / `UNAVAILABLE`

managed と bundled の LOSAT 0.1.0 は、3 つの program をすべて持っている（§2.5）。program による platform の違いはない。

**Bundled binary**

- track している dev build（`3d357ae9…`）を、lock の Linux binary（`fb01e0f6…`）に置き換える（判断 D9）。
- これで CI と tutorial が利用者と同じ binary を使うようになる。
- `tools/build_metazoan_mitochondria_comparison_fixture.py:22-23, 193-208` の sha と v1 flag も直す。新しい binary でも 3 表が byte 一致することは確認済み。

**Job plan**

- Python の `plan_losat_jobs` は Web の `planLosatSourceJobs` と同じ規則に従う（source file = genome、file 内は query record を除く、self は要求時だけ、`searchContext`）。
- `tests/fixtures/losat_job_plan_cases.json` を Python と `tests/web/linear-sources.test.mjs` の両方で実行する。
- FASTA の抽出（canonical hash の入力）も `tests/fixtures/losat_fasta_extraction_cases.json` で固定する。
- これで CLI と Web の raw cache key が一致する。CLI で保存した session を Web で開くと、検索せずに再利用できる。逆も同じ。
- LOSATP にも同じ scope を使う（判断 D7）。そうすれば `LOSAT_DATABASE_SCOPE_NOTE` の CLI についての注記を削除できる。

**Thread と決定性**

- CLI は job を順に実行し、各 job に `--losat_threads` を渡す。thread 数を変えても出力が同じことは確認済みで、test でも固定する。
- raw TSV は LOSAT の出力順のまま残す。
- 同点の行の順序は runtime によって違ってよい。同等性は row の多重集合で判定する。描画の結果が同じことは §2.5 で確認済み。

**Cache**

- CLI の 1 回の run の中では、同じ raw key の検索を 1 回だけにする（CW-02）。collinear の順方向と逆方向の job などがこれにあたる。
- run をまたぐ再利用は既存の 2 つで足りる。session replay の `losatCache` と、`--losat_output_dir` の出力を `--comparisons_table` / `--conservation_table` に渡す方法。
- 長く残る on-disk cache は作らない（CW-06）。
- derived payload cache は Web の対話専用のままにする。

**Runtime の記録**

- `losatCache` entry に、key に入らない field として `runtime: {kind: "losat"|"ncbi-blast", version, source}` を加える。session 44 は dev 専用なので in-place で足せる。
- Web は WASM の build identity を同じ field に書く。これで Run Info と Session が「program と version」を記録できる（`comparison-programs-thresholds-and-results.md:62-65`）。

**WASM**

- Web の WASM を LOSAT v0.1.0（candidate `6bfb1b09…`）から build し直して pin する（判断 D11）。
- それまでは、`dc-megablast` は Web だけの選択肢として残る。CLI の choices は `megablast` と `blastn`。

### 3.7 Error、Session、Source recipe / Run Info

**Diagnostic**

- producer が code を持つ（R6、`gbdraw/web/CLAUDE.md:267-285`）。
- wording は `gbdraw/web/js/services/error-normalization.js` に足す。adapter の whitelist（`error_adapter.py:45-53`）と、`tests/test_web_error_producer_coverage.py` の部分集合 test も同時に更新する。

| code / reason | いつ出るか | context |
|---|---|---|
| `COMPARISON_INPUT` / `LOSAT_OPTION_PROGRAM` | ほかの program 用の option を指定した（例: `--losatp_mode` と `--losat losatn`、`--losat_gencode` と losatn、`--collinear_*` を collinear 以外で指定） | `field`（option 名）、`program` |
| `COMPARISON_INPUT` / `LOSAT_PLAN` | `--losat` と `-b` を併用した。table に `source=losat` の行があるのに `--losat` がない。逆に `--losat` があるのに使われる edge がない。selected の edge と LOSATP の非 pairwise を組み合わせた（Web の `selected-losat-requires-pairwise`）。record が 1 つしかない | `row`（table 行） |
| `COMPARISON_INPUT` / `RING_LOSAT_PROGRAM`、`RING_LOSAT_INPUT` | Circular で losatp を指定した。`--conservation_blast` と `--losat` を併用した。`--losat` があるのに reference が `query` | `field` |
| `LOSAT_RUNTIME` / `UNAVAILABLE`、`FAILED`、`OUTPUT` | runtime がない。終了コードが 0 以外。出力の列数や ID が不正 | `program`、`exitCode`、`row` |

- CLI の syntax error（choices 違反、retired flag）は argparse で exit 2 にする。retired flag は §3.5 の table から置き換え先を出す。
- それ以外は `GbdrawError` として typed 層から出す。CLI は message を表示し、0 以外で終了する。

**Session**（`--save_session` / `--session_output`）

- LOSATN/TLOSATX:
  - edge ごとに search frame の raw TSV を `nucleotideBlast` resource として保存する。
  - `losatCache` には schema 2 の entry を保存し、key は Web と同じにする。
  - Circular は `conservationBlastFiles`、`conservationFastaFiles`、`flow: 'circular-conservation'` の entry を保存する。
  - replay に LOSAT は要らない。Web が書く形と同じ。
- LOSATP: 現行どおり `generatedProteinComparison`（intent と settings）と `losatCache` を保存する。
- `cliInvocation.args` には新しい flag を書く。

**`--losat_output_dir DIR`**

- Linear nucleotide: edge ごとの raw TSV と `comparisons.tsv`（`blast`、`query`、`subject`、`#N`）。ファイル名は Web の `*.losatn.tsv` / `*.tlosatx.tsv`（`run-analysis.js:1493-1497`）に揃える。`comparisons.tsv` はそのまま `--comparisons_table` に渡せる。
- Circular: ring ごとの TSV と `conservation.tsv`（`blast`、`comparison_fasta`、`label`、`color`）。そのまま `--conservation_table` に渡せる。
- LOSATP: `losatp.raw.tsv`。中身は現行の `--protein_blastp_output` と同じ形式。

**Source recipe / Run Info（Web）**

- 当面は今の table 形式（lossless）を保つ。
- WASM を v0.1.0 に pin して shared vector が一致した後に、`--losat ...` 形式の recipe を出すようにする（判断 D14）。このとき LOSATP も recipe を出せるようになる。
- scope の parity（判断 D7）を入れた PR で、`LOSAT_DATABASE_SCOPE_NOTE` から CLI についての文を削除する。

### 3.8 Architecture の evidence（ratchet）

**PR ごとの owner と path**

- **Runtime の owner**: `protein_colinearity.py` の runtime 部分を `gbdraw/comparisons/losat_runtime.py` へ移す。移動は 1 対 1 で、元の場所から削除する。delta(OE)=0。
- **Native の nucleotide 検索**: 既存の Python runtime owner と job plan owner に program を data として足す。新しい path は作らない。
- **Job plan**: JS（`planLosatSourceJobs`）と Python の 2 か所にある。Python 側には今も protein の pair 走査があるので、この状態は既存のもの。Python 側は 1 module にまとめ、protein の pair 走査（`collinearity.py:236-251` の search pair 生成と `protein_colinearity.py:6561-6760` の走査）を置き換える。
  - capability を「LOSAT job plan」1 つとみなせば delta(OE)=0。
  - nucleotide の plan を別の capability とみなすと +1 になる。その場合は exception packet と maintainer の判断が要る（判断 D15）。
  - Web は Settings の見積もりを Python なしで出す必要がある（CW-01、`ARCHITECTURE_FITNESS_FUNCTION_RATCHET.md:400`）。そのため JS 側を Python に統合できない。
- **互換**: 既存の legacy argv canonicalizer の table に行を足すだけ。wire 名は変えないので、新しい persisted compatibility path は 0。delta(CB)=0。

**Computation**

- FASTA の抽出は record と purpose ごとに 1 回（CW-02）。
- raw の検索は key ごとに 1 回。
- 長く残る cache は足さない（CW-06）。

---

## 4. 実装計画（PR の順序）

CI class は `tools/ci-impact-policy.mjs:145-197` の分類による。`gbdraw/analysis/` と `gbdraw/bin/` は FULL_BY_DEFAULT。

### PR-1 LOSAT runtime owner を program 共通にする（利用者から見える変化なし）

- **範囲**
  - `gbdraw/comparisons/losat_runtime.py` を新設し、§3.6 の内容を移す。
  - argv と cache-args の対応 table を作る。
  - NCBI fallback を program ごとにする。
  - `losatCache` entry に `runtime` を記録する。
  - bundled binary を v0.1.0 の binary に置き換える（D9）。
  - fixture tool を直す。
- **Owner file**: `gbdraw/comparisons/losat_runtime.py`（新）、`gbdraw/analysis/protein_colinearity.py`（削除側）、`gbdraw/bin/linux-x86_64/losat`、`tools/build_metazoan_mitochondria_comparison_fixture.py`。
- **先に失敗させる test**（`tests/test_losat_runtime.py`）:
  - 3 program × 解決順序
  - LOSAT v2 と NCBI の argv golden
  - cache args の v1 golden
  - pinned binary が lambda-de3 の losatn を byte、tlosatx を multiset、mtDNA の 3 表を byte で再現すること
  - thread 数で結果が変わらないこと
- **Docs**: `docs/CLI_Reference.md:1386-1391` の bundled についての記述を直す。`docs/REFERENCE/command-line.md:19-56` の解決順序。
- **CI**: full と losat-integration。

### PR-2 LOSATP の改名（CLI、入門 API、typed API）

- **方針**: 振る舞いは変えない。古い flag は Retired inputs にする。
- **範囲**
  - §3.2 の LOSATP 系 flag。`--losat` の choices はまだ `losatp` だけ。
  - `--losatp_member_max_hits` と `--collinear_infer_orthogroups {on,off}` を新設する（D8）。
  - `LosatSearchOptions` / `LosatRuntimeOptions` を作る。
  - 入門 API を改名する。
  - codec に wire 名 table を作る。
  - legacy argv の書き換え table を作る（拒否にも使う）。
  - error field map を更新する（`error_adapter.py:33,42,218`、`error-normalization.js:12,21,145,203,213,220`）。
  - `cli.py:105,116-117` を更新する。
- **先に失敗させる test**
  - 各旧 flag が exit 2 になり、置き換え先を示すこと。
  - `0.13.0:gbdraw/web/gallery/sessions/BGC0000708-BGC0000713.gbdraw-session.json`（release tag の positive fixture）を replay した SVG が変わらないこと。
  - schema 7 の fixture（`tests/fixtures/sessions/*v42*`）を decode した結果が新しい typed field に入り、re-encode した wire key が同じであること。
  - 入門 API の新しい field。
- **Docs / Gallery**
  - CLI_Reference の help を再生成する（`tools/update_cli_reference_help.py`、`tests/test_cli_reference.py:37`）。
  - T-CLI-08、T-CLI-10、T-PY-05、`docs/RECIPES.md:322-340`、`docs/REFERENCE/{command-line,python-api,comparison-programs-thresholds-and-results}.md`、`docs/SESSION_COMPATIBILITY.md:236-255`。
  - `docs/recipes/run_cli_scenarios.py:1170-1198`、`docs/scenarios/manifest.json:1086-1088`、`docs/internal/SCENARIO_EVIDENCE.md`（H-CLI-06/07/08）。
  - Gallery の `examples.json` にある 5 つの command と 5 つの session を再生成する（`tools/prepare_interactive_gallery_assets.py`、`tools/refresh_gallery_sessions.py`）。SVG は byte 一致すること。
  - `tools/verify_losat_installation.py`、`tools/reproduce_examples_manifest.py`、`tools/prepare_issue597_s01_fixture.py`。
- **CI**: python-core、session-persistence、gallery、documentation、web-runtime（playwright-functional を含む）。

### PR-3 Linear LOSATN / TLOSATX

- **範囲**
  - `--losat losatn|tlosatx`、`--losatn_task`、`--losat_gencode`、records_table の `losat_gencode` 列、comparisons_table の `source` 列、`--losat_output_dir`。
  - Python の `plan_losat_jobs`（source batching）と planner での解決。
  - session に resource と cache を書く。
  - 新しい diagnostic を Web の wording と合わせて足す。
- **先に失敗させる test**
  - `gbdraw linear ... --losat losatn`（T-CLI-07 の threshold を指定）が `docs/images/t-cli-07/lambda-de3-losatn.svg` と byte 一致すること。
  - tlosatx で record ごとの gencode が効くこと。
  - LOSAT と table の edge を混在できること。
  - `LOSAT_OPTION_PROGRAM` / `LOSAT_PLAN` が出ること。
  - LOSAT がない環境で保存した session を replay すると同じ SVG になること。
  - shared vector（job plan、FASTA、cache key）を Python と `tests/web/linear-sources.test.mjs` の両方で実行すること。
- **Docs**
  - T-CLI-07 と Python 版を直接実行の手順にする。table の手順は同じ page に残す（新しい page は作らない）。
  - `comparison-programs-thresholds-and-results.md:7-29` の capability matrix。
  - `input-formats-and-tsv-schemas.md`、FAQ の該当項目、CLI_Reference の prose。
- **CI**: losat-integration、python-core、web-runtime、documentation。

### PR-4 Circular ring を LOSAT で作る

- **範囲**
  - circular の `--losat losatn|tlosatx`、`--conservation_losat_gencode`、conservation_table の `losat_gencode` 列。
  - `ComparisonRingOptions.losat` と typed の `losat_search`。
  - reference を `subject` に固定する。
- **先に失敗させる test**
  - §3.2 例 5 の TLOSATX command（T-CLI-09 と同じ label、color、threshold を指定）が `docs/images/t-cli-09/precomputed_circular_rings.svg` と byte 一致すること。
  - losatp を拒否すること。
  - `--conservation_blast` との併用と reference `query` を拒否すること。
  - session を replay できること。
- **Docs**
  - T-CLI-09 と Python 版に直接実行の手順を足す。
  - `docs/CLI_Reference.md:543` の記述を直す。
  - FAQ。
- **CI**: losat-integration、python-core、web-runtime、documentation。

### PR-5 LOSATP の E-value scope を Web に揃える（D7 を採るとき）

- **範囲**: LOSATP の job plan を source batching に切り替える。
- **削除するもの**: `run-info.js:1690-1692` と `comparison-programs-thresholds-and-results.md:97-99` の CLI についての文。
- **先に失敗させる test**: 複数 record の file で、Web の vector と同じ job と cache key になること。
- **Gallery**: 出力が変わる例だけを再生成し、差分を review する。
- **CI**: losat-integration、gallery、web-runtime。

### PR-6 Web 側の後続（D11、D14 を採るとき）

- **範囲**
  - LOSAT v0.1.0 から WASM を build して pin する。identity を記録する。
  - `dc-megablast` を扱う（v0.1.0 が対応しなければ選択肢から外す）。
  - ring 行の "Subject gencode" を "Comparison gencode" に改める。
  - Source recipe で `--losat` 形式を出す。
- **CI**: losat-integration、web-runtime、gallery。

**順序の理由**

- PR-1 と PR-2 は振る舞いを変えない。SVG と Gallery が byte 一致することで正しさを確かめられる。
- PR-3 と PR-4 は新しい capability。既存の reference SVG を acceptance に使える。
- PR-5 と PR-6 は科学的な出力や Web を変える。PR-1〜4 と切り離して review する。

---

## 5. Owner の判断が要る点（推奨付き）

| # | 判断 | 推奨 |
|---|---|---|
| D1 | CLI の spelling | **案 A**（`--losat`、`--losatp_mode`、table の `source` 列） |
| D2 | program の token | **`losatn`/`tlosatx`/`losatp`**（Web の label）。`blastn` 系の token は NCBI BLAST+ と紛らわしい |
| D3 | LOSATP 表示の token と既定値 | **`similarity_groups`/`collinear`/`pairwise`。既定は `similarity_groups`**（Web と同じ）。wire の `orthogroup` は変えない |
| D4 | 旧 CLI flag の扱い | **fresh run では拒否して置き換え先を示す。legacy session（27–30）は書き換える**（Retired inputs の前例）。代わりに 1 release だけ hidden alias と warning を残す案もある |
| D5 | 旧 Python field の扱い | **alias なしで改名**（標準の `TypeError`、Retired inputs 表に記載） |
| D6 | Persisted wire 名 | **変えない**（`generatedProteinComparison`、settings key、`mode:"orthogroup"`、Web state）。変えると schema 7 の reader が増え、CB +1 の exception になる |
| D7 | E-value の database scope | **CLI / API も Web と同じ source-file batching にする（LOSATP を含む）**。複数 record の file では CLI の LOSATP の E-value が変わる |
| D8 | LOSATP の既定値の差（Collinear の上限と inference） | **この campaign では既定値を変えない**。`--losatp_member_max_hits` と `--collinear_infer_orthogroups` を新設し、Web の設定を CLI で表せるようにする |
| D9 | bundled binary と CLI_Reference の誤記 | **track している dev build を v0.1.0 の binary に置き換え、記述を直す** |
| D10 | NCBI BLAST+ fallback | **全 program で自動の PATH fallback にする**（blastp の現行と同じ）。runtime は Session に記録する |
| D11 | Web WASM の pin | **LOSAT v0.1.0 から build して pin する**（別 PR）。`dc-megablast` は v0.1.0 が対応しないので外す |
| D12 | Circular の範囲 | **LOSATN/TLOSATX だけ。入力は FASTA だけ。reference は subject に固定**（Web と同じ） |
| D13 | raw の出力 | **`--losat_output_dir`**（edge ごとの TSV と、そのまま再利用できる manifest）。`--protein_blastp_output` はこれに置き換える |
| D14 | Web の Source recipe | **D11 が済むまでは table 形式。済んだら `--losat` 形式** |
| D15 | Job plan の owner | **Python と JS の 2 実装を shared vector で固定する**（既存の OE を増やさない解釈）。Python に統合すると CW-01 に反する |
| D16 | 周辺の名前 | **`--align_orthogroup_feature` → `--similarity_alignment_feature`** は今回変える。`--conservation_*` と `--show_labels orthogroup_top` は今回変えない |
| D17 | TLOSATX gencode の既定値 | **1（Web と同じ）**。`/transl_table` からの推定はしない |

---

## 6. 影響範囲

件数は `origin/dev` `be1bbd46` の tracked file に対する `git grep -o -F` で数えた。o は出現数、f は file 数。

**Python package（`gbdraw/`、web を除く）**

| 名前 | 件数 |
|---|---|
| `protein_blastp_mode` | 75 o / 12 f |
| `losatp_bin` | 54 o / 9 f |
| `losatp_threads` | 66 o / 10 f |
| `ncbi_blastp_bin` | 53 o / 8 f |
| `protein_blastp_max_hits` | 24 o |
| `protein_blastp_candidate_limit` | 25 o |
| `align_orthogroup_feature` | 15 o / 4 f |

主な file:

- `gbdraw/linear.py`、`gbdraw/cli.py`、`gbdraw/api/options.py`、`gbdraw/api/diagram.py:1573-1596`、`gbdraw/api/requests.py`、`gbdraw/api/request_render.py`、`gbdraw/interface.py`
- `gbdraw/session_request_codec.py`、`gbdraw/session_io.py`（`3911-3948` は Web state から CLI argv を作る）
- `gbdraw/web_support/error_adapter.py`

**Web**

- 改名に関わるのは error field map だけ（`error-normalization.js`）。wire と state は変えない。
- PR-3 以降で wording を足し、shared vector を test する。
- PR-6 で WASM と Source recipe を変える。

**Docs**

- `docs/CLI_Reference.md`: 約 55 か所（生成した help と prose）。
- `docs/REFERENCE/comparison-programs-thresholds-and-results.md`: capability matrix（7-16）、25-29、69、83-89、97-99、175。
- `docs/REFERENCE/command-line.md`: 11、19-56、141-175。
- `docs/REFERENCE/python-api.md`: 97-121、203-215。
- `docs/SESSION_COMPATIBILITY.md`: 91、236-255。
- `docs/RECIPES.md`: 322-340。`docs/INSTALL.md`: 63、88-95。
- FAQ の 38-55、143-166。
- `docs/DOCS.md` には flag 名がないので変えない。page は追加しない。既存の owner page を直す。

**Tutorial**

| ID | file | 対応する PR |
|---|---|---|
| T-CLI-07 | `docs/TUTORIALS/CLI/compare-genomes-losatn.md` | PR-3 で直接実行の手順を足す |
| T-CLI-08 | `docs/TUTORIALS/CLI/compare-proteins-losatp.md:83-156` | PR-2 |
| T-CLI-09 | `docs/TUTORIALS/CLI/add-precomputed-circular-comparison-rings.md:86-112` | PR-4 |
| T-CLI-10 | `docs/TUTORIALS/CLI/compare-proteins-losatp-collinear.md:77-107` | PR-2 |
| T-PY-05 | `docs/TUTORIALS/PYTHON/compare-proteins-losatp.md:51-156`（typed の `protein_blastp_mode`、`losatp_threads` を使う） | PR-2 |
| T-PY-07 | 入門 API の `protein_mode`、`threads` を使う | PR-2 |
| Python 版 LOSATN / rings | PYTHON 配下 | PR-3、PR-4 |

GUI の tutorial は UI の用語しか使っていないので、変えない。

**Recipe / scenario**

- `docs/recipes/run_cli_scenarios.py:1170-1198`（`_assert_pinned_losat` が旧 flag を探している）。
- `docs/scenarios/manifest.json:1086-1088`。
- `docs/internal/SCENARIO_EVIDENCE.md`（H-CLI-06/07/08）。これを `tests/test_cli_comparison_how_to_recipe_contracts.py:150-175` が検査している。
- `recipe/meta.yaml` は変えない（`losat ==0.1.0`）。

**Gallery**

- `gbdraw/web/gallery/examples.json` の 5 entry（`hepatoplasmataceae_collinear`、`vibrio-harveyi-group-collinear`、`hepatoplasmataceae_orthogroup`、`BGC0000708-BGC0000713`、`majanivirus_orthogroup`）の `command`。
- session 5 本（plain 1 本、gz 4 本）の `cliInvocation.args`。
- どちらも generator-owned なので再生成する。SVG は byte 一致が条件。

**Tools**

- `tools/prepare_interactive_gallery_assets.py:152-182, 388`
- `tools/reproduce_examples_manifest.py:1311-1474`
- `tools/verify_losat_installation.py:109, 172-176, 342-343`（これを `losat-distribution.yml` が実行する）
- `tools/prepare_issue597_s01_fixture.py:86-117`
- `tools/build_metazoan_mitochondria_comparison_fixture.py`

**Test**

- Python: option 名や field 名を参照する file が 21、関数が約 93。多いのは次の file。
  - `test_protein_colinearity.py`: 31
  - `test_session_request_codec.py`: 10
  - `test_collinearity.py`: 9
  - `test_similarity_alignment_reference.py`: 6
  - `test_api_library_usage.py`: 6
  - `test_session_io.py`: 5
- Python で doc の flag 文字列を検査する test: `test_cli_comparison_how_to_recipe_contracts.py:168-175`、`test_python_tutorial_recipe_contracts.py:85`、`test_cli_reference.py:37`。
- Web: option 名や camelCase の名前を参照する `test()` block は約 11。`linearComparisonPlan` と `generatedProteinComparison` は変えないので、ほかの約 60 block は影響を受けない。
- fixture: gz session 8 本と plain JSON 2 本。wire を変えないので、書き換えは `cliInvocation.args` だけ（0.13.0 由来の fixture は reader の positive fixture として残す）。

**CI workflow**

- workflow の YAML は直接変えない。
- `losat-distribution.yml` は `tools/verify_losat_installation.py` を通して影響を受ける。
- `test.yml` の LOSAT job（938-1007）は、置き換えた bundled binary を使う。
