<!-- Raw design report of workstream W2 (sub-agent output, 2026-09-30, base origin/dev 4c89bab1). Evidence for 01_REMEDIATION_PROPOSAL.md, not Product or architecture authority. -->

# W2 修正提案: record discovery、入力と record の構成、annotation の対象、record identity、座標表示

対象は dev `4c89bab1`（読み取りのみ）。パスはリポジトリ直下から書き、JS は `gbdraw/web/js/...` と省略せずに書く。

## 自分で確かめた事実（scratch: `.../scratchpad/design/w2/`）

- **IN-02**
  - BioPython 1.85 は、空の ACCESSION/VERSION（Prokka 形式）で `contig_1` を返す。ACCESSION だけがあり VERSION が空なら `NC_1` を返す。
  - 現在の fast path は、どちらの場合も `KEYWORDS` を返す。改行をまたがない `[ \t]+` にすると、11 本のベクタすべてで BioPython と一致した。
- **IN-02（新規）** `gbdraw/web/js/app/record-discovery.js:84` の `/^ {2}ORGANISM\s+([^\r\n]+)/m` も改行をまたぐ。
  - ORGANISM 行が空だと、JS は organism を `"Unclassified."`、inferredDefinition を `<i>Unclassified.</i>` とする。Python は `''` を返す。
- **IN-02（新規）** 同じ誤った正規表現が、あと 2 か所にある。
  - `gbdraw/web/js/app/run-analysis.js:772-774`（`parseGenbankRecordsFast`。LOSAT の FASTA 抽出と record 選択 `:735-749, :847-866` に使う）
  - `gbdraw/web/js/app/match-sequences.js:820-822`
  - record-discovery だけを直すと、Prokka 形式の入力で LOSAT 抽出が `Record selector 'contig_1' did not match any record ID.` で失敗する（コードからの判断。未実行）。
- **IN-02 BOM**
  - `TextDecoder` は既定で BOM を取り除くので、JS は 1 record と数える。
  - Pyodide/CLI の `SeqIO.parse` は 0 record を返す。
- **IN-03** 原因は fast path だけではない。
  - Worker helper の `list_gff_fasta_records`（`gbdraw/web/js/app/python-helpers.js:1645-1660`）も FASTA を列挙するだけである。`09c99fb8`（2026-07-29）で `load_gff_fasta` 呼び出しから置き換えられた。
  - canonical loader の record 集合は `merge_gff_fasta_records`（`gbdraw/io/genome.py:105-135`）が決める。BCBio が作る GFF record（feature 行の seqid と、埋め込み `##FASTA` の配列）のうち FASTA と照合できたものを、FASTA 順に並べたものになる。
  - 実例:
    - `gffpair` → `chr, plasmid2`
    - `##sequence-region` だけの行は無視される
    - 埋め込み `##FASTA` がある場合 → `chr, plasmid1(0 features), plasmid2`
  - JS の試作（約 20 行、判断できない入力は decline）は、この 6 ベクタで Python と一致した。
- **FE-05** 型付き API で `plan_circular_batch_request` を確かめた。
  - multi.gb（2 record）で target の record が null → `ValidationError: Annotation target without a record selector is ambiguous for multiple records.`
  - `record_id=TESTB` を明示 → 成功し、item 1 だけに解決する。
  - したがって監査の回避策「target record を消す」は、batch では成り立たない。

---

## IN-02（P1）空の ACCESSION/VERSION で record ID が `KEYWORDS` になる

1. **分類**: IMPLEMENT_EXISTING_AUTHORITY。
   - 根拠 1: PD-OI-044 の mustPreserve「exact source-bound identity」と truthful discovery。
   - 根拠 2: web CLAUDE.md の「Python owns … loading」。record ID は BioPython が決める事実である。
   - BOM は二つに分かれる。
     - (A) fast path を decline させ、Python と同じ失敗にそろえる → 既存 authority。
     - (B) Python の loader（`gbdraw/io/genome.py`）で BOM を受け付ける → 受け付ける入力が変わるので PRODUCT_DECISION_REQUIRED。
     - 推奨は、今は A、B は別途判断。
2. **根本原因**
   - `gbdraw/web/js/app/record-discovery.js:137-138` の `\s+` が改行をまたぎ、次の行の見出しを拾う。
   - 同じ原因が `:84`（ORGANISM）、`:135`・`:144`（LOCUS）、`gbdraw/web/js/app/run-analysis.js:772-774`、`gbdraw/web/js/app/match-sequences.js:820-822` にもある。
   - 起源は `0ce21a32`（2026-08-10、main に含まれる）。監査の「#597 regression」は誤り。
3. **修正案（コード）**
   - owner は `gbdraw/web/js/app/record-discovery.js` とする。
     - GenBank 見出しの正規表現を `[ \t]+` にする。
     - `genbankHeaderIds(chunk)` → `{locus, accession, version, recordId}` を export し、`parseGenBankRecordText` で使う。
   - 重複の削除:
     - `gbdraw/web/js/app/run-analysis.js:772-775` は `genbankHeaderIds(chunk).recordId` に置き換える。
     - `gbdraw/web/js/app/match-sequences.js:820-825` は alias をこの関数の出力から作る。
   - BOM:
     - `discoverSequenceRecords`（`:186-211`）で、native File の先頭 3 byte が `EF BB BF` なら fast path を decline して Worker に回す。
     - production の呼び出し側が渡している冗長な `readText: readFileText`（`gbdraw/web/js/app/app-setup.js:1098-1107, 2818-2819`、`gbdraw/web/js/app/run-analysis.js:1793-1798`）は削除する。読み方の決定を owner に集める。
   - 追加しないもの: 新しいモジュール、Python の変更、error mapping（W3）。
4. **上位の対応（JS の軽量 parser を残すか）**
   - 結論: **残す。ただし authority ではなく最適化と位置づけ、「完全一致か decline か」を契約にする。**
   - 残す理由:
     - PD-OI-044 は「valid native upload の自動 record discovery」「既存 parser/helper 境界」「preview-only Load の Python Worker 0」の保持を求めている。
     - lazy-Worker 規則は discovery のための Worker 起動を許すが、義務づけてはいない。
     - fast path をなくすと、upload のたびに Pyodide の cold start がかかる（監査は Worker 起動に約 8 秒を観測）。これは product effect の変更で、新しい Product Decision が必要になる。
   - 不採用の案:
     - Worker だけで discovery する案: PD が必要になる。
     - fast path と Worker の両方で調べる案: upload 時に Worker を起動するので、fast path の意味がなくなる。
   - parity の強制:
     - 共有ベクタを一つ持つ。新しい fixture は作らず、`tests/fixtures/record_metadata_inference_cases.json` に `discovery` 節を足す。
     - 同じベクタを次の 3 つが実行する。
       - (a) canonical loader（`load_gbks`/`load_gff_fasta`）
       - (b) Worker helper（`python-helpers.js` の Python 部分を exec する。`tests/test_web_feature_metadata.py:28-37` と同じ方法）
       - (c) JS fast path（期待値と一致するか、decline して helper を呼ぶことを assert）
     - web CLAUDE.md の Computation ownership に invariant を 1 行加える。
5. **テスト**
   - Python（`tests/test_record_metadata.py`）: (a)(b) が expected と一致する。
   - Node（`tests/web/record-metadata-inference.test.mjs`）: (c)。
   - ベクタ: Prokka 形式、ACCESSION だけ、GI 付き VERSION、複数 accession、空の ORGANISM、BOM（期待: loader が error、JS は decline）。
   - Playwright: 単一 record の Prokka ファイルで Circular Generate が成功し、出力名が `contig_1.svg` になる。Linear の LOSAT 抽出も成功する。
   - 検出できた assertion: `recordId === 'contig_1'`。
6. **規模・PR**
   - production 3 ファイル、churn 約 35、net 約 +8。
   - Ordinary profile。Review REQUIRED（新しい export）。
   - IN-03 と同じ PR にする（fixture と harness を共有するため）。
7. **依存・リスク**
   - 3 か所の parser を同時に直さないと、LOSAT が壊れる。
   - main にも同じコードがあるので、backport するかを判断する必要がある。

## IN-03（P1）GFF3+FASTA で、GFF feature のない FASTA 配列が NO_MATCH になる

1. **分類**
   - discovery を loader に合わせる部分: IMPLEMENT_EXISTING_AUTHORITY（PD-OI-044 と、Python が loading を持つ規則）。
   - FASTA だけにある配列を、feature のない record として CLI/Web で描くか: 別の PRODUCT_DECISION_REQUIRED。CLI の科学的出力が変わる。PD-OI-018(1) の「every record … of paired GFF3/FASTA source」はどちらにも読める。今は現状（落とす）を推奨する。
2. **根本原因**
   - discovery の 2 経路（JS の FASTA のみ `:213-229` と、Worker helper `python-helpers.js:1645-1660`）が、FASTA の見出しを列挙する。
   - loader は GFF record と FASTA の積を返す。
   - 監査の「Python drops featureless records」は不正確。落ちるのは GFF 行を一つも持たない FASTA 配列で、埋め込み `##FASTA` や `region` 行があれば 0 feature でも残る。
3. **修正案（コード）**
   - 必須: `list_gff_fasta_records` を `load_gff_fasta([gff],[fasta])` の結果（id、`len(seq)`、FASTA 順）にする。`09c99fb8` 以前の意味に戻す。
   - `discoverGffFastaRecords`（`gbdraw/web/js/app/record-discovery.js:213-229`）:
     - GFF の text も読む。`##FASTA` より前の非コメント行の 1 列目と、`##FASTA` 以降の `>` id を集め、FASTA 順で積をとる。
     - 次の場合は decline する: 8 列未満、seqid が `.`、start/end が欠けている、FASTA の id が重複、GFF の seqid が FASTA にない。
   - 代替: GFF+FASTA の fast path を削除する（約 10 行減）。ただし upload 時の遅延が変わるので PD が必要。
   - 追加しないもの: 新しいモジュール、FASTA だけの配列の扱いの変更。
4. **上位の対応**: IN-02 と同じ（exact-or-decline の契約と 3 つの実行者）。
5. **テスト**
   - ベクタ: gffpair、sequence-region だけの行、埋め込み `##FASTA`、`region` 型、GFF の seqid が FASTA にない、GFF の並び順が違う。
   - 3 つの実行者で確かめる。
   - Playwright: `gffpair` で Circular の All records と Linear の自動展開が Generate まで通る。
6. **規模・PR**
   - production 2 ファイル（record-discovery.js 約 +22、python-helpers.js 約 ±8）。
   - IN-02 と合わせて 4 ファイル、churn 約 110、net 約 +35。Ordinary。
7. **依存・リスク**
   - fallback helper で BCBio の parse が走る。大きな GFF で所要時間を測り、PR に記録する。
   - 保存済み Session に FASTA だけの record を指す行があっても、もともと Generate できなかったので、互換性の影響はない。

## IN-04（P2）Linear で crop 済みファイルを複数 record のファイルに置き換えると、古い crop と定義が残る

1. **分類**
   - 複数 record のファイルに置き換えた場合: IMPLEMENT_EXISTING_AUTHORITY。
     - PD-OI-018(1)(2): 「replacing a source updates all of its records」。
     - PD-OI-044: 先頭 record の自動選択は受け入れない。
     - `gbdraw/web/js/services/config.js:2881-2886`: 選択なしの multi-record と crop の組み合わせは常に不正。
   - 単一 record から単一 record への置き換えで状態を残すか: 任意の判断で、outcome は二つ。
     - (A) 置き換えは常に新規 upload と同じ扱い（rotation drafts の purge と一貫する）。
     - (B) 現状どおり残す。
     - 推奨は B。IN-04 の範囲では変えない。
2. **根本原因**
   - `gbdraw/web/js/app/app-setup.js:4117-4123` が record 単位の状態を引き継ぐ。
   - `expandDiscoveredLinearRecords` の `:1049` が region を理由に展開を飛ばす。
   - さらに展開時に `...source` をコピーするので（`:1052-1057`）、definition、record_subtitle、crop が全 record に広がる。
3. **修正案（コード）**
   - owner は `expandDiscoveredLinearRecords`（`gbdraw/web/js/app/app-setup.js:1044-1062`）。
   - 早期 return は `region_record_id` がある場合だけにする。
   - 展開する行は source 単位の値だけで作る: gb/gff/fasta、depth、losat_gencode、file_definition/file_subtitle、`region_record_id: record.value`。
   - record 単位の値（definition、record_subtitle、region_start/end、region_reverse）は捨てる。
   - 追加しないもの: 新しい状態。メッセージ（IN-07/X-01 は W3）。
4. **上位の対応**: 不要。横断的な提案 8 を参照。
5. **テスト**
   - Playwright（`tests/web/linear-multi-record.playwright.spec.js`）:
     - crop と定義を設定したあと 2 record のファイルに置き換える → `linearSeqs.length === 2`、各行の region が null、definition が `''`、Generate が成功する。
     - 単一 record から単一 record への置き換えでは crop が残る（回帰の確認）。
6. **規模・PR**: 1 ファイル、churn 約 12。Ordinary。IN-06 と同じ PR にする（同じ関数群）。
7. **依存・リスク**: 展開前に利用者が選択なしで入れた crop は失われるが、その crop はもともと有効にならない。

## IN-05（P2）Circular→Linear→Circular で Multi-Record Canvas の並び順が消える

1. **分類**: IMPLEMENT_EXISTING_AUTHORITY。OIPC-C06「inactive modes … do not imply deletion」と、コードのコメント `gbdraw/web/js/app/run-analysis.js:1830`。
2. **根本原因**
   - `:1748` の `hasActiveInput` にモードの判定が入っている。
   - そのため Linear で watcher（`gbdraw/web/js/app/watchers.js:586-603`）が動くと、`:1767-1780` が list と `multi_record_positions` を消す。
   - Circular に戻ると、空の配列と merge される。
3. **修正案（コード）**
   - `runCircularRecordRefresh` の先頭に `if (mode.value !== 'circular') return;` を置く。await 後の guard は残す。
   - `:1737` の `mode.value === 'circular' &&` を削除する。
   - `hasActiveInput` をやめて `hasCompleteInput` にする。削除するのは入力がなくなったときだけにする。
   - watcher は変更しない。
4. **上位の対応**: 不要。
5. **テスト**
   - Node（`tests/web/run-analysis-simple-path.test.mjs`。既存の harness `runner.refreshCircularRecordOrder` を使う）。
   - 位置を並べ替える → `mode='linear'` で refresh → `mode='circular'` で refresh。
   - assert: `adv.multi_record_positions` が `#2@1,#1@1,#3@1` のまま、helper と fast path の呼び出しが 0 回。
6. **規模・PR**: 1 ファイル、churn 約 8、net −2。Ordinary。Review 対象外。単独 PR。
7. **依存・リスク**: 「mode 切り替えで list が消える」ことを前提にしたテストがあれば直す（OIPC-C08）。

## IN-06（P2）1 ファイル内の複数 record すべてに、1 番目の record の organism が付く

1. **分類**: PRODUCT_DECISION_REQUIRED。現状は、別の生物名を付けるので科学的に容認できない。
   - (A) record の inferredDefinition がすべて同じときだけ、ファイルの既定値に入れる。違えば空にする（CLI と同じで organism の行なし）。
   - (B) 異なる場合は、各 record の `definition` に自分の値を入れる。表示は正しくなるが override 扱いになり、「Using file default」と Reset の意味が崩れる。
   - (C) record ごとの推定値の層を新設する。Session field が必要になるので YAGNI。
   - **推奨は A**。docs の原則（`docs/REFERENCE/web-app.md:51-53`「a subtitle names one replicon rather than the whole file」、spec `:4168` のコメント）と同じ考え方になる。
2. **根本原因**
   - `gbdraw/web/js/app/app-setup.js:1070-1081` が `records[0].inferredDefinition` をファイルの既定値にする。
   - 監査が挙げた `gbdraw/web/js/app/linear-sources.js:247-255, 271-276` は存在しない（ファイルは 235 行）。該当するのは `:107-145`。
3. **修正案（コード）（A）**: `:1073-1079` で、`new Set(records.map(r => r.inferredDefinition || ''))` の要素が 1 つで空でないときだけ set する。
4. **上位の対応**: 採用後に `docs/REFERENCE/web-app.md:51` を「全 record で共通のとき」と 1 行直す。
5. **テスト**
   - Playwright（`tests/web/linear-multi-record.playwright.spec.js` の `:4168` 付近）。
   - organism が異なる 3 record のファイル → 既定値が `''`。
   - 同じ organism の複数 record → 従来どおり入る。
6. **規模・PR**: 1 ファイル、churn 約 6。IN-04 と同じ PR（PD が決まってから）。
7. **依存・リスク**: 既存の Session の保存値は変わらない（移行不要）。

## IN-08（P2）最初の Generate の前に Save Session が UNKNOWN で失敗する

1. **分類**: IMPLEMENT_EXISTING_AUTHORITY。
   - 保存を拒むこと自体は仕様である: `docs/REFERENCE/web-app.md:210`「Invalid drafts cannot be saved as valid Sessions」、`docs/SESSION_COMPATIBILITY.md` の Session 42 の節。
   - 不具合は案内だけ。案内の扱いは PD-OI-046 に従う。
2. **根本原因**
   - (a)(b) Linear の catalog の issue（`gbdraw/web/js/app/annotations/record-catalog.js:100`）が mapping されていない → W3 の X-01。
   - (c) Save（`gbdraw/web/js/services/config.js:4078-4102`）は、Generate の事前検証（`gbdraw/web/js/app/run-analysis.js:2692-2696`）を通らない。
     - そのため内部 invariant の `Canonical resource record-1-genbank is missing.`（`gbdraw/web/js/services/session-request.js:641`）がそのまま利用者に出る。
3. **修正案（コード）**
   - `gbdraw/web/js/services/config.js` に、active mode の入力の有無を確かめる関数を 1 つ置いて export する。
     - 文言は mapping 済みの既存のもの: `Please upload a GenBank file.` / `GFF3 and FASTA are required.`
   - Save（`:4078` の直前）と Generate（`run-analysis.js:2692-2696` を置き換え）の両方から呼ぶ。
   - `addFile` の例外は内部 invariant として残す。
4. **上位の対応**: 不要。X-01（W3）に依存する。
5. **テスト**
   - Node または Playwright: Circular を active にし入力なし、Linear は読み込み済みで Save → `{status:'error', code:'INPUT_REQUIRED'}`。
   - (a) は W3 の mapping テストで確かめる。
6. **規模・PR**: 2 ファイル、churn 約 16。Ordinary。Review REQUIRED（新しい export）。W3 の X-01 PR と同時か、その後にする。
7. **依存・リスク**: Generate のエラーの順番（depth より先に入力を確かめる）は変えない。

## FE-05（P2）Circular の複数 record（single/batch）で「Selected features」が Generate を止める

1. **分類**: IMPLEMENT_EXISTING_AUTHORITY。
   - Python の typed request が validation の authority である。
   - PD-OI-041 の mustPreserve が「multi-record の record 省略は fatal」と定めている。
   - OIPC-C03 にも反する。
2. **根本原因**
   - `gbdraw/web/js/app/annotations/record-catalog.js:145` の `allowExplicitSelectors = multiRecordCanvas || records.length <= 1` は、Circular batch 導入（`09c99fb8`）より前の規則（`5eff9481`）のまま残っている。
   - batch では、Python が必要とする明示 target を禁止している。record 省略は Python で ambiguity エラーになる（確認済み）ので、**有効な設定が一つもない**。
   - single では、request の record は 1 つなのに catalog が全 record を持っている。
   - `target-actions.js:91-106` が record を結び付けること自体は正しい。
3. **修正案（コード）**
   - `gbdraw/web/js/services/session-request.js:1017-1036` にある Circular の record 集合の解決（`singlePresentationRequested`、`resolveDisambiguatedRecordSelection` を使う）を、純関数として 1 つ export する。buildRecords はその関数を使い、元の inline は削除する。
   - `gbdraw/web/js/app/app-setup.js:1111-1135`（getAnnotationRecordCatalog）はこの関数で request の record 集合を得て、catalog に渡す。
   - `buildCircularCatalog` では `requiresSelection = records.length > 1` とする。
   - 削除:
     - `allowExplicitSelectors` という概念そのもの（`record-catalog.js:128, 145, 152`）
     - `gbdraw/web/js/app/annotations/validation.js:53-55`
     - `gbdraw/web/js/app/annotations/record-selector.js:99-100, 113-114, 137, 147-153`（「Automatic (current output record)」を含む）
     - `gbdraw/web/js/services/error-normalization.js:65, 304` の TARGET_MODE（W3 と調整）
   - CW-01 に従い、catalog は request を組み立てない。
4. **上位の対応**: 不要。Circular の grouping の owner は session-request のまま、catalog がそれに依存する（DIP）。
5. **テスト**
   - Node（`tests/web/annotations.test.mjs:392-402` を反転。OIPC-C08）:
     - batch では target が null だと「Choose a target record」、明示なら `''`。
     - single の catalog は選択した record だけを持つ。
   - Playwright（`tests/web/annotation-resolution.playwright.spec.js`）: single と batch で「Selected features」を使い Generate が成功する。
   - Python: batch で null → ValidationError、明示 → その item にだけ結び付く。
6. **規模・PR**: 6 ファイル、churn 約 70、net 約 −10。Ordinary。Review REQUIRED（新しい export）。単独 PR。
7. **依存・リスク**
   - batch で target が null の既存 Session は、もともと Generate できなかった。
   - single で明示 target を持つ Session は、使えるようになる。

## FE-09（P3）Linear で record ID が重複すると「This feature only」の色が効かない

1. **分類**: PRODUCT_DECISION_REQUIRED。Web が Python の照合できない値を出していること自体は、OIPC-C03 違反で許されない。
   - (A) 常に stable hash を出す。同一の複製 record の両方に効く。
   - (B) Python が rendered instance id（`<hash>_record_N`）で照合できるようにする。色テーブルで受け付ける値が増えるので CLI にも影響する。
   - (C) 衝突する feature では「This feature only」を無効にし、理由を示す。
   - **推奨は A**。コードの削除だけで済み、Python の identity 設計とも一致する。
2. **根本原因**: `gbdraw/web/js/app/feature-editor/rule-actions.js:462-470` が衝突時に svg instance id を出す。`gbdraw/features/selector_values.py:289-292` は stable hash としか照合しない。
3. **修正案（コード）（A）**: `:465-468` の分岐と `generationHashCounts`（`:447-460`）を削除する。
4. **上位の対応**: 不要。
5. **テスト**
   - Node: 重複 record では `getFeatureQualifier` が stable hash を返す。
   - `tests/web/python-rule-parity.playwright.spec.js`: live と Generate の両方で色が反映される。
6. **規模・PR**: 1 ファイル、net 約 −18。Ordinary。
7. **依存・リスク**: A では 2 つの複製の両方の色が変わる。

## FE-11（P3）原点をまたぐ feature や分割された feature の位置と長さが誤る。座標が 0 始まりになる

1. **分類**: IMPLEMENT_EXISTING_AUTHORITY。INSDC の 1 始まりの閉区間という決まりがあり、Python の `location_parts[].display`（`gbdraw/web_support/feature_metadata.py:196-212`）が formatter の owner である。
2. **根本原因**
   - `gbdraw/web/js/app/feature-editor/svg-actions.js:130-136, 175-181` が start/end の外包区間から位置と長さを計算している。
   - drawer（`gbdraw/web/index.html:6606`）は 0 始まりの生の値を出す。
   - 検索の `Start`/`End`（`gbdraw/web/js/app/feature-search/search-core.js:356-357`）も生の値である。
   - 正しい実装は `search-core.js:251-266` にすでにある。
   - 監査にない範囲（コードからの判断）: Interactive SVG の `gbdraw/web/js/services/standalone-interactivity-assets.js:2621-2624, 4475-4482, 5765-5770` も同じ誤り。
3. **修正案（コード）**
   - `search-core.js` の `buildFeatureLocation` を export し、part の合計で長さを返す関数を横に置く。
   - svg-actions は自前の 2 関数を削除して、これを使う。
   - drawer は、既存の feature-editor API から公開した formatter を使う。
   - 検索の Start/End の項目は削除する（Location と重複し、値も誤っている）。
   - standalone は別 runtime で import できないので、同じ規則で直し、共有ベクタで拘束する。
4. **上位の対応**: 不要。
5. **テスト**
   - Node:
     - parts `[(3900,4000),(0,200)]` → `3901..4000, 1..200`、`300 bp`
     - 0 始まりの `300` → `301..600`
     - Location に `300` を入れても `301..600` に一致しない
   - Interactive SVG の Playwright で popup の文字列を確かめる。
6. **規模・PR**
   - Web 側 4 ファイル、churn 約 40。Review REQUIRED（新しい export）。
   - standalone は別 PR にする。Gallery の例の SVG と docs の画像にこの script が埋め込まれているため、Gallery refresh が必要になる。
7. **依存・リスク**: 旧 catalog で `display` のないものは、既存の fallback に任せる。

## SE-10（P3）Reset Settings が Linear の record ごとの状態を残す（仕様が曖昧）

1. **分類**: PRODUCT_DECISION_REQUIRED。
   - (A) Linear の record 単位の表示状態（definition、record_subtitle、region_start/end、region_reverse）と alignment plan を初期化する。
     - 残すもの: files、展開された行の `region_record_id`、file defaults、depth。
     - Circular の form の初期化、および `recordDisplayDrafts` の初期化（`gbdraw/web/js/services/reset.js:155`）と一貫する。
   - (B) Circular 側も record 単位の値を残すようにして、両モードをそろえる。
   - (C) 現状のまま docs に明記する。
   - **推奨は A**。
2. **根本原因**: `gbdraw/web/js/services/reset.js:146-166` は `form` と `adv` だけを置き換え、`linearSeqs` と `similarityAlignmentPlan` には触れない。
3. **修正案（コード）（A）**
   - `gbdraw/web/js/app/app-setup.js:3224-3240` の Reset の checkpoint 内で、既存の `applyLinearSeqMutation` を通して上の項目を初期化する。
   - alignment は owner の遷移 `clearCommittedPlan`（`gbdraw/web/js/app/similarity-alignment.js:940-955`）で消す。ref への直接代入はしない。
4. **上位の対応**: 採用後、docs の Reset の範囲を 1 行更新する。
5. **テスト**: Playwright で Gallery の BGC を読み込み Reset → A の項目が空になり、files と selector が残り、Undo で戻る。
6. **規模・PR**: 1〜2 ファイル、churn 約 15。PD が決まってから。
7. **依存・リスク**: Reset 後の alignment と comparison `none` の組み合わせで Generate がどうなるかは未確認。

---

## 横断的な提案

1. **fast path は最適化であって authority ではない。** Python が持つ事実を JS で先に出すなら、canonical の結果と完全に一致するか、decline するかのどちらかにする。1 つの共有ベクタ（`tests/fixtures/record_metadata_inference_cases.json` を拡張）を、loader、Worker helper、JS の 3 者で実行して強制する。web CLAUDE.md の Computation ownership に 1 行加える（IN-02、IN-03、BOM、ORGANISM）。
2. **JS の GenBank 見出し parser は 1 つにする**（record-discovery）。`tests/web/architecture-contracts.test.mjs` に、`/^VERSION` と `/^ACCESSION` の正規表現が他のファイルに現れないことを確かめる assertion を置く。
3. **「Generate が読み込むもの」に答える Worker helper は、canonical loader を呼ぶ。** `SeqIO` で意味を作り直さない。helper の出力 == loader の出力を Python テストで確かめる（IN-03）。
4. **行単位のファイル形式では、`m` flag の正規表現に `\s` を使わず `[ \t]` を使う。** 空の見出し欄をベクタに必ず入れる。
5. **Web の事前検証は、canonical request の record 集合（session-request.js）と Python の binding 規則から導く。** Web だけの可否 flag は削除する。Circular の single/grid/batch と target の null/明示の組み合わせを、request→Python planner の parity テストで確かめる（FE-05）。
6. **Web が出す selector の値は、Python が照合できるものに限る**（OIPC-C03）。`tests/web/python-rule-parity.playwright.spec.js` を、出力するすべての selector の種類に広げる（FE-09）。
7. **表示座標は Python の `location_parts[].display` だけから作る。** JS の formatter は 1 つにし、standalone はベクタで拘束する。表示用モジュールに `start + 1` や `end - start` の外包計算が残らないことを test で確かめる（FE-11）。
8. **record に結び付く draft は、モードごとに owner の遷移表（upload、replace、discover、mode 切り替え、Reset、Load）で扱う。** 状態ごとに残すか消すかを明記する。OIPC-C06 は、per-mode の全 field の往復テストで強制する（IN-04、IN-05、SE-10）。
9. **Save と Generate は同じ入力前提の owner を使う。** 内部 invariant の文言は利用者に届かないようにし、その mapping は W3 に任せる（IN-08）。
10. **file の既定値は、source 全体で成り立つ事実からだけ推定する**（IN-06。subtitle ですでに文書化されている原則を organism にも適用する）。

## 監査の分類・根拠で誤りや不足がある点

- **IN-03**
  - Worker helper（`python-helpers.js:1645-1660`）も同じ原因なので、fast path だけを直しても解決しない。
  - 「Python drops featureless records」は不正確（上記参照）。
  - 「#597 regression」は誤り。起源は `09c99fb8` と `0ce21a32` で、どちらも main に含まれる。
- **IN-02**
  - 範囲が足りない: 同じ正規表現が `run-analysis.js:772-774` と `match-sequences.js:820-822` にある。
  - 同じ種類の新たな不具合がある: ORGANISM の `:84`。
  - 起源は #597 ではなく `0ce21a32`。
  - BOM は Web 固有の読み込みの不具合ではない。CLI も失敗し、不具合は「1 record inspected」と Python と違う結果を出すことにある。
- **IN-06**: 行番号 `linear-sources.js:247-255, 271-276` は存在しない。この organism の推定は Web だけの機能で、CLI の Linear は organism を推定しない。
- **FE-05**
  - 回避策「target record を消す」は batch では無効（Python が ambiguity で拒否することを確認した）。
  - 原因は record の結び付け（`target-actions.js`）ではなく、catalog の規則（`record-catalog.js:145`）である。
- **IN-08**: Save を拒むことは文書化された仕様である。不具合は X-01 の案内と、(c) の内部 invariant が表に出ることに限られる。Save の機能自体の不具合ではない。
- **FE-11**: Interactive SVG の出力と Gallery の成果物にも同じ誤りがある（コードからの判断）。
- **FE-09**: 照合の不具合というより、Python が消費できない値を Web が受け付けていること（OIPC-C03）である。仕様上の穴として扱うのが妥当。
