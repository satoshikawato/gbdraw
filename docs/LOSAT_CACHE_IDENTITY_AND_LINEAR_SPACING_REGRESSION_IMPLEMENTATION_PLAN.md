# LOSAT Cache Identity and Linear Spacing Regression Implementation Plan

- 作成日: 2026-07-21
- 対象: Web Linear の LOSATP cache identity、session schema v2 読み込み、Linear record spacing
- 現行基準: session version 34、`renderRequest.schema == 3`、raw LOSAT cache schema 2、derived cache schema 1
- 変更後: session version 35、`renderRequest.schema == 3`、raw protein LOSATP cache schema 3、derived cache schema 2
- 非対象cache: nucleotide LOSAT は raw schema 2 を維持し、session内のschema 2/3混在を正式に扱う
- 状態: 実装済み、最終検証中
- 関連設計: [`LINEAR_TRACK_OCCUPANCY_LAYOUT_IMPLEMENTATION_PLAN.md`](LINEAR_TRACK_OCCUPANCY_LAYOUT_IMPLEMENTATION_PLAN.md)、[`PYTHON_SESSION_CANONICAL_REQUEST_PLAN.md`](PYTHON_SESSION_CANONICAL_REQUEST_PLAN.md)

## 0. 結論

今回確認された現象は、原因の異なる二つの回帰である。

1. **Workstream A — LOSAT identity/cache**
   - Web LOSATP の QUERY/SUBJECT ID が、ファイルの `name`、`size`、`lastModified` から作る
     record hash に依存している。
   - schema v2 session の復元で `File` metadata が変わると、配列が同じでも FASTA header と
     FASTA 全文 hash が変わり、保存済み raw cache を再利用できない。
   - QUERY/SUBJECT は実レコードと実 feature を示す安定 ID に変更し、cache fingerprint と
     user-visible ID と runtime transport ID の責務を分離する。
2. **Workstream B — Linear spacing**
   - definition、record body、comparison exclusion の独立した占有制約を `canvas_band` へ union
     した後、さらに `comparison_height` を加えている。
   - definition は左列、comparison は plot 列にあるため、同じ Y 投影だけで加算すると余白を
     二重計上する。
   - 独立した collision domain ごとに必要な axis 間隔を求め、最後に `max()` で合成する。

両 Workstream は Phase 0 の実在 session fixture と最終 browser flow だけを共有する。LOSAT を
直すために layout を変更せず、layout を直すために cache/session identity を変更しない。

次の対症療法は採用しない。

- `Auto 60 px` を小さくして過剰 spacing を隠す。
- hash 検証を省略して schema v2 cache を無条件に再利用する。
- `lastModified` だけを固定し、metadata 依存の cache identity を残す。
- 旧レイアウトを長期 flag として併存させる。

## 1. 再現結果と固定する baseline

対象は、履歴上の実在する BGC0000708–BGC0000713 schema v2 session とする。現行 schema 3
session を schema 2 と書き換えただけの synthetic fixture は、migration の証拠にしない。

| 観測項目 | 保存済み schema v2 | 現行コードで Load → Generate | 修正後の条件 |
|---|---:|---:|---:|
| LOSATP comparison pair | 25 | 25 | 25 |
| raw cache hit | — | 0 | 25 |
| raw cache miss | — | 25 | 0 |
| LOSAT job | — | 25 | 0 |
| record axis 間隔 | 105 px | 168.450520833 px | 独立制約の `max()` |
| 最初の comparison path の縦幅 | 50 px | 128.950520833 px | exclusion edge 間の契約値 |
| SVG viewBox height | 850.922 px | 1045.277 px | painted content を過不足なく包含 |

105 px 自体は font metric や表示設定に依存するため、修正後の magic number にはしない。
受入判定は band から計算した式、非衝突、非 clipping で行う。

## 2. 根因

### 2.1 LOSATP cache miss

現行 Web path は概ね次の連鎖を持つ。

```text
File {name, size, lastModified}
  -> recordInstanceKey
  -> p_r_<hash>_<start>_<end>_<strand>_<aa12>
  -> FASTA header を含む全文 hash
  -> raw LOSAT cache key
```

`p_r_...` は cache key そのものではなく、LOSAT に渡す内部 protein ID である。実際の
`protein_id`、`locus_tag`、`featureSvgId`、location は別の `protein_map` にしか存在しない。
raw TSV だけでは、この内部 ID を実 feature へ戻せない。

schema v2 promotion/materialization は埋め込み resource から `File` を再構築する。元 metadata を
完全には保存できず、特に `lastModified: 0` を `entry.lastModified || Date.now()` で現在時刻へ
置換する path もある。このため biological content が同じでも `recordInstanceKey` が変わる。

### 2.2 Linear row spacing

comparison がある場合、現行 `_build_linear_record_vertical_plans()` は
`comparison_exclusion_band` を起点に definition band を `canvas_band` へ union する。その後、
one-record-per-row path は次の式を使う。

```text
axis gap = current.canvas_bottom_extent
         + comparison_height
         + next.canvas_top_extent
```

BGC case では definition が各 canvas extent を支配するため、概ね次になる。

```text
54.2252604165 + 60 + 54.2252604165 = 168.450520833 px
```

既存の Linear occupancy 設計は、次を別々に保持すると定めている。

- `record_body_band`: track、feature、label、annotation の配置用。
- `comparison_exclusion_band`: comparison endpoint と corridor の配置用。
- `canvas_band`: definition を含む最終 canvas enclosure 用。

現行実装はこの境界を row spacing の段階で失っている。

## 3. 変更後の不変条件

### 3.1 Identity layer を分離する

| Layer | 役割 | 許可する依存 |
|---|---|---|
| Protein-set identity | LOSATへ渡す解析対象protein集合のcontent identity | canonical順のfeature analysis IDと完全なAA digest |
| Record analysis identity | 選択された biological record の解析同一性 | record ID/accession、選択領域、protein-set identity |
| Record instance identity | 同じrecordを複数配置した時のsession内一意性 | 永続canonical record key。File metadataには依存しない |
| Feature analysis identity | 一つの CDS feature をrecord内で一意に示すmachine key | feature type、全 location parts、strand、必要時の永続同位置ordinal |
| Source/display identity | 人が読める feature 名。machine identityには使わない | `protein_id`、`locus_tag`、GFF `ID`、location label |
| LOSAT transport ID | FASTA/LOSAT が受理する一意 token | escaped display identity、record-instance identity、feature analysis ID |
| Raw cache fingerprint | LOSAT 実行結果の再利用可否 | ordered protein-set identity、transport mapを含むstable record-instance binding、program、意味のある LOSAT args |
| Derived cache fingerprint | genomic link/orthogroup/描画 payload の再利用可否 | raw key、feature mapping、view transform、filter/converter settings |
| UI pair identity | editor と SVG の comparison 行を結ぶ | canonical record key と pair direction |

次を biological identity、transport ID、raw cache fingerprint に含めない。

- upload filename
- resource filename/resource ID
- file size を単独の identity とした値
- `lastModified`
- session を開いた時刻
- temporary Pyodide path
- 永続化されない現在の配列 index だけから作る token

### 3.2 QUERY/SUBJECT ID contract

新規 LOSATP job の QUERY/SUBJECT は、machine identity、display alias、transport encoding を次の
三段階で作る。qualifierが一意かどうかによってmachine identityの作り方を変えない。

1. `featureAnalysisId` は次のcanonical JSONのSHA-256を`f_<64 lowercase hex>`で表したものとする。

   ```json
   {
     "featureType": "CDS",
     "locationOperator": "join",
     "locationParts": [[7438, 8458, 1]],
     "strand": 1,
     "sameLocationOrdinal": 1
   }
   ```

   canonical JSONはkeyを辞書順にし、空白なしのUTF-8、integer座標でserializeする。
   `locationParts`は全partを正規化した`(start, end, strand)`順で保存し、location operatorも別fieldで
   hashへ含める。座標は0-based、end-exclusiveに統一する。`sameLocationOrdinal`は同じ
   type/location operator/location parts/strandを持つfeature群を、正規化した
   source ID、完全AA digest、永続化するsource feature positionの順でsortして割り当てる。ordinalは
   manifestへ保存し、現在の配列indexだけから再計算しない。
2. `displayAlias` は、空でない `protein_id`、`locus_tag`、GFF `ID`、
   `CDS:<full-location-parts>:<strand>` の順で最初の値を選ぶ。重複していても選択規則を変えず、
   一意性は`featureAnalysisId`が保証する。
3. transport ID は次の固定grammarで作る。

   ```text
   <record-source-id>@<record-instance-id>|<display-alias>~<feature-analysis-id>
   ```

   各text fieldはUnicode NFCへ正規化したUTF-8 byte列とし、`[A-Za-z0-9._-]`以外を大文字16進の
   `%HH`でpercent-encodeする。したがって区切り文字`@`、`|`、`~`、`%`、空白、control characterは
   field内で必ずescapeされる。lossyな置換や切り詰めは行わず、生成tokenがLOSATの非空白ID制約を
   満たすことをvalidatorで確認する。このgrammarとcanonical JSONはPythonのpure functionを唯一の
   ownerとし、WebはPyodide helper経由で同じ実装を呼ぶ。golden fixtureは仕様の代用ではなく、
   Python/Web境界でのserialization driftを検出するために使う。

同じaccession/record IDを複数ロードした場合は、File metadataではなくsessionに保存されたcanonical
record keyを`record-instance-id`として使う。欠損・重複の解決にupload timestampや永続化されない
配列indexを使わない。

例:

```text
BGC0000708@record-1|CAG38695.1~f_2f...
BGC0000708@record-1|CDS%3A7438-8458%3A%2B~f_91...
```

raw TSV の QUERY/SUBJECT だけでrecord instanceと人が読めるsource identityを判別できることを
必須にする。完全なlocation、元の未escape値、`featureSvgId`、AA digestのauthorityはmanifestとし、
transport IDへ同じ情報を重複保存しない。combined `protein_map` のkeyは全recordを通じて一意で
なければならず、全transport IDをmanifestから一意なfeatureへ解決できなければならない。

### 3.3 Raw cache と derived cache の invalidation

| 変更 | Raw LOSAT | Derived payload |
|---|---|---|
| filename / `lastModified` / resource rename のみ | reuse | reuse |
| save → load で同じ biological input | reuse | reuse |
| record の表示順・定義文・subtitle | reuse | rebuild if indices/presentation are embedded |
| display reverse-complement のみ | reuse | rebuild for view transform |
| product/gene/note だけの変更 | reuse | rebuild |
| `protein_id` / `locus_tag` / GFF `ID` の変更 | miss because transport binding changes | rebuild |
| selector/region/visibility により protein set が変化 | miss | rebuild |
| feature analysis identity または AA 配列が変化 | miss | rebuild |
| genetic code により翻訳結果が変化 | miss | rebuild |
| LOSAT program/outfmt/search args が変化 | miss | rebuild |
| post-filter/orthogroup/collinearity option だけの変更 | reuse | miss |
| query/subject direction が変化 | direction-aware lookup | rebuild |

raw reuse は「設定名が同じ」ではなく、canonical protein set と search semantics が同じことから
判断する。false hit より明示的な miss を優先する。

### 3.4 Linear collision domain

汎用2D collision solverは導入せず、隣接row間のclearanceだけを解くpure helperを一つ置く。
solverの入力は次の最小構造とする。

```text
CollisionBand {
  kind: "body" | "comparison" | "definition",
  x_start,
  x_end,
  top_y,
  bottom_y
}
```

隣接するcurrent/next rowについて、X区間が正の幅で交差し、下表のpair policyを持つband pairだけを
比較する。kind pairごとの必要間隔と最終axis gapは次で固定する。

```text
pairs = {
  (current_band, next_band)
  where intervals_overlap(current_band.x, next_band.x)
    and clear_gap(current_band.kind, next_band.kind, boundary) is not None
}

axis_gap = max(
  [0] + [
    current_band.bottom_y
    - next_band.top_y
    + clear_gap(current_band.kind, next_band.kind, boundary)
    for each pair in pairs
  ]
)
```

`clear_gap(kind pair, boundary)`は次の一か所で定義する。

| Current kind | Next kind | clear gap | Active condition |
|---|---|---:|---|
| `body` | `body` | `ordinary_row_gap` | 常に |
| `comparison` | `comparison` | `comparison_height` | そのrow境界を実際に跨ぐcomparisonがある時だけ |
| `definition` | `definition` | `definition_clear_gap` | 常に |
| `definition` | `body` | `definition_clear_gap` | Xが交差する時 |
| `body` | `definition` | `definition_clear_gap` | Xが交差する時 |

`definition_clear_gap = max(1 px, 0.5 * vertical_padding)`とする。body reserveはcomparison paintを
包含するため、`definition`対`comparison`と`body`対`comparison`のcross-kind pairは追加せず、
同じ物理制約を重複計算しない。各recordで`comparison_exclusion_band`が`record_body_band`に包含
されることをsolver入力validationで必須にする。

- `comparison_height` は comparison exclusion edge 間の**最小 clear corridor**とする。
- `record_body_band`と`comparison_exclusion_band`は各recordまたは連続するplot X区間ごとのband set
  としてrowへ集約し、離れたX区間を一つのbounding intervalへ潰さない。comparisonのない境界では
  `plot-comparison` bandを渡さない。
- definitionの`local_band`と`row_band`は、それぞれ実際のX範囲を持つ`definition` bandとして保持する。
  multi-recordでも一つの非負canvas extentへ潰さず、Xが交差する全pairの最大値を取る。
- 左definition列とplot列のようにXが交差しないpairはspacingへ影響させない。multi-recordのlocal
  headerのようにdefinitionとplotのXが交差するpairは`definition`/`body` policyでclearanceを保証する。
- `canvas_band`は配置後の全painted contentのenclosureにのみ使い、axis gapの入力にしない。
- single-rowとmulti-record-rowは同じ`required_axis_gap()`を呼び、別のspacing式を持たない。

## 4. Target cache/session architecture

### 4.1 Protein set、record analysis、record-instance binding

protein content、record固有のanalysis context、session内の配置instanceを別objectとして保存する。
`proteinSets`はcontent-addressed objectであり、record ID/accessionを所有しない。

```json
{
  "proteinSets": {
    "sha256:...": {
      "schema": 1,
      "proteins": [
        {
          "featureAnalysisId": "f_2f...",
          "sameLocationOrdinal": 1,
          "sourceFeaturePosition": 17,
          "sourceProteinId": "CAG38695.1",
          "locusTag": null,
          "gffId": null,
          "featureSvgId": "fed46a3a6",
          "featureType": "CDS",
          "locationOperator": "join",
          "locationParts": [[7438, 8458, 1]],
          "strand": 1,
          "aaSha256": "..."
        }
      ]
    }
  },
  "recordAnalyses": {
    "sha256:...": {
      "schema": 1,
      "recordSourceId": "BGC0000708",
      "selector": null,
      "region": null,
      "proteinSetHash": "sha256:..."
    }
  },
  "recordInstances": {
    "record-1": {
      "recordAnalysisId": "sha256:...",
      "transportIds": {
        "f_2f...": "BGC0000708@record-1|CAG38695.1~f_2f..."
      }
    }
  }
}
```

`proteinSetHash` は`featureAnalysisId`でsortしたcanonical protein listの`featureAnalysisId`と完全な
AA SHA-256から作る。`recordAnalysisId`はrecord source identity、selector/region、
`proteinSetHash`のcanonical JSONから作る。いずれもFASTA description、File metadata、record
instance key、JSON object insertion orderに依存させない。

protein-set manifest はpairごとに複製せずhashでdeduplicateする。同じprotein setを持つ異なる
record analysisも同じcontent objectを参照でき、record固有情報を上書きしない。同じrecordを二回
配置した場合は、一つのrecord analysisへ二つのrecord-instance bindingを作る。binding hashは
canonical `recordInstanceKey`、`recordAnalysisId`、全`transportIds` mapから作る。raw entryとraw cache keyはquery/subjectの
protein-set hashに加えてこのbinding hashを含み、保存TSV内のtransport IDがそのbindingと一致する
ことを検証する。初回実装では別record instance間でraw TSVを共有しない。これによりcombined
`protein_map`のtransport IDを全recordで一意に保ち、別rowのIDを誤利用しない。

### 4.2 Schema version

変更後の version owner は次のようにする。

| Owner | Current | New | 方針 |
|---|---:|---:|---|
| Session envelope | 34 | 35 | Python/Web writer を同時更新 |
| Canonical `renderRequest` | 3 | 3 | render request semantics は変更しない |
| Raw protein LOSATP cache | 2 | 3 | protein manifest と protein-set hash を導入 |
| Raw nucleotide LOSAT cache | 2 | 2 | 現行contractを維持し、protein entryと混在可能にする |
| Derived LOSATP cache | 1 | 2 | mapping/view/converter fingerprint を分離 |
| Protein identity manifest | — | 1 | protein set、record analysis、record-instance bindingを分離して永続化 |
| Legacy protein cache candidate envelope | — | 1 | 旧sessionからimportした未検証schema 2 entryだけを隔離して保持 |

current writer はversion 35のみを書く。protein LOSATP writerはschema 3だけを書き、nucleotide
LOSAT writerはschema 2を維持する。readerはversion 27–34の既存対応を維持し、protein schema 2
raw cacheを`legacyArtifacts.proteinRawCandidates`へ移してlazy migrationの候補として保持する。
v35 writerはschema 2 entryをcurrent `losatCache`へ書かないが、未検証candidateはlegacy envelope
schema 1としてlosslessにround-tripする。entryの`program`/`flow`/identity kindでschema ownerを
判別し、unknown/newer schemaを推測で読み替えない。

### 4.3 全 session file の更新

同梱 current session は一部だけでなく全件を version 35 へ再生成する。

| Scope | 件数 | 更新方法 |
|---|---:|---|
| `gbdraw/web/gallery/sessions/*.gbdraw-session.json[.gz]` | 11 | `tools/refresh_gallery_sessions.py` の canonical path |
| `tests/test_inputs/*.gbdraw-session.json` | 2 | 各 fixture の canonical save/generator path |
| 固定 legacy migration fixture | 1 | schema v2 のまま保持し、current inventory から除外 |

更新後は次を満たす。

- Gallery と通常の test input session はすべて `version == CURRENT_SESSION_VERSION == 35`。
- Linear session は現行の spacing、track slot、definition、comparison 設定を欠落なく保持する。
- LOSATP sessionはprotein raw schema 3、derived schema 2、manifest schema 1を保持する。
- canonical pathで再生成したcurrent sessionは`legacyArtifacts.proteinRawCandidates`を持たない。
- nucleotide LOSAT/BLAST cacheはraw schema 2を維持し、protein manifest migrationの対象にしない。
- `.json` と `.json.gz` は同じ validator を通す。
- file list の重複定義と実ファイルの差を inventory test で検出する。

layout 変更も生成 SVG に影響するため、最終更新は session JSON だけで終えない。review 済みの
session から Gallery source SVG、example SVG、thumbnail、`examples.json` を公式 tool で更新する。
チュートリアル screenshot は実際に旧 spacing が写っているものだけを、Gallery media skill の
手順で再撮影する。

外部ユーザーが保有する旧sessionは一括更新できないため、version 35 readerのverified lazy
migrationで開く。repository内の旧session envelope/render schemaとprotein raw schema 2は
migration fixture以外に残さない。current v35 session内のnucleotide raw schema 2はこの制約の
例外とする。外部旧sessionをGenerate前にv35として再保存した場合だけ、schema 2 entryはcurrent
cacheではなく明示的なlegacy candidate envelope内に残ることを許可する。

## 5. Schema 2 protein raw cache の verified lazy migration

import 時点では Pyodide protein extraction が完了していないため、ID変換は import 中ではなく、
最初の Generate で protein manifest を作った後に行う。schema 2 entryをcurrent cache mapへ混在
させず、`legacyArtifacts.proteinRawCandidates`のread-only candidateとして隔離する。

candidateは次の小さい状態機械で扱う。

| State | 意味 | Save時の扱い |
|---|---|---|
| `pending` | 未検証。Generate前または検証待ち | original entryをlegacy envelopeへlosslessに保存 |
| `promoted` | schema 3 copyの作成とvalidationが完了 | schema 3 entryを保存し、legacy candidateは次snapshotから除去 |
| `rejected` | 現在のinputでは検証不能または不一致 | 理由とoriginal entryをlegacy envelopeへ保存し、current cache hitには使わない |

```json
{
  "legacyArtifacts": {
    "proteinRawCandidates": {
      "schema": 1,
      "entries": [
        {
          "state": "pending",
          "originalEntry": {"schema": 2, "kind": "raw-losat", "key": "...", "text": "..."},
          "rejectionReason": null
        }
      ]
    }
  }
}
```

Load直後にGenerateせずSaveしても、`pending` candidateを失ってはならない。SaveのためだけにPyodide
抽出を起動せず、Saveを禁止もしない。current `losatCache`のvalidator/serializerはschema 3 protein
entryとschema 2 nucleotide entryだけを扱い、legacy schema 2 proteinの解釈はmigration moduleだけが
所有する。

1. schema 2 entry を import 時に破棄せず、`pending` legacy candidate として保持する。
2. 現在の record から schema 1 protein manifest を抽出する。
3. legacy QUERY/SUBJECT の `record token + start/end/strand/aa12` を解析する。
4. `losatDerivedCache` / `orthogroupState` に保存された `recordIndex`、`sourceProteinId`、
   `featureSvgId`、location を補助 evidence として使う。
5. 現在の feature と legacy ID が一対一対応することを確認する。
6. legacy record token を使って旧 FASTA header を再構成し、保存済み
   `queryCanonicalHash` / `subjectCanonicalHash` と完全一致することを確認する。
7. program、outfmt、direction、search args も一致した entry だけ、TSV 第1・第2列を新transport ID
   へ変換する。
8. schema 3 key、protein-set hash、manifest reference で新 entry をcopy-on-writeで作る。
9. derived schema 1 は blind reuse せず破棄し、再利用した raw TSV と現在の manifest から
   derived schema 2 を再構築する。
10. schema 3 entryとderived schema 2のvalidationが成功した後だけcandidateを`promoted`にする。
11. 全検証が終わるまで旧 entry を削除・上書きしない。失敗時は`rejected`として理由を保持する。

次の場合は、その pair だけ cache miss として LOSAT を再実行し、理由を UI/console に残す。

- legacy ID が欠損または解析不能。
- empty result で record token を他 entry/metadata から一意に復元できない。
- compound location または重複 feature が曖昧。
- AA digest、FASTA hash、program、args のいずれかが不一致。
- manifest が同じ legacy ID を複数 feature へ解決する。
- session cache が破損している。

BGC schema v2 fixture は保存済み derived/orthogroup metadata が十分であり、25 pair 全件を migration
できることを acceptance test にする。全 legacy session の無条件 hit は保証しない。

別のv35 readerで`pending`または`rejected` candidateを再読込できること、Load → Save → Load →
Generateでも元のcandidateと25/25 reuseを維持することをsession contractに含める。repository内で
canonical生成するcurrent sessionにはlegacy candidateを残さない。

`lastModified: 0` は有効値として `??` または明示判定で復元し、`Date.now()` に置換しない。ただし
この修正は migration の補助であり、新 cache identity 自体は `lastModified` 非依存とする。

## 6. 実装 phase

### Phase 0: Characterization と failing regression tests

- 履歴上の実在 BGC schema v2 session を最小化・gzip 化し、専用 legacy fixture に固定する。
- fixture metadata/testで `source commit=c64ff8c`、`session.version=33`、
  `renderRequest.schema=2`、protein raw schema 2、derived schema 1を固定する。
- 保存raw entryは34件、対象Generate pairは25件であることを別々にassertする。
- current schema 3 を schema 2 と relabel する既存 test と区別する。
- Load → Generate の LOSAT timing を固定する。
  - `totalPairs=25`
  - `cacheHits=0`
  - `cacheMisses=25`
  - `uniqueJobs=25`
- 保存 SVG と現行 SVG から axis Y、comparison path Y、viewBox を抽出する regression test を追加する。
- Phase 0 では production code、current session、reference SVG を変更しない。

### Phase A1: Stable protein identity と manifest

- Python の CDS extraction を identity の single owner にする。
- Web helper 内に複製された `p_r_<metadata hash>...` 生成を削除する。
- canonical location、同位置ordinal、feature analysis ID、display alias、percent encodingを一つの
  pure Python moduleで実装する。source qualifierの重複有無でmachine IDの分岐を作らない。
- `ProteinSet`をrecord非依存のcontent objectにし、`RecordAnalysis`と`RecordInstanceBinding`へ
  record固有情報とsession配置情報を分離する。
- Web/CLIとも全input rowのrecord-instance keyを先に確定し、それをshared extractorへ渡す。
  Webのper-record extractionとCLIのbatch extractionでuniqueness scopeを変えない。
- Webは既存row `uid` / canonical `renderRequest.records[].recordKey`、Pythonは既存
  `gbdraw_record_key` / canonical record keyを再利用し、別のrandom ID体系を追加しない。
- identical recordを複数配置したfixtureでtransport IDは別、protein-set hashは同じ、raw cache
  keyはbindingごとに別になることを確認する。
- 異なるrecord analysisが同じprotein setを持つfixtureで、一つの`ProteinSet`を共有しつつ
  `recordAnalysisId`とbindingを上書きしないことを確認する。
- Python ownerの出力をWebがbyte-identicalにserializeするgolden testを追加する。Webに第二の
  ID/hash実装を作らない。
- reserved delimiter、空白、control character、Unicode、duplicate qualifier、compound location、
  同位置featureを含むproperty/golden testを追加する。
- filename、mtime、resource name を変えても manifest と protein-set hash が変わらないことを確認する。

### Phase A2: Raw/derived cache schema と lazy migration

- cache key construction/lookup/promotion を `run-analysis.js` から focused module へ移す。
- protein raw schema 3、nucleotide raw schema 2、derived schema 2のdiscriminated validator、
  serializer、pruningを実装する。
- protein schema 2はdual-readし、proteinの新規書き込みはschema 3に限定する。nucleotideの
  新規書き込みはschema 2を維持する。
- legacy protein schema 2 entryはcurrent cache mapへ入れず、legacy candidate envelope
  `pending/promoted/rejected`の状態機械で管理する。
- schema 3 raw keyにquery/subjectのstable binding hashを含め、cache hit時にTSV内の全transport IDが
  現在のbindingに属することを検証する。
- `(query, subject)` はdirectional keyとし、初回実装ではreverse pairをdirect hit扱いしない。
  将来swap reuseを行う場合もTSV列swapと全metadata検証を別contractにする。
- verified migration を pure mapping/validation と runtime orchestration に分ける。
- Load → Save → LoadをGenerate前に行っても`pending` candidateがlosslessに残ることを固定する。
- promotionはcopy-on-writeとし、schema 3 rawとschema 2 derivedのvalidation完了後だけlegacy
  candidateを除去する。
- migration failure は pair-local miss にし、cache 全体を破棄しない。
- blastn、tblastx、circular conservation の schema 2 path を regression test で維持する。

### Phase A3: Version 35 codec/reader と TSV export

- Python/Webにversion 35 codec/validator/reader testを追加するが、current writer constantはまだ
  切り替えない。
- protein set、record analysis、record-instance binding、legacy candidate envelopeをsession artifact
  としてserialize/restore/reset/snapshotする。
- v35 current cacheへprotein schema 2を書かず、import由来の未検証entryだけをlegacy envelopeで
  round-tripすることをvalidatorで分離する。
- `Save Raw LOSAT TSV` の QUERY/SUBJECT が normative encodingに従うreadable transport IDになることを確認する。
- manifest から全 ID を実 feature へ100%解決できることを検証する。
- derived payload、orthogroup member/edge/path、selection、editor override、feature/result metadataの
  全protein referenceを新IDへ解決するinventory/rewriteを追加する。
- orthogroupのname/description overrideとselected alignment featureは、旧ID文字列ではなく
  resolved feature identityにより意味を保持する。
- version 27–34 read、version 35 save/load、gzip round-trip を固定する。
- compatibility matrix と release note に cache artifact の境界を記録する。

### Phase B1: Pure vertical clearance solver

- `kind`、X interval、signed Y edgeを持つ`CollisionBand`と、Xが交差するeligible kind pairの
  最大clearanceを返す`required_axis_gap()`をlayout層へ追加する。
- comparison/body/definition がそれぞれ支配する case を unit test する。
- asymmetric Above/Below、片側 definition、comparison なしを扱う。
- axis の上側だけにある definition band を非負 extent の和として過大評価しない。
- body/comparison/definitionのkind-pair gap policyを一つのtable/functionで管理する。
- comparison exclusionがrecord body reserveに包含されない入力を拒否する。
- 左definition列のようなX非交差pairを無視し、multi-record local headerのようなX交差する
  definition/body pairへdefinition clear gapを適用する。
- `canvas_band` を spacing helper の parameter にしない。

### Phase B2: Single-row と multi-record の統合

- one-record-per-row loop を pure clearance helper へ切り替える。
- `LinearRecordMeasurement`はenclosure用canvas extentとspacing用`CollisionBand`を別fieldで保持する。
- horizontal placementを先に確定し、multi-record solverはrowごとのbody/comparison band setと、
  実X範囲を持つ`local_band` / `row_band`のdefinition band setを同じhelperへ渡す。
- comparison corridor は、その境界を跨ぐ explicit/generated comparison がある場合だけ予約する。
- placement 後に translated `canvas_band` から content top/bottom と viewBox を求める。
- comparison endpoint は引き続き `comparison_exclusion_band` から求める。
- diagnostics metadata に各 band と、選択された constraint 値を出せるようにする。

### Phase C: 統合、全 session 再生成、文書

1. Workstream A/B の focused tests を通す。
2. browser wheel を current Python code から準備する。
3. Python/Webのcurrent writer constantを同じchangeで35へ一度だけ切り替える。
4. Gallery refresh mergeを更新し、refreshed側のversionとLOSAT artifactsを採用できるようにする。
5. 全11 Gallery session と全2 current test input session を version 35 へ更新する。
6. legacy fixtureはcurrent writer/refresh toolへ一度も通さない。
7. session diff で version、render schema、cache schema、manifest、設定値を review する。
8. Gallery source/example SVG の Y geometry diff を review する。
9. review 後に thumbnail と必要な tutorial media だけ更新する。
10. dedicated browser acceptance runnerでLoad → Save → Load → Generateと、Generate → Save → Reload →
    Generateを確認する。
11. browser acceptanceをrequired gateとして通し、skipを成功扱いしない。
12. non-slow suite、ruff、reference comparison を通す。

## 7. ファイル別の変更計画

### 7.1 LOSAT identity/cache

- `gbdraw/analysis/protein_colinearity.py`
  - canonical feature analysis ID、display alias、transport encoding、protein set、record analysis、
    manifest/hashのsingle owner。
  - schema 2 ID remap と Python/CLI raw cache schema 3。
- `gbdraw/web/js/app/python-helpers.js`
  - metadata hash ID の複製を削除し、Python owner を呼ぶだけにする。
- `gbdraw/web/js/app/run-analysis.js`
  - File fingerprint 由来の protein identity を廃止する。
  - schema 3 lookup、lazy migration、current manifestをorchestrationする。
  - LOSAT executionを小さいexecutor interface経由にし、acceptance testでworker callをcountできるようにする。
- `gbdraw/web/js/app/losat-cache.js`（新規、pure function のみ）
  - protein schema 3 / nucleotide schema 2のdiscriminated cache payload、validator、lookup、
    promotion、legacy candidate状態機械、copy-on-write mapping。
  - 現行の共有`LOSAT_CACHE_SCHEMA`をprotein/nucleotideの別ownerへ分け、program/flowだけに依存する
    暗黙判定を散在させない。
- `gbdraw/web/js/services/config.js`
  - raw/derived/manifest/legacy candidate envelope serialization、dual-read、session version 35。
  - `lastModified: 0` を保持する。
- `gbdraw/web/js/services/session-authority.js`
  - `legacyArtifacts`をartifact authorityとして宣言し、canonical request promotionから分離する。
- `gbdraw/web/js/state.js`、`gbdraw/web/js/app/app-setup.js`、reset owner
  - protein set、record analysis、record-instance binding、legacy candidateをcache stateとして
    初期化・snapshot・resetする。
- `gbdraw/web/js/services/gallery-session-migration.js`
  - render schema 2→3 promotion 中に legacy cache/evidence を落とさない。
  - protein ID の runtime 変換はここで行わない。
- `gbdraw/web/js/services/session-request.js`
  - canonical resource metadata を deterministic に materialize する。
- `gbdraw/linear.py`、`gbdraw/session_io.py`、必要に応じて `gbdraw/session.py`
  - CLI/Web の cache schema、manifest、session versionを一致させ、統合gateでversion 35へ切り替える。

### 7.2 Linear layout

- `gbdraw/layout/linear.py`
  - `CollisionBand`、X overlap、kind-pair clear-gap policy、`required_axis_gap()`のpure owner。
- `gbdraw/layout/linear_multi_record.py`
  - X範囲付きrow measurementとspacingを同じdomain contractへ変更する。
- `gbdraw/diagrams/linear/assemble.py`
  - definition geometry を canvas enclosure と row constraint に別々に渡す。
  - `canvas extent + comparison_height + canvas extent` を削除する。
- `gbdraw/diagrams/linear/track_slots.py`
  - `record_body_band`、`comparison_exclusion_band`、`canvas_band` の既存 semantics を維持する。
  - 新しい spacing 式を重複実装しない。

### 7.3 Sessions、tests、docs、generated assets

- `tools/refresh_gallery_sessions.py`
  - 全Gallery session inventoryとtransactional refreshを使用する。
  - `_merge_refreshed_gallery_artifacts()` はpromoted側のcanonical `renderRequest` / resource authorityを
    維持しつつ、refreshed側の`version`、`losatCache`、`losatDerivedCache`、新manifest/artifactを採用する。
  - staged validatorでversion 35、protein raw schema 3、derived schema 2、manifest schema 1を確認し、
    旧artifactへ戻ったsessionをcommit前に拒否する。
  - 旧orthogroup/editor overrideはfeature analysis identityで新artifactへ再適用する。
- `gbdraw/web/gallery/sessions/`
  - 11 session を version 35 へ再生成する。
- `tests/test_inputs/`
  - 2 current session fixture を version 35 へ再生成する。
- `tests/fixtures/sessions/`（新規候補）
  - 最小化した実在 schema v2 migration fixture を current inventory と分離して保持する。
  - fixtureと同じ場所にexpected metrics JSONを置き、Node/Python adapter共通のacceptance oracleにする。
- `tests/web/losat-cache.test.mjs`（新規）
  - pure cache validator、legacy状態遷移、manifest reference、Save-before-Generate payloadを検証する。
- `tests/web/losat-cache-migration.playwright.spec.js`（新規）
  - 実在fixture、counting LOSAT executor、structured telemetryを使うbrowser acceptance本体。
- `tests/run_losat_cache_browser_acceptance.py`（新規）
  - Node `@playwright/test`を優先し、なければPython Playwrightで同じfixture/expected metricsを実行する。
    両方利用不能、assertion未実行、browser起動失敗はexit 0やskipにせず失敗させる。
- `.github/workflows/test.yml`のbrowser test job
  - 上記runnerをrequired gateとして実行し、Chromiumとbrowser wheelを明示的に準備する。
- `docs/PYTHON_SESSION_COMPATIBILITY_MATRIX.md`
  - version 35、raw/derived/manifest schema、legacy policy を追記する。
- `docs/LINEAR_TRACK_OCCUPANCY_LAYOUT_IMPLEMENTATION_PLAN.md`
  - 実装後、constraint composition regression の修正結果を記録する。
- `docs/TUTORIALS/4_Protein_Comparisons.md`、`docs/TUTORIALS/7_Linear_Layout.md`、
  `docs/TUTORIALS/8_Interactive_SVG_Sessions.md`、`docs/FAQ.md`、release notes
  - user-visible ID、cache reuse、Auto spacing の説明が必要な箇所だけ更新する。
- Gallery source/example SVG、thumbnail、tutorial media
  - reviewed geometry change に限定して公式生成 path で更新する。

## 8. Test plan

### 8.1 LOSAT pure/Python tests

- source `protein_id` をdisplay aliasに採用するが、machine feature IDには使わない。
- missing/duplicate `protein_id`でも同じcanonical feature ID規則を使い、display aliasの重複は
  feature analysis ID suffixで一意化する。
- reserved delimiter、空白、control character、Unicodeをpercent-encodeし、transport tokenに
  非空白文字だけが残る。
- 同位置featureの永続ordinal、compound location、0-based/end-exclusive座標規則を固定する。
- 同じrecord source IDを複数配置してもcombined protein mapが衝突しない。
- record-local feature ID uniquenessとrecord-instance prefixがWeb/CLIで一致する。
- 同じprotein setでも別record-instance bindingのraw cache entryをdirect hitにしない。
- 異なるrecord analysisが同じprotein setを共有してもrecord固有metadataを上書きしない。
- compound location、同一 AA の別 feature、GFF3+FASTA を扱う。
- filename/mtime/resource rename で ID/protein-set hash が不変。
- `protein_id` / `locus_tag` / GFF `ID`の変更でfeature analysis IDとprotein-set hashは不変だが、
  transport bindingとraw keyが変わる。
- AA、visible CDS、genetic code、search args の変更で raw miss。
- product/gene/note、view reverse、filter option の変更で raw hit + derived miss。
- schema 2 TSV の全 query/subject を schema 3 ID に変換する。
- ambiguous/missing/corrupt legacy mapping を拒否する。
- query/subject reverse pairをdirect cache hitにしない。
- Python ownerのmanifest/hashとCLI/Web boundary serializationがgolden fixtureで一致する。

### 8.2 Web unit/session tests

- protein schema 3 raw entry、nucleotide schema 2 raw entry、legacy candidate envelope schema 1、
  schema 1/2 derived entryのdiscriminated validator。
- protein set/record analysis/record-instance bindingのdeduplication、round-trip、size/pruning。
- mixed protein/nucleotide cache を誤分類しない。
- schema 2→3 render promotion が legacy cache/evidence を保持する。
- Load → Save → LoadをGenerate前に行っても`pending` candidateがbyte-identicalに残り、current
  `losatCache`へprotein schema 2が混入しない。
- migration成功時だけ`pending`から`promoted`へ進み、失敗時は理由付き`rejected`を保持する。
- `lastModified: 0` が save/load 後も 0。
- version 27–34 reader と version 35 writer。
- current session inventory は legacy fixture を除き全件 version 35。
- Gallery session file list と refresh tool/test inventory が一致する。
- Gallery refresh後のversion/cache/manifestがrefreshed artifact側から採用され、旧schemaへ戻らない。
- current v35 protein artifacts内に未解決の旧`p_r_...`参照がない。
- derived、orthogroup member/edge/path、selection、editor overrideの全protein IDがmanifest/bindingへ
  一意に解決し、名前・説明・選択の意味が保持される。

### 8.3 Layout unit/integration tests

- body、comparison、definition の各制約が支配する case。
- axis gap が三制約の `max()` と一致し、合計値にならない。
- Xが交差しないdefinition/body bandはspacingを増やさず、Xが交差するdefinition/bodyおよび
  definition/definition pairの最大値だけを採用する。
- multi-record local headerと隣接row bodyがX方向に重なるcaseでdefinition clear gapを保証する。
- comparison exclusion bandがrecord body reserveに包含されない入力を拒否する。
- large definition font/line count は必要な時だけ definition constraint を増やす。
- comparison なしでは 60 px を予約しない。
- 一部の row 境界だけに comparison がある場合、comparison のない境界では 60 px を予約しない。
- axis 上側だけにある local header と row definition を signed band edge で処理する。
- dense labels/feature lanes/depth/annotation に comparison が侵入しない。
- default/custom slots、Middle/Above/Below、ribbon/curve。
- one-record-per-row と multi-record-row。
- definition on/off、片側 definition、非対称 track stack。
- canvas が全 paint を含み clip しない。
- X/bp mapping、feature ID、色、record order は不変。

### 8.4 Browser acceptance

fixtureと同じ場所のexpected metrics JSONをbrowser acceptance contractのsingle ownerとする。
`tests/web/losat-cache-migration.playwright.spec.js`はNode adapterとして、実在schema v2 fixture、
counting LOSAT executor、structured telemetryを使い、console文字列の解析に依存せず次を確認する。

```text
totalPairs=25
cacheHits=25
cacheMisses=0
uniqueJobs=0
LOSAT worker calls=0
```

さらに次を確認する。

- filename/mtime だけを変えても 0 job。
- 一つの AA/protein set を変えると影響する pair だけ job が走る。
- Load → Save → Load → Generateでもlegacy candidateを失わず0 job。
- Generate → Save → Reload → Generateを二回行っても0 job。
- QUERY/SUBJECT は `p_r_<metadata-hash>...` ではなく readable stable ID。
- 全 QUERY/SUBJECT が manifest から一意な feature へ戻る。
- derived/orthogroup/editor/result内の全protein referenceが同じmanifest/bindingへ解決する。
- orthogroup名・説明overrideと選択中featureがsave/load後も同じ対象を指す。
- BGC layout の axis gap は band の `max()` 式と一致する。
- definition bbox、record body、comparison corridor が各 domain 内で交差しない。
- preview と exported SVG の geometry が一致する。

`tests/run_losat_cache_browser_acceptance.py`はNode `@playwright/test`を検出できれば上記specを実行し、
なければPython Playwright adapterで同じfixtureとexpected metrics JSONを検査する。両adapterは
browser操作だけを所有し、期待値を複製しない。両方利用不能、Chromium
起動失敗、assertionが一件も実行されない場合はnon-zeroで終了する。CIではこのrunnerをrequired gate
とし、skipを成功扱いしない。

### 8.5 Verification commands

実装時は変更範囲に合わせ、最低限次を実行する。

```bash
node --check gbdraw/web/js/app/run-analysis.js
node --check gbdraw/web/js/app/losat-cache.js
node --check gbdraw/web/js/services/config.js
node tests/web/losat-cache.test.mjs

pytest tests/test_protein_colinearity.py -v
pytest tests/test_linear_vertical_layout.py -v
pytest tests/test_linear_multi_record_layout.py -v
pytest tests/test_linear_definition_alignment.py -v
pytest tests/test_linear_track_layout.py -v
pytest tests/test_session_io.py tests/test_api_session.py -v
pytest tests/test_web_packaging.py -v

python tests/run_losat_cache_browser_acceptance.py

pytest tests/test_output_comparison.py::TestOutputComparison -v
pytest tests/ -v -m "not slow"
ruff check gbdraw/
```

browser acceptance runnerはNode Playwrightを優先し、利用できなければPython Playwrightで同じ
assertionを行う。runner自体のskipは許可しない。reference outputは通常testで先に差分を確認し、
geometry changeをreviewした後だけ`--update-reference-outputs`で更新する。

session と Gallery asset の最終更新は transactional tool で行う。

```bash
# schema/cache/session差分だけを先に確認
python tools/refresh_gallery_sessions.py --no-assets

# layoutを含む最終承認後、sessionとvisible assetsをまとめて更新
python tools/refresh_gallery_sessions.py
```

## 9. Rollout order

1. Phase 0 fixture と failing tests を review する。
2. Workstream A1/A2を実装し、protein set、record analysis、binding、transport IDのcontractを成立させる。
3. schema v2 BGC fixtureのLoad → Save → Load → Generateと25/25 verified reuseを成立させる。
4. Workstream B1/B2 を実装し、single/multi-row の constraint solver を統一する。
5. Python/Web session version を 35 へ一度だけ切り替える。
6. Gallery merge toolを新artifact authorityに対応させる。
7. 全 current session を一括再生成し、legacy fixture だけを旧session/protein形式で残す。
8. SVG/reference/Gallery visual diff を review する。
9. docs、release notes、必要な media を更新する。
10. full verification後にmetadata-dependent ID generator、legacy current-cache write path、旧spacing式を削除する。

protein schema 2 reader/migratorはversion 35 releaseでは残す。protein schema 2 entryをcurrent cacheへ
書くwriterとmetadata-dependent ID generatorは残さない。import由来candidateのlossless legacy envelope
writerとnucleotide schema 2 writer/readerは維持する。rollback時もlegacy protein schema 2 entryを
破壊しないよう、lazy migrationはcopy-on-successとする。

## 10. Risks and mitigations

| Risk | Mitigation |
|---|---|
| legacy cache の false hit | full FASTA hash、AA digest、program/args、一対一 mapping をすべて検証する |
| 同じprotein setを持つ異なるrecordのmetadata上書き | `ProteinSet`からrecord固有情報を外し、record analysis/bindingで参照する |
| missing/duplicate source ID | source IDをdisplay aliasだけに使い、canonical feature analysis IDで一意化する |
| transport IDのescape衝突 | normative NFC/UTF-8 percent encodingとreserved-character/property test |
| Generate前の再保存でlegacy cacheを失う | current cacheと分離したlegacy candidate envelopeをlosslessにround-tripする |
| manifest による session 肥大化 | protein setをhashでdeduplicateし、instance bindingは1 rowにつき1件、pair entryはreferenceのみ保持 |
| empty legacy result を誤って再利用 | record token を証明できなければ pair-local miss |
| Web/Python ID drift | Pythonをsingle ownerにし、Web boundary serializationのbyte-identical golden test |
| derived payload に古い record index が残る | schema 1をblind reuseせずcurrent manifestからschema 2を再構築 |
| layout を詰め過ぎてtrack collisionが戻る | body/comparison/definition の独立 non-overlap test |
| canvas clipping | spacingとは別にtranslated canvas bandsからviewBoxを計算 |
| multi-recordだけ別式が残る | single/multiで同じ pure clearance helperを使用 |
| definition bandの水平domain誤判定 | 実X区間とkind-pair policyを使い、local header/bodyのcross-kind testを固定 |
| 広範なSVG差分 | X mappingとstyleを固定し、Y/viewBox/comparison pathの差だけをreview |
| session一括更新の途中失敗 | transactional refreshと更新後全件validation |
| stale tutorial screenshot | visible旧spacingを含むmediaだけ実UIから再撮影 |
| browser acceptanceがskipされる | Node/Python fallback runnerをrequired gateにし、未実行をnon-zeroにする |

## 11. Non-goals

- LOSAT/BLAST の検索アルゴリズム、score、orthogroup inference semantics を変更しない。
- 外部 database で `protein_id` を再照合しない。
- 全ファイル形式を横断する普遍的 biological deduplication system を作らない。
- record-instance間のraw TSV neutralization/rebinding最適化を初回実装に含めない。
- nucleotide LOSAT cache を protein manifest へ無理に統合しない。
- canonical `renderRequest` schema 3 全体を再設計しない。
- Circular layout を同時に変更しない。
- 汎用2D collision solver、DOM bbox反復solver、scene graphを導入しない。
- 旧105 pxを固定値として復元しない。
- 検証不能な旧cacheをhit扱いしない。

## 12. Acceptance criteria / Definition of Done

次をすべて満たした時だけ完了とする。

1. 実在 BGC schema v2 fixture が `25 hits / 0 misses / 0 jobs` で Generate できる。
2. 検証不能または破損した legacy pair は誤利用せず、その pair だけ再実行する。
3. feature analysis IDはtype、全location parts、strand、永続同位置ordinalから一意に決まり、source
   qualifierの一意性に依存しない。
4. new QUERY/SUBJECTは規範的percent encodingに従うreadable transport IDで、全件をmanifestから
   一意に実record/featureへ解決できる。
5. 同じsource recordを複数配置してもtransport IDとcombined protein mapが衝突せず、別bindingの
   raw TSVをdirect hitとして誤利用しない。
6. 異なるrecord analysisが同じprotein setを共有してもrecord固有metadataとbindingを上書きしない。
7. filename、mtime、resource rename、save/loadだけではraw cache keyが変わらない。
8. AA、protein membership、意味のあるLOSAT argsの変更では該当cacheが確実にmissする。
9. source/display qualifierの変更はfeature analysis IDを変えないが、transport bindingを変えて
   raw missとderived rebuildになる。
10. reverse query/subject pairを未検証のdirect hitとして再利用しない。
11. raw cacheとderived cacheのinvalidation boundaryがunit testで固定される。
12. Pythonのsingle ownerがID/manifest/hashを生成し、CLIとWeb boundaryがbyte-identicalにserializeする。
13. current schema 3 artifact内の全protein referenceがmanifest/bindingへ解決し、旧`p_r_...`参照が
    残らない。import由来の未検証entryだけはlegacy envelope内に隔離される。
14. Load → Save → LoadをGenerate前に行ってもlegacy candidateがlosslessに残り、current cacheへ
    protein schema 2 entryが混入しない。
15. orthogroup名・説明overrideと選択状態がID移行後も同じfeatureを指す。
16. axis gapはXが交差するeligible body/comparison/definition kind pairの必要間隔の`max()`であり、
    二重加算しない。
17. `comparison_height=60`はexclusion edge間の最小clear corridorとして保証される。
18. comparisonのないrow境界は`comparison_height`を予約しない。
19. Xが交差しないdefinition/body pairはspacingを増やさず、Xが交差するlocal header/body pairは
    definition clear gapを満たす。
20. definition、non-overlay body、comparisonが各collision domainで交差しない。
21. off-axis definitionを非負extentの和で過大評価しない。
22. canvasが全painted contentを含み、clipしない。
23. single/multi-row、default/custom、ribbon/curve、sparse depthの回帰testが通る。
24. Python/Webのcurrent session versionは35で一致し、通常の同梱session全13件がversion 35になる。
25. repository内の旧session envelope/render schemaとprotein raw schema 2は専用migration fixtureだけに残る。
26. current sessionではprotein raw schema 3とnucleotide raw schema 2の混在を正しく検証できる。
27. Gallery refresh後もversion/cache/manifestが旧値へ戻らない。
28. 全Gallery session、SVG、thumbnail、必要なtutorial mediaがreview済みcurrent出力と一致する。
29. dedicated browser acceptanceがskipなしでrequired gateを通り、focused tests、non-slow suite、ruff、
    reference comparisonが成功する。

## 13. Implementation progress handoff（2026-07-21）

### 完了した実装

- Python single ownerにstable feature/record/protein-set identity、readable transport ID、schema-1 manifest、
  directional schema-3 protein raw keyを実装した。WebはPyodide helper経由で同じownerを使う。
- session version 35、`renderRequest.schema == 3`、protein raw schema 3、nucleotide raw schema 2、
  derived schema 2、manifest schema 1の境界をPython/Webに実装した。readerはversion 27～35を受理する。
- legacy protein raw/derived artifactをcurrent cacheから隔離し、Save-before-Generateで保持される
  schema-1 envelopeと、検証済みpairだけをcopy-on-writeでpromotionするpathを実装した。
- X範囲付き`CollisionBand`、kind-pair policy、`required_axis_gap()`を実装し、single/multi-rowを
  同じsolverへ統合した。body/comparison/definitionはeligible制約の`max()`で合成し、
  comparisonのない境界にcorridorを予約しない。

### 確認済みの範囲

- protein identity/cache、legacy fixture、session codec/API、Linear vertical layout、Web pure cache/session の
  focused testおよびJavaScript syntax checkで実装単位の動作を確認した。
- 統合focused suiteは`190 passed, 1 skipped`で完了した。repository-wide test suiteではなく、
  protein/session/API/refresh/Linear vertical layoutに限定した結果である。
- 履歴上のBGC0000708–BGC0000713 schema-v2 sessionをgzip fixtureとして固定し、25 pairの
  legacy candidateがPython側の検証対象になることを確認した。
- 同BGC sessionのcanonical Python replayとstaged validatorは、version 35、protein raw schema 3が
  25件、derived schema 2が1件、manifest schema 1、legacy candidateなし、旧`p_r_...`参照なしで通過した。
- browser acceptance runner/spec/CI gateを追加し、fixture SHA/version/schema/count、test discovery、
  runner syntaxを確認した。通常sandboxではlocal serverの`listen`が拒否されるため、browser実行には
  sandbox escalationが必要である。

### 未完了・未確認のgate

- 実browserでのLoad → Save → Load → Generateと`25 hits / 0 misses / 0 jobs`のrequired acceptance。
  escalated Node runは実在fixtureの1 testを開始し、Save-before-Generateで見つかった全null Depth行の
  round-trip不具合を修正した後も同地点を通過したが、最終assertion完了前に時間制限で中断した。
- current Gallery/test sessionのtransactional refresh、session/SVG/thumbnail差分のreview、必要なmediaの判定。
- 意図したLinear Y/viewBox reference差分の更新と再比較、non-slow suite、repository-wide ruff、build。

### 次回の実行順

1. browser wheelはcurrent Python codeから準備済みである。sandbox escalation付きで
   `python tests/run_losat_cache_browser_acceptance.py`を完走し、続けて
   `python tests/run_losat_cache_browser_acceptance.py --python`でfallback adapterも確認する。
2. `python tools/refresh_gallery_sessions.py --no-assets`でartifact差分を確認した後、
   `python tools/refresh_gallery_sessions.py`でsessionとvisible assetをtransactionalに更新する。
3. output comparisonで差分を確認し、意図したLinear referenceだけを公式コマンドで更新する。
4. focused/full validationと最終diff reviewを完了するまで、Definition of Done全体を完了と扱わない。
