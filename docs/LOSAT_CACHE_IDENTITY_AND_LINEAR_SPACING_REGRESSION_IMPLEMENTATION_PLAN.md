# LOSAT Runtime Handle Cache Identity and Linear Spacing Regression Implementation Plan

- 作成日: 2026-07-21
- 更新日: 2026-07-25
- 対象: Web/CLI Linear の LOSATP cache identity、session version 27–35 読み込み、Linear record spacing
- 現行基準: session version 36、`renderRequest.schema == 3`、raw protein LOSATP cache schema 4、
  derived cache schema 3、protein identity manifest schema 2
- 非対象cache: nucleotide LOSAT は raw schema 2 を維持し、session内のprotein raw 4 /
  nucleotide raw 2混在を正式に扱う
- 状態: compact runtime handle設計への更新、実装、最終検証完了
- 関連設計: [`LINEAR_TRACK_OCCUPANCY_LAYOUT_IMPLEMENTATION_PLAN.md`](LINEAR_TRACK_OCCUPANCY_LAYOUT_IMPLEMENTATION_PLAN.md)、[`PYTHON_SESSION_CANONICAL_REQUEST_PLAN.md`](PYTHON_SESSION_CANONICAL_REQUEST_PLAN.md)
- 移行設計: [`LOSAT_COMPACT_RUNTIME_HANDLE_MIGRATION_IMPLEMENTATION_PLAN.md`](LOSAT_COMPACT_RUNTIME_HANDLE_MIGRATION_IMPLEMENTATION_PLAN.md)

本書の旧version 35設計にあった「readable long transport IDを内部raw TSVへ反復保存する」方針は
廃止する。version 35 / raw schema 3 / derived schema 2 / manifest schema 1は移行入力と履歴として
のみ扱い、current artifactのcontractはversion 36 / raw 4 / derived 3 / manifest 2とする。
Linear spacingのcollision-domain設計はこの更新では変更しない。

## 0. 結論

今回確認された現象は、原因の異なる二つの回帰である。

1. **Workstream A — LOSAT identity/cache**
   - Web LOSATP の QUERY/SUBJECT ID が、ファイルの `name`、`size`、`lastModified` から作る
     record hash に依存している。
   - schema v2 session の復元で `File` metadata が変わると、配列が同じでも FASTA header と
     FASTA 全文 hash が変わり、保存済み raw cache を再利用できない。
   - version 35の長いreadable transport IDはmetadata依存を解消したが、raw/derived payloadへ
     64桁feature hashを反復し、実データではsession size gateを超えた。
   - `featureAnalysisId`をmanifest authority、短い`runtimeHandle`をFASTA/raw/derivedの内部参照、
     `exportId`をdownload時にhydrationするuser-visible IDとして分離する。
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
- size上限を満たすためにvalidなraw/derived entryをpruneする。
- 内部runtime handleを未解決のままユーザー向けTSVへ出力する。
- 旧レイアウトを長期 flag として併存させる。

## 1. 再現結果と固定する baseline

対象は、履歴上の実在する BGC0000708–BGC0000713 schema v2 session とする。version 35の
schema-3 sessionをschema 2と書き換えただけのsynthetic fixtureは、migrationの証拠にしない。

| 観測項目 | 保存済み schema v2 | 旧コードで Load → Generate | 修正後の条件 |
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

compact identityのsize regressionは、59 raw entryを持つVibrio Gallery sessionを追加baselineとする。
entryを欠落させず、gzip `< 100,000,000 bytes`、expanded JSON `< 536,870,912 bytes`をhard gateとする。

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

`p_r_...` は cache key そのものではなく、LOSAT に渡す旧内部 protein ID である。実際の
`protein_id`、`locus_tag`、`featureSvgId`、location は別の `protein_map` にしか存在しない。
raw TSV だけでは、この内部 ID を実 feature へ戻せない。

schema v2 promotion/materialization は埋め込み resource から `File` を再構築する。元 metadata を
完全には保存できず、特に `lastModified: 0` を `entry.lastModified || Date.now()` で現在時刻へ
置換する path もある。このため biological content が同じでも `recordInstanceKey` が変わる。

version 35では、安定したfeature identityを次のreadable transport IDへ含めることでこのmissを
解消した。

```text
<record-source-id>@<record-instance-id>|<display-alias>~<feature-analysis-id>
```

しかし、この長いIDをraw TSVの全hit行とderived payloadへ反復したため、Vibrio sessionでは
QUERY/SUBJECT列だけで約136 MBとなった。version 36では長いIDをcurrent artifactへ保存せず、
manifest schema 2で一度だけ保持するmachine/display metadataを短いsession-global handleから
参照する。

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
| Runtime handle | FASTA、内部raw TSV、`protein_map`、derived payloadの短い参照 | record-instance identity、feature analysis identity |
| Export ID | `Save Raw LOSAT TSV`で人に見せるID | display alias、同一record instance内のduplicate ordinal |
| Runtime binding hash | Raw LOSAT実行結果のbinding identity | protein-set identity、record-instance identity、runtime-handle map |
| Display binding hash | export/UI/derived metadataのinvalidation | record/display metadata、export mapping |
| Raw cache fingerprint | LOSAT 実行結果の再利用可否 | ordered protein-set identity、runtime binding、program、outfmt、意味のある LOSAT args |
| Derived cache fingerprint | genomic link/orthogroup/描画 payload の再利用可否 | raw key、runtime/display mapping、view transform、filter/converter settings |
| UI pair identity | editor と SVG の comparison 行を結ぶ | canonical record key と pair direction |

次を biological identity、runtime handle、runtime binding hash、raw cache fingerprint に含めない。

- upload filename
- resource filename/resource ID
- file size を単独の identity とした値
- `lastModified`
- session を開いた時刻
- temporary Pyodide path
- 永続化されない現在の配列 index だけから作る token

### 3.2 Runtime handle と QUERY/SUBJECT contract

新規LOSATP jobの内部QUERY/SUBJECTは、machine identityとdisplay identityを混ぜず、次の二段階で
作る。qualifierが一意かどうかによってmachine identityの作り方を変えない。

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
2. `runtimeHandle`は次のdomain-separated payloadからPythonのpure functionだけで作る。

   ```json
   {
     "featureAnalysisId": "f_...",
     "recordInstanceKey": "record-1"
   }
   ```

   key辞書順・空白なしのUTF-8 canonical JSONへ
   `b"gbdraw-runtime-handle-v1\0"`を前置してSHA-256を計算し、先頭16 bytesをRFC 4648 Base32の
   lowercase、paddingなしで表し、`h_`を付ける。

   ```text
   h_[a-z2-7]{26}
   ```

同じaccession/record IDを複数ロードした場合は、File metadataではなくsessionに保存されたcanonical
record keyを使うためhandleはinstanceごとに異なる。同じinstance/featureからは常に同じhandleを作り、
alias、product、gene、note、filename、mtimeでは変えない。manifest validatorは全record instanceを
通じた重複をfail closedにし、衝突時の非決定的suffixを作らない。

内部raw TSVの第1・第2列、LOSATP FASTA header、combined `protein_map`、derived payload内の
protein referenceはruntime handleへ統一する。完全なfeature identityとdisplay metadataのauthorityは
manifest schema 2とし、raw/derivedの各hit/member/edgeへ反復保存しない。

`Save Raw LOSAT TSV`ではdownload直前にmanifest resolver/hydratorを通し、第1・第2列だけを
`exportId`へ置換する。`exportId`は空でない`protein_id`、`locus_tag`、GFF `ID`、完全location
fallbackの順で選ぶ。同一record instance内でaliasが重複する場合だけfeature analysis ID順の
`~1`、`~2`を付ける。aliasはUnicode NFC化後にpercent encodeする。

hydratorは非comment行の厳密な12列、row順、3～12列、数値表現、改行を保存する。一件でも
unknown/wrong-binding/重複解決handleがあればdownload全体を失敗させ、内部handleをfallback出力
しない。50 MiB確認は内部cache文字数ではなくhydration後のUTF-8 byte数で判定する。ユーザーが
Pairwise Comparisonsへ投入したTSVは変更しない。

### 3.3 Raw cache と derived cache の invalidation

| 変更 | Runtime handle | Raw LOSAT | Derived payload | Export ID |
|---|---|---|---|---|
| filename / `lastModified` / resource rename | same | reuse | reuse | same |
| save → loadで同じbiological input | same | reuse | reuse | same |
| `protein_id` / `locus_tag` / GFF `ID`のみ | same | reuse | display metadataを含む時rebuild | change |
| product/gene/noteのみ | same | reuse | embedded時rebuild | usually same |
| record sourceの表示labelのみ | same | reuse | embedded時rebuild | UI/filenameのみchange |
| display reverse-complementのみ | same | reuse | view transformをrebuild | same |
| AA配列変更 | feature IDがsameの場合あり | protein-set hashでmiss | miss | same/updated |
| selector/region/visibilityによるmembership変更 | affected handles/set change | miss | miss | change |
| feature location/strand/ordinal変更 | change | miss | miss | change |
| record instance key変更 | change | miss | miss | duplicate scopeを再計算 |
| LOSAT program/outfmt/search args変更 | same | miss | miss | same |
| post-filter/orthogroup/collinearity optionのみ | same | reuse | miss | same |
| query/subject direction変更 | same handles | separately cachedでなければdirectional miss | miss | columns reverse |

raw reuse は「設定名が同じ」ではなく、canonical protein set と search semantics が同じことから
判断する。display binding hashをraw keyへ含めず、runtime binding hashとdisplay binding hashを
分離する。false hit より明示的な miss を優先する。

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

### 4.1 Protein set、record analysis、runtime/display binding

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
      "schema": 2,
      "recordAnalysisId": "sha256:...",
      "runtimeBindingHash": "sha256:...",
      "displayBindingHash": "sha256:...",
      "runtimeIds": {
        "f_2f...": "h_..."
      },
      "featureMetadata": {
        "f_2f...": {
          "displayAlias": "CAG38695.1",
          "featureSvgId": "fed46a3a6"
        }
      }
    }
  }
}
```

`proteinSetHash` は`featureAnalysisId`でsortしたcanonical protein listの`featureAnalysisId`と完全な
AA SHA-256から作る。`recordAnalysisId`はrecord source identity、selector/region、
`proteinSetHash`のcanonical JSONから作る。いずれもFASTA description、File metadata、record
instance key、JSON object insertion orderに依存させない。

protein-set manifestはpairごとに複製せずhashでdeduplicateする。同じprotein setを持つ異なる
record analysisも同じcontent objectを参照でき、record固有情報を上書きしない。同じrecordを二回
配置した場合は、一つのrecord analysisへ二つのrecord-instance bindingを作る。

runtime binding hashはcanonical `recordInstanceKey`、protein-set hash、全`runtimeIds` mapから作り、
record source label、display alias、feature display metadataを含めない。display binding hashは
`recordAnalysisId`、record source/display metadata、export mappingを別にcanonicalizeする。raw keyは
query/subjectのprotein-set hashとruntime binding hashを含み、保存TSV内の全handleがそのbindingに
属することを検証する。derived keyはraw keyに加えてruntime/display mappingを含む。これにより
aliasだけの変更でLOSATを再実行せず、UI/export metadataは必要に応じて再構築できる。

### 4.2 Schema version

変更後の version owner は次のようにする。

| Owner | v35 / migration input | v36 / current | 方針 |
|---|---:|---:|---|
| Session envelope | 35 | 36 | Python/Web writerを同時更新 |
| Canonical `renderRequest` | 3 | 3 | render request semantics は変更しない |
| Raw protein LOSATP cache | 3 | 4 | compact runtime handleを使用 |
| Raw nucleotide LOSAT cache | 2 | 2 | 現行contractを維持し、protein entryと混在可能にする |
| Derived LOSATP cache | 2 | 3 | protein referenceをruntime handleへ統一 |
| Protein identity manifest | 1 | 2 | runtime/display bindingを分離 |
| v27–34 protein raw candidate envelope | 1 | 1 | schema-2 verified lazy migrationを維持 |
| v35 protein raw candidate envelope | — | 1 | schema-3＋manifest-1をlosslessに隔離 |

current writerはversion 36だけを書く。protein LOSATP writerはschema 4、derived writerはschema 3、
manifest writerはschema 2だけを書き、nucleotide LOSAT writerはschema 2を維持する。readerは
version 27–34のprotein raw schema 2を`legacyArtifacts.proteinRawCandidates`へ、version 35の
protein raw schema 3＋manifest schema 1を`legacyArtifacts.proteinRawV35Candidates`へ隔離する。
どちらもGenerate前のSave → Loadでlosslessにround-tripする。current `losatCache`へ旧protein
schemaを混在させず、unknown/newer schemaを推測で読み替えない。

### 4.3 全 session file の更新

同梱current sessionは一部だけでなく全件をversion 36へ再生成する。

| Scope | 件数 | 更新方法 |
|---|---:|---|
| `gbdraw/web/gallery/sessions/*.gbdraw-session.json[.gz]` | 11 | `tools/refresh_gallery_sessions.py` の canonical path |
| `tests/test_inputs/*.gbdraw-session.json` | 2 | 各 fixture の canonical save/generator path |
| 固定 legacy migration fixture | 2 | schema v2 / version 35のまま保持し、current inventoryから除外 |

更新後は次を満たす。

- Gallery と通常の test input session はすべて `version == CURRENT_SESSION_VERSION == 36`。
- Linear session は現行の spacing、track slot、definition、comparison 設定を欠落なく保持する。
- LOSATP sessionはprotein raw schema 4、derived schema 3、manifest schema 2を保持する。
- canonical pathで再生成したcurrent sessionは旧protein raw/derived candidateを持たない。
- nucleotide LOSAT/BLAST cacheはraw schema 2を維持し、protein manifest migrationの対象にしない。
- `.json` と `.json.gz` は同じ validator を通す。
- file list の重複定義と実ファイルの差を inventory test で検出する。
- Vibrio sessionは59 raw entryを欠落なく保持し、gzip/expanded hard gateを満たす。

layout 変更も生成 SVG に影響するため、最終更新は session JSON だけで終えない。review 済みの
session から Gallery source SVG、example SVG、thumbnail、`examples.json` を公式 tool で更新する。
チュートリアル screenshot は実際に旧 spacing が写っているものだけを、Gallery media skill の
手順で再撮影する。

外部ユーザーが保有する旧sessionは一括更新できないため、version 36 readerのverified lazy
migrationで開く。repository内の旧session envelope、protein raw schema 2/3、derived schema 1/2、
manifest schema 1は専用migration fixture以外に残さない。current v36 session内のnucleotide raw
schema 2はこの制約の例外とする。外部旧sessionをGenerate前に再保存した場合だけ、旧artifactが
current cacheではなく明示的なlegacy candidate envelope内に残ることを許可する。

## 5. Legacy protein artifact の verified lazy migration

import 時点では Pyodide protein extraction が完了していないため、ID変換は import 中ではなく、
最初のGenerateでmanifest schema 2を作った後に行う。version 27–34のschema-2 entryとversion 35の
schema-3 entryをcurrent cache mapへ混在させず、それぞれ専用のread-only candidateとして隔離する。

candidateは次の小さい状態機械で扱う。

| State | 意味 | Save時の扱い |
|---|---|---|
| `pending` | 未検証。Generate前または検証待ち | original entryをlegacy envelopeへlosslessに保存 |
| `promoted` | schema 4 copyの作成とvalidationが完了 | schema 4 entryを保存し、legacy candidateは次snapshotから除去 |
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

Load直後にGenerateせずSaveしても`pending` candidateと、version 35ではsource manifest-1/derived
evidenceを失ってはならない。SaveのためだけにPyodide抽出を起動せず、Saveを禁止もしない。current
validator/serializerはschema-4 protein raw、schema-3 derived、manifest-2、schema-2 nucleotide raw
だけを扱い、旧protein artifactの解釈はmigration moduleだけが所有する。

### 5.1 Version 27–34 / protein raw schema 2

1. schema-2 entryをimport時に破棄せず、`proteinRawCandidates`へ保持する。
2. 現在のrecordからmanifest-2とruntime handle mapを抽出する。
3. legacy `p_r_` QUERY/SUBJECTのrecord token、location、strand、AA digestを解析する。
4. derived/orthogroup metadataのrecord index、source ID、feature SVG ID、locationを補助evidenceに使う。
5. legacy IDと現在のfeatureが一対一対応し、旧FASTA hash、AA digest、program、outfmt、args、
   directionが一致することを確認する。
6. TSV第1・第2列だけをmanifest-2のruntime handleへ変換し、row順と第3～12列を保持する。
7. schema-4 keyを再計算し、copy-on-successでcurrent cacheへ昇格する。
8. derived schema 1はdirect hitにせず、promoted rawとcurrent manifestからderived schema 3を再構築する。

### 5.2 Version 35 / raw schema 3 / manifest schema 1

1. version 35 validatorでmanifest-1、schema-3 key/binding、TSV全行を検証し、
   `proteinRawV35Candidates`へsource manifestとともにlosslessに隔離する。
2. 旧long transport IDをmanifest-1から`featureAnalysisId`へ逆引きする。
3. current manifest-2の`runtimeIds`へ置換する。
4. protein set、AA digest、record instance、direction、program、outfmt、argsを再検証する。
5. row順と第3～12列をbyte-equivalentに保ってschema 4 keyを再計算する。
6. alias/display metadataだけが変わっていてもruntime identityが一致すれば昇格を許可する。
7. derived schema 2はdirect hitにせず、raw schema 4からderived schema 3を再構築する。
8. orthogroup state、editor override、selectionはartifact inventory resolverでruntime handleへ移す。

両経路ともschema-4 rawとderived schema 3のvalidationに成功した後だけcandidateを除去する。全検証が
終わるまで旧entryを削除・上書きしない。genericな文字列置換やschema relabelでcurrent artifactを
作らない。

次の場合は、その pair だけ cache miss として LOSAT を再実行し、理由を UI/console に残す。

- legacy ID が欠損または解析不能。
- empty result で record token を他 entry/metadata から一意に復元できない。
- compound location または重複 feature が曖昧。
- AA digest、FASTA hash、program、args のいずれかが不一致。
- manifest が同じ legacy ID を複数 feature へ解決する。
- session cache が破損している。

BGC schema v2 fixture は保存済み derived/orthogroup metadata が十分であり、25 pair 全件を migration
できることを acceptance test にする。全 legacy session の無条件 hit は保証しない。

version 36 readerで`pending`または`rejected` candidateを再読込できること、Load → Save → Load →
Generateでも元のcandidateと25/25 reuseを維持することをsession contractに含める。version 35
schema-3 sessionもLOSAT worker 0件でschema 4へ移行する。repository内でcanonical生成するcurrent
sessionにはlegacy candidateを残さない。

`lastModified: 0` は有効値として `??` または明示判定で復元し、`Date.now()` に置換しない。ただし
この修正は migration の補助であり、新 cache identity 自体は `lastModified` 非依存とする。

## 6. 実装 phase

### Phase 0: Characterization と failing regression tests

- 履歴上の実在 BGC schema v2 session を最小化・gzip 化し、専用 legacy fixture に固定する。
- fixture metadata/testで `source commit=c64ff8c`、`session.version=33`、
  `renderRequest.schema=2`、protein raw schema 2、derived schema 1を固定する。
- 保存raw entryは34件、対象Generate pairは25件であることを別々にassertする。
- version 35 fixtureはraw schema 3、derived schema 2、manifest schema 1、保存entry数を固定し、
  current inventoryから除外する。
- version 35 raw schema 3をschema 2とrelabelするsynthetic testと区別する。
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
- canonical location、同位置ordinal、feature analysis ID、runtime handle、export ID mappingを一つの
  pure Python moduleで実装する。source qualifierの重複有無でmachine IDの分岐を作らない。
- `ProteinSet`をrecord非依存のcontent objectにし、`RecordAnalysis`と`RecordInstanceBinding`へ
  record固有情報とsession配置情報を分離する。
- Web/CLIとも全input rowのrecord-instance keyを先に確定し、それをshared extractorへ渡す。
  Webのper-record extractionとCLIのbatch extractionでuniqueness scopeを変えない。
- Webは既存row `uid` / canonical `renderRequest.records[].recordKey`、Pythonは既存
  `gbdraw_record_key` / canonical record keyを再利用し、別のrandom ID体系を追加しない。
- identical recordを複数配置したfixtureでruntime handleは別、protein-set hashは同じ、raw cache
  keyはbindingごとに別になることを確認する。
- 異なるrecord analysisが同じprotein setを持つfixtureで、一つの`ProteinSet`を共有しつつ
  `recordAnalysisId`とbindingを上書きしないことを確認する。
- runtime binding hashからdisplay metadataを除き、display binding hashを別ownerにする。
- Python ownerのfeature/record/protein-set identity、runtime/display binding、runtime handle、
  raw key、manifest、hydration出力をWebがbyte-identicalにserializeするgolden testを追加する。
  Webに第二のidentity/raw-key実装を作らない。derived keyはsurface-localとし、両surfaceで同じ
  意味のinvalidation inputsをすべて含める。
- 128-bit handle golden vector、duplicate handle rejection、Unicode/NFC、duplicate alias、
  compound location、同位置featureを含むproperty/golden testを追加する。
- filename、mtime、resource name を変えても manifest と protein-set hash が変わらないことを確認する。

### Phase A2: Raw/derived cache schema と lazy migration

- cache key construction/lookup/promotion を `run-analysis.js` から focused module へ移す。
- protein raw schema 4、nucleotide raw schema 2、derived schema 3、manifest schema 2の
  discriminated validatorとserializerを実装する。size対策のentry pruningは実装しない。
- protein schema 2/3はmigration candidateとしてdual-readし、新規書き込みはschema 4に限定する。
  nucleotideの
  新規書き込みはschema 2を維持する。
- legacy protein schema 2 entryとv35 schema-3 entryはcurrent cache mapへ入れず、別のlegacy
  candidate envelopeで`pending/promoted/rejected`を管理する。
- schema 4 raw keyにquery/subjectのruntime binding hashを含め、cache hit時にTSV内の全handleが
  指定bindingに属することを検証する。display binding hashはraw keyへ含めない。
- LOSAT FASTA、raw text、combined `protein_map`、converter input、derived protein referenceを
  runtime handleへ統一する。
- derived schema 3 keyにraw key、runtime/display mapping、view/filter/converter settingを含める。
- `(query, subject)` はdirectional keyとし、初回実装ではreverse pairをdirect hit扱いしない。
  将来swap reuseを行う場合もTSV列swapと全metadata検証を別contractにする。
- verified migration を pure mapping/validation と runtime orchestration に分ける。
- Load → Save → LoadをGenerate前に行っても`pending` candidateがlosslessに残ることを固定する。
- promotionはcopy-on-writeとし、schema 4 rawとschema 3 derivedのvalidation完了後だけlegacy
  candidateを除去する。
- migration failure は pair-local miss にし、cache 全体を破棄しない。
- blastn、tblastx、circular conservation の schema 2 path を regression test で維持する。

### Phase A3: Version 36 codec/reader と TSV export hydration

- Python/Webにversion 36 codec/validator/reader testを追加するが、current writer constantはまだ
  切り替えない。
- protein set、record analysis、runtime/display binding、v27–34/v35 legacy candidate envelopeを
  session artifactとしてserialize/restore/reset/snapshotする。
- v36 current cacheへprotein schema 2/3を書かず、import由来の旧entryだけをlegacy envelopeで
  round-tripすることをvalidatorで分離する。
- `Save Raw LOSAT TSV`のpair/bulk downloadをPython hydratorへ通し、QUERY/SUBJECTを通常aliasへ
  置換する。一意aliasへrecord/hash suffixを付けず、duplicate時だけ決定的ordinalを付ける。
- hydratorが12列、row順、第3～12列、comments、empty result、改行を保持し、unknown/wrong-binding
  handleでfail closedになることを確認する。
- manifestから全runtime handleを実featureへ100%解決できることを検証する。
- derived payload、orthogroup member/edge/path、selection、editor override、feature/result metadataの
  全protein referenceをruntime handleへ解決するinventory/rewriteを追加する。
- orthogroupのname/description overrideとselected alignment featureは、旧ID文字列ではなく
  resolved feature identityにより意味を保持する。
- version 27–35 read、version 36 save/load、gzip round-tripを固定する。
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
3. Python/Webのcurrent writer constantを同じchangeで36へ一度だけ切り替える。
4. Gallery refresh mergeを更新し、refreshed側のversionとLOSAT artifactsを採用できるようにする。
5. 全11 Gallery session と全2 current test input session を version 36 へ更新する。
6. legacy fixtureはcurrent writer/refresh toolへ一度も通さない。
7. session diffでversion、render schema、raw-4/derived-3/manifest-2、legacy absence、設定値をreviewする。
8. Gallery source/example SVG の Y geometry diff を review する。
9. review 後に thumbnail と必要な tutorial media だけ更新する。
10. dedicated browser acceptance runnerでLoad → Save → Load → Generateと、Generate → Save → Reload →
    Generateを確認する。
11. Vibrio sessionの59 raw entryとgzip/expanded size hard gateを確認し、entry pruningがないことを
    inventoryで固定する。
12. browser acceptanceをrequired gateとして通し、skipを成功扱いしない。
13. non-slow suite、ruff、reference comparisonを通す。

## 7. ファイル別の変更計画

### 7.1 LOSAT identity/cache

- `gbdraw/analysis/protein_colinearity.py`
  - canonical feature analysis ID、runtime handle、export ID hydration、protein set、record analysis、
    runtime/display binding、manifest/hashのsingle owner。
  - schema-2/3 ID remapとPython/CLI raw cache schema 4。
- `gbdraw/web/js/app/python-helpers.js`
  - metadata hash IDの複製を削除し、Python ownerのextraction/migration/hydrationを呼ぶ。
- `gbdraw/web/js/app/run-analysis.js`
  - File fingerprint 由来の protein identity を廃止する。
  - schema 4 lookup、v27–34/v35 lazy migration、manifest-2をorchestrationする。
  - FASTA、raw、converter、derived/UI referenceをruntime handleへ統一する。
  - LOSAT executionを小さいexecutor interface経由にし、acceptance testでworker callをcountできるようにする。
- `gbdraw/web/js/app/losat-cache.js`（新規、pure function のみ）
  - protein schema 4 / nucleotide schema 2 / derived schema 3 / manifest schema 2の
    discriminated payload、validator、lookup、promotion、legacy candidate状態機械、copy-on-write mapping。
  - 現行の共有`LOSAT_CACHE_SCHEMA`をprotein/nucleotideの別ownerへ分け、program/flowだけに依存する
    暗黙判定を散在させない。
- `gbdraw/web/js/services/config.js`
  - raw/derived/manifest/v27–34/v35 legacy candidate envelope serialization、dual-read、
    session version 36。
  - `lastModified: 0` を保持する。
- `gbdraw/web/js/services/session-authority.js`
  - `legacyArtifacts`をartifact authorityとして宣言し、canonical request promotionから分離する。
- `gbdraw/web/js/state.js`、`gbdraw/web/js/app/app-setup.js`、reset owner
  - protein set、record analysis、runtime/display binding、legacy candidateをcache stateとして
    初期化・snapshot・resetする。
- `gbdraw/web/js/services/gallery-session-migration.js`
  - render schema 2→3 promotion 中に legacy cache/evidence を落とさない。
  - protein ID の runtime 変換はここで行わない。
- `gbdraw/web/js/services/session-request.js`
  - canonical resource metadata を deterministic に materialize する。
- `gbdraw/linear.py`、`gbdraw/session_io.py`、必要に応じて `gbdraw/session.py`
  - CLI/Webのcache schema、manifest、session versionを一致させ、統合gateでversion 36へ切り替える。
- raw TSV download owner
  - pair/bulk download直前にPython hydratorを使い、hydrate後byte数で50 MiB確認を判定する。
  - user-uploaded TSVにはこの変換を適用しない。

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
  - staged validatorでversion 36、protein raw schema 4、derived schema 3、manifest schema 2を確認し、
    旧artifactへ戻ったsessionをcommit前に拒否する。
  - 旧orthogroup/editor overrideはfeature analysis identityで新artifactへ再適用する。
  - Vibrioの59 raw entry、gzip/expanded hard gateを確認し、size対策のsilent pruningを拒否する。
- `gbdraw/web/gallery/sessions/`
  - 11 sessionをversion 36へ再生成する。
- `tests/test_inputs/`
  - 2 current session fixtureをversion 36へ再生成する。
- `tests/fixtures/sessions/`（新規候補）
  - 最小化した実在 schema v2 migration fixture を current inventory と分離して保持する。
  - fixtureと同じ場所にexpected metrics JSONを置き、Node/Python adapter共通のacceptance oracleにする。
- `tests/web/losat-cache.test.mjs`（新規）
  - pure cache validator、v27–34/v35 legacy状態遷移、manifest-2 reference、
    Save-before-Generate payload、export hydration境界を検証する。
- `tests/web/losat-cache-migration.playwright.spec.js`（新規）
  - 実在fixture、counting LOSAT executor、structured telemetryを使うbrowser acceptance本体。
- `tests/run_losat_cache_browser_acceptance.py`（新規）
  - Node `@playwright/test`を優先し、なければPython Playwrightで同じfixture/expected metricsを実行する。
    両方利用不能、assertion未実行、browser起動失敗はexit 0やskipにせず失敗させる。
- `.github/workflows/test.yml`のbrowser test job
  - 上記runnerをrequired gateとして実行し、Chromiumとbrowser wheelを明示的に準備する。
- `docs/PYTHON_SESSION_COMPATIBILITY_MATRIX.md`
  - version 36、raw-4/derived-3/manifest-2、v35 legacy policyを追記する。
- `docs/LINEAR_TRACK_OCCUPANCY_LAYOUT_IMPLEMENTATION_PLAN.md`
  - 実装後、constraint composition regression の修正結果を記録する。
- `docs/TUTORIALS/4_Protein_Comparisons.md`、`docs/TUTORIALS/7_Linear_Layout.md`、
  `docs/TUTORIALS/8_Interactive_SVG_Sessions.md`、`docs/FAQ.md`、release notes
  - runtime handleを露出せず、hydrated export ID、cache reuse、Auto spacingの説明が必要な箇所だけ更新する。
- Gallery source/example SVG、thumbnail、tutorial media
  - reviewed geometry change に限定して公式生成 path で更新する。

## 8. Test plan

### 8.1 LOSAT pure/Python tests

- 128-bit runtime handleのgolden vectorと`h_[a-z2-7]{26}` grammar。
- source `protein_id`をdisplay aliasに採用するが、machine feature ID/runtime handleには使わない。
- missing/duplicate `protein_id`でも同じcanonical feature ID規則を使い、duplicate export aliasだけ
  feature analysis ID順の短いordinalで一意化する。
- Unicode/NFC、reserved delimiter、空白、control characterを含むaliasをexport時にpercent-encodeする。
- 同位置featureの永続ordinal、compound location、0-based/end-exclusive座標規則を固定する。
- 同じrecord source IDを複数配置してもsession-global handleとcombined protein mapが衝突しない。
- manifest validatorがruntimeIdsの完全coverageとglobal uniquenessをfail closedで検証する。
- 同じprotein setでも別record-instance bindingのraw cache entryをdirect hitにしない。
- 異なるrecord analysisが同じprotein setを共有してもrecord固有metadataを上書きしない。
- compound location、同一 AA の別 feature、GFF3+FASTA を扱う。
- filename/mtime/resource rename で ID/protein-set hash が不変。
- `protein_id` / `locus_tag` / GFF `ID`の変更でfeature analysis ID、runtime handle、runtime binding、
  raw keyは不変だが、display bindingとexport IDが変わる。
- AA、visible CDS、genetic code、search args の変更で raw miss。
- product/gene/note、view reverse、filter option の変更で raw hit + derived miss。
- schema 2の`p_r_` TSVとschema 3のlong-ID TSVをschema 4 runtime handleへ変換する。
- 12列TSVの第1・第2列だけをexport IDへhydrateし、row順、第3～12列、comments、empty result、
  末尾改行を保持する。
- unknown/wrong-binding handle、ambiguous/missing/corrupt legacy mappingを拒否する。
- query/subject reverse pairをdirect cache hitにしない。
- Python ownerのmanifest/hash/hydrationとCLI/Web boundary serializationがgolden fixtureで一致する。

### 8.2 Web unit/session tests

- protein schema 4 raw entry、nucleotide schema 2 raw entry、v27–34/v35 legacy candidate envelope
  schema 1、derived schema 3、manifest schema 2のdiscriminated validator。
- protein set/record analysis/runtime・display bindingのdeduplication、round-trip、size inventory。
- mixed protein/nucleotide cache を誤分類しない。
- schema 2→3 render promotion が legacy cache/evidence を保持する。
- Load → Save → LoadをGenerate前に行っても`pending` candidateがbyte-identicalに残り、current
  `losatCache`へprotein schema 2/3が混入しない。
- migration成功時だけ`pending`から`promoted`へ進み、失敗時は理由付き`rejected`を保持する。
- `lastModified: 0` が save/load 後も 0。
- version 27–35 readerとversion 36 writer。
- v35 manifest-1＋raw-3をlossless quarantineし、Generate時にworker call 0件でraw-4へ昇格する。
- current session inventoryはlegacy fixtureを除き全件version 36。
- Gallery session file list と refresh tool/test inventory が一致する。
- Gallery refresh後のversion/cache/manifestがrefreshed artifact側から採用され、旧schemaへ戻らない。
- current v36 raw/derived内に旧long transport ID、`p_r_...`、full feature hashの反復参照がない。
- raw/derived内の全handleがmanifest-2へ解決し、wrong-binding referenceを拒否する。
- derived、orthogroup member/edge/path、selection、editor overrideの全protein IDがmanifest/bindingへ
  一意に解決し、名前・説明・選択の意味が保持される。
- pair/bulk downloadにruntime handleが残らず、hydrate後byte数でsize confirmationを判定する。
- user-uploaded comparison TSVはbyte-for-byte変更しない。
- Vibrio sessionが59 raw entryを保持してgzip/expanded hard gateを満たす。

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
- alias/display metadataだけを変えてもrawは0 jobで、downloadは新しいdisplay aliasを使う。
- 一つの AA/protein set を変えると影響する pair だけ job が走る。
- Load → Save → Load → Generateでもlegacy candidateを失わず0 job。
- Generate → Save → Reload → Generateを二回行っても0 job。
- v34 schema-2とv35 schema-3の両fixtureがraw schema 4へ0 jobで移行する。
- current internal QUERY/SUBJECTはruntime handleで、manifest-2から一意なfeature/bindingへ戻る。
- current raw textに`@...|...~f_<64hex>`と`p_r_...`が残らない。
- pair/bulk downloadのQUERY/SUBJECTは通常aliasで、`h_...`、`p_r_...`、64桁feature hashを出さない。
- duplicate aliasだけが決定的ordinalを持ち、未解決handleはdownload全体をfailさせる。
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
2. Workstream A1/A2を実装し、protein set、record analysis、runtime/display binding、
   runtime handleのcontractを成立させる。
3. schema v2 BGC fixtureのLoad → Save → Load → Generateと25/25 verified reuseを成立させる。
4. Workstream B1/B2 を実装し、single/multi-row の constraint solver を統一する。
5. version 35 schema-3/manifest-1 fixtureの0-job migrationとexport hydrationを成立させる。
6. Python/Web session versionを36へ一度だけ切り替える。
7. Gallery merge toolをraw-4/derived-3/manifest-2 authorityとsize inventoryに対応させる。
8. Vibrioの59 raw entryを保持したcompact sessionでsize hard gateと0-job reuseを確認する。
9. 全current sessionを一括再生成し、legacy fixtureだけを旧session/protein形式で残す。
10. SVG/reference/Gallery visual diffをreviewする。
11. docs、release notes、必要なmediaを更新する。
12. full verification後にmetadata-dependent ID generator、legacy current-cache write path、旧spacing式を削除する。

protein schema 2/3 reader/migratorはversion 36 releaseで残す。protein schema 2/3 entryをcurrent
cacheへ書くwriterとmetadata-dependent ID generatorは残さない。import由来candidateのlossless
legacy envelope writerとnucleotide schema 2 writer/readerは維持する。rollback時もlegacy protein
artifactを破壊しないよう、lazy migrationはcopy-on-successとする。

## 10. Risks and mitigations

| Risk | Mitigation |
|---|---|
| legacy cache の false hit | full FASTA hash、AA digest、program/args、一対一 mapping をすべて検証する |
| 同じprotein setを持つ異なるrecordのmetadata上書き | `ProteinSet`からrecord固有情報を外し、record analysis/bindingで参照する |
| missing/duplicate source ID | source IDをdisplay aliasだけに使い、canonical feature analysis IDで一意化する |
| 128-bit handleのdigest collision | manifest全体でglobal uniquenessを検証し、重複をfail closedにする |
| runtime handleがUIへ露出 | user-visible境界をinventory化し、manifest resolver/hydratorを必須にする |
| export IDのescape/重複衝突 | normative NFC/UTF-8 percent encodingと決定的duplicate ordinalをtestする |
| Generate前の再保存でlegacy cacheを失う | current cacheと分離したlegacy candidate envelopeをlosslessにround-tripする |
| v35 raw schema 3をcurrent hitと誤認 | version-36 validatorとraw-4 discriminatorで専用candidateへ隔離する |
| manifestによるsession肥大化 | protein setをhashでdeduplicateし、instance bindingは1 rowにつき1件、pair entryは短いhandleだけを反復する |
| empty legacy result を誤って再利用 | record token を証明できなければ pair-local miss |
| Web/Python identity drift | feature/manifest/runtime/raw identityはPython single ownerとWeb boundaryのbyte-identical golden testで固定し、surface-local derived keyは同等のinvalidation inputsをtestする |
| alias変更で不要なLOSAT再実行 | runtime binding hashとdisplay binding hashを分離する |
| derived payloadに旧long IDが残る | recursive inventoryとmanifest resolutionでderived schema 3を検証する |
| downloadが内部handleのまま | fail-closed hydratorとpair/bulk browser download assertion |
| compact化してもsize gateを超える | raw/derived reference inventoryを測定し、entryをpruneせず59件を固定する |
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
- runtime handleを外部database IDまたはuser-visible IDとして扱わない。
- user-uploaded TSVを正規化・renameしない。
- size対策としてvalidなcache entryをpruneしない。
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
4. new FASTA、raw schema 4、derived schema 3のprotein referenceはsession-global runtime handleを使う。
5. runtime handleはalias、filename、mtime、record/display metadataに依存せず、manifest validatorが
   全件を一意なrecord instance/featureへ解決する。
6. 同じsource recordを複数配置してもruntime handleとcombined protein mapが衝突せず、別bindingの
   raw TSVをdirect hitとして誤利用しない。
7. 異なるrecord analysisが同じprotein setを共有してもrecord固有metadataとbindingを上書きしない。
8. raw schema 4の全QUERY/SUBJECTが指定query/subject runtime bindingに属する。
9. filename、mtime、resource rename、save/load、alias/display metadataだけではraw cache keyが
   変わらない。display変更はexport IDと必要なderived metadataだけを更新する。
10. AA、protein membership、feature identity、record instance、意味のあるLOSAT argsの変更では
    該当raw cacheが確実にmissする。
11. reverse query/subject pairを未検証のdirect hitとして再利用しない。
12. raw cacheとderived cacheのinvalidation boundaryがunit testで固定される。
13. `Save Raw LOSAT TSV`は通常aliasを出力し、duplicate時だけ決定的ordinalを付ける。
14. download TSVは12列、row順、第3～12列、comments、empty result、改行を保持し、未解決handleを
    ユーザーへ出力しない。
15. pair/bulk downloadに`h_...`、`p_r_...`、64桁feature hashが残らず、user-uploaded TSVを変更しない。
16. Pythonのsingle ownerがfeature/record/protein-set identity、runtime/display binding、
    runtime handle、raw key、manifest、hydrationを生成し、CLIとWeb boundaryがbyte-identicalに
    serializeする。surface-localなderived keyは同じ意味のinvalidation inputsをすべて含む。
17. current raw/derived artifact内の全protein referenceがmanifest-2へ解決し、old readable long
    transport ID、`p_r_...`、full feature hashの反復参照が残らない。
18. Load → Save → LoadをGenerate前に行ってもv27–34/v35 legacy candidateがlosslessに残り、
    current cacheへprotein schema 2/3 entryが混入しない。
19. version 35 schema-3/manifest-1 sessionがLOSAT worker 0件でraw schema 4へ移行する。
20. orthogroup名・説明override、selection、editor stateがID移行後も同じfeatureを指す。
21. axis gapはXが交差するeligible body/comparison/definition kind pairの必要間隔の`max()`であり、
    二重加算しない。
22. `comparison_height=60`はexclusion edge間の最小clear corridorとして保証される。
23. comparisonのないrow境界は`comparison_height`を予約しない。
24. Xが交差しないdefinition/body pairはspacingを増やさず、Xが交差するlocal header/body pairは
    definition clear gapを満たす。
25. definition、non-overlay body、comparisonが各collision domainで交差しない。
26. off-axis definitionを非負extentの和で過大評価しない。
27. canvasが全painted contentを含み、clipしない。
28. single/multi-row、default/custom、ribbon/curve、sparse depthの回帰testが通る。
29. Python/Webのcurrent session versionは36で一致し、通常の同梱session全13件がversion 36になる。
30. repository内の旧session envelope、protein raw schema 2/3、derived schema 1/2、manifest schema 1は
    専用migration fixtureだけに残る。
31. current sessionではprotein raw schema 4とnucleotide raw schema 2の混在を正しく検証できる。
32. Vibrio sessionが全59 raw entryを保持し、gzip `< 100,000,000 bytes`、expanded JSON
    `< 536,870,912 bytes`を満たす。size対策のsilent pruningを行わない。
33. Gallery refresh後もversion/cache/manifestが旧値へ戻らず、identity migrationだけでSVG geometryを
    変更しない。全Gallery session、SVG、thumbnail、必要なtutorial mediaがreview済みcurrent出力と
    一致する。
34. dedicated browser acceptanceがskipなしでrequired gateを通り、focused tests、non-slow suite、ruff、
    reference comparisonが成功する。

## 13. Final implementation and verification results（2026-07-25）

### 完了した実装

- Python single ownerにstable feature/record/protein-set identity、128-bit runtime handle、
  runtime/display binding hash、manifest schema 2、export ID map/TSV hydrator、directional raw schema 4
  keyを実装した。WebはPyodide helper経由で同じownerを使う。
- session version 36、`renderRequest.schema == 3`、protein raw schema 4、nucleotide raw schema 2、
  derived schema 3、manifest schema 2の境界をPython/Webに実装した。readerはversion 27～35を
  migration inputとして受理する。
- v27–34の`p_r_` raw artifactとv35 raw schema 3 / manifest schema 1 / derived schema 2をcurrent
  cacheから隔離し、Save-before-Generateでlosslessに保持して、検証済みpairだけをcopy-on-successで
  runtime handleへpromotionするpathを実装した。
- LOSATP FASTA、raw text、combined `protein_map`、derived referenceをruntime handleへ統一し、
  pair/bulk downloadをmanifestからのexport hydrationへ切り替えた。
- orthogroup、selection、editor/result metadataの旧long IDと`p_r_` referenceをcurrent runtime
  handleへ解決するartifact inventory/rewriteを実装した。
- X範囲付き`CollisionBand`、kind-pair policy、`required_axis_gap()`を実装し、single/multi-rowを
  同じsolverへ統合した。body/comparison/definitionはeligible制約の`max()`で合成し、
  comparisonのない境界にcorridorを予約しない。

### 最終検証結果

- 履歴上のBGC0000708–BGC0000713 schema-v2 sessionとv35 schema-3/manifest-1 sessionをgzip
  fixtureとして固定した。実browserのLoad → Save → Load → Generateで、25 pairが
  `25 hits / 0 misses / 0 jobs`のままcurrent schemaへpromotionされ、pair/bulk downloadが内部IDを
  出力しないことを確認した。
- browser acceptanceはPython Playwright＋Chromiumで45,094 assertionsを通過した。cancel、
  renderer error、rollback後の再Generateも含め、migration stateとcacheが成功時だけcommitされる。
  Nodeの`@playwright/test`がない環境だったため、repository guidanceどおりPython adapterを
  required browser gateとして使用した。
- 全11 Gallery sessionをtransactional refreshし、version 36、protein raw 4、derived 3、
  manifest 2へ統一した。source/result/thumbnailを同期し、generated SVG IDを正規化した
  visible geometry比較で全11例が一致した。
- Vibrio artifactは59 directional raw entryのexact setとderived coverage 59/59を保持し、
  Generate相当のcache lookupは59 hits、LOSAT worker 0 jobsだった。最終サイズはgzip
  88,887,768 bytes、expanded JSON 377,965,278 bytesで、hard gateとoperational targetを
  ともに満たした。
- current derived validatorはstrict zero-hit payloadを許可しつつ、unknown runtime reference、
  compound referenceのunknown component、malformed payloadをPython/Web双方で拒否する。
  derived keyはsurface-localだが、LOSAT mode、raw key、全collinearity parameterを含む同等の
  invalidation boundaryをfocused testsで固定した。
- Linear layoutとidentity migrationのreference output comparisonは16件すべて通過した。
- Web packaging/focused integrationは100 passed、6 skipped、post-refresh Gallery gateは3 passed、
  full non-slow suiteは1,789 passed、19 skipped、10 deselectedだった。
- `ruff check gbdraw/`は成功し、current browser wheelの準備後、isolated buildで
  `gbdraw-0.14.0b0.tar.gz`と`gbdraw-0.14.0b0-py3-none-any.whl`を正常生成した。

以上により、§12のDefinition of Doneをすべて完了とする。
