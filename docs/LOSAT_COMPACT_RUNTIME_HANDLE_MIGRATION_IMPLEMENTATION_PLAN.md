# LOSAT compact runtime handle migration implementation plan

## 0. 文書の位置づけ

- 作成日: 2026-07-24
- 状態: **次セッション向け実装計画（実装未着手）**
- 対象:
  - Web/CLI の protein LOSATP FASTA ID
  - raw LOSAT cache
  - derived LOSATP cache
  - protein identity manifest
  - session migration
  - `Save Raw LOSAT TSV`
- 関連文書:
  - `docs/LOSAT_CACHE_IDENTITY_AND_LINEAR_SPACING_REGRESSION_IMPLEMENTATION_PLAN.md`

本書は、上記既存計画のうち「readable transport IDをraw TSVの全hit行へ保存する」方針を
置き換える。linear spacing側の設計・実装は変更しない。

次セッションでは、容量上限でcache entryを間引く対症療法ではなく、machine identity、
runtime reference、user-visible IDを分離する本書の方針を実装する。

## 1. 結論

protein identityは次の三層に分離する。

1. `featureAnalysisId`
   - CDS featureの完全なmachine identity。
   - 現行の`f_<64 lowercase hex>`を維持する。
   - manifestに一度だけ保存し、raw/derivedの反復参照には使わない。
2. `runtimeHandle`
   - FASTA、LOSAT内部raw TSV、`protein_map`、derived payloadで使う短いsession-global ID。
   - record instanceとfeature analysis IDから決定的に作る。
   - display alias、filename、mtime、record source labelには依存しない。
3. `exportId`
   - `Save Raw LOSAT TSV`でユーザーへ見せる通常のprotein/feature名。
   - manifestのdisplay metadataからdownload時に生成する。
   - 内部cacheには各hit行ごとに反復保存しない。

job-localな`q0`/`s0`方式は採用しない。session-global handleなら、現行のglobal
`protein_map`とconverterを維持でき、rawだけでなくderived payload内の反復IDも短くできる。

## 2. 方針変更が必要な理由

現行のversion 35 / protein raw schema 3は、QUERY/SUBJECTに次の文字列を保存する。

```text
<record-source-id>@<record-instance-id>|<display-alias>~<feature-analysis-id>
```

例:

```text
BGC0000708@record-1|CAG38695.1~f_2f...64-hex...
```

このIDは一意で検証可能だが、同じ文字列を大量のhit行へ反復する。Vibrio gallery sessionの
実測は次の通り。

| 項目 | 実測値 |
|---|---:|
| expanded session JSON | 673,521,877 bytes |
| gzip session | 134,533,144 bytes |
| raw cache JSON | 171,600,482 bytes |
| raw entries | 59 |
| raw TSV rows | 657,743 |
| manifest transport IDs | 24,027 |
| QUERY/SUBJECT列だけの合計 | 136,078,956 bytes |

session-globalな短い連番tokenでread-only試算すると、raw cacheは
171,600,482 bytesから40,735,044 bytesへ減った。一方、rawとmanifestだけを短くしても
expanded sessionは約540 MB、gzipは約104 MBであり、現行gateをまだ超える。

したがって、短いIDはrawだけでなくderived payload、orthogroup reference、
`protein_map`にも一貫して使う必要がある。64桁feature hashはmanifestのauthorityとして
保持するが、hit行やderived edgeごとには反復しない。

## 3. 変更後のidentity model

| Layer | 役割 | 保存場所 | 反復参照 |
|---|---|---|---|
| Protein-set identity | feature＋AA集合のcontent identity | manifest | hash参照のみ |
| Feature analysis identity | record内CDSの完全なmachine identity | manifest | しない |
| Record instance identity | 同じrecordを複数配置した時の区別 | manifest/session | binding metadataのみ |
| Runtime handle | FASTA/raw/derivedの短い参照 | manifest map＋payload参照 | 許可 |
| Display metadata | protein ID、locus tag、GFF ID、location fallback | manifest | しない |
| Export ID | ユーザー向けTSVのQUERY/SUBJECT | download時に導出 | sessionへ保存しない |
| Runtime binding hash | raw cacheのbinding identity | manifest/raw key | hash参照のみ |
| Display binding hash | export/UI/derived metadataのinvalidation | manifest/derived key | hash参照のみ |

### 3.1 維持するfeature analysis ID

`featureAnalysisId`の現行contractは変更しない。

- feature type
- location operator
- 全location parts
- strand
- 永続same-location ordinal

をcanonical JSON化し、SHA-256を`f_<64 lowercase hex>`で保存する。source qualifierの
一意性には依存させない。

### 3.2 Runtime handle contract

Pythonを唯一のownerとして、次のdomain-separated payloadからhandleを作る。

```json
{
  "featureAnalysisId": "f_...",
  "recordInstanceKey": "record-1"
}
```

規範的な生成手順:

1. keyを辞書順にした空白なしUTF-8 canonical JSONを作る。
2. `b"gbdraw-runtime-handle-v1\0"`を前置する。
3. SHA-256を計算する。
4. 先頭16 bytes（128 bits）をRFC 4648 Base32でlowercase化し、paddingを除く。
5. `h_`を前置する。

grammar:

```text
h_[a-z2-7]{26}
```

必須条件:

- 同じ`recordInstanceKey`＋`featureAnalysisId`から常に同じhandleを作る。
- alias、product、gene、note、filename、mtimeを変えてもhandleは変わらない。
- 同じrecordを別instanceとして配置するとhandleは変わる。
- manifest validatorは全record instanceを通じたhandle重複を検出し、衝突時はfail closedにする。
- handleを途中で切り詰めたり、衝突時に非決定的suffixを追加したりしない。

### 3.3 Runtime binding hash

runtime binding hashはraw searchのbindingだけを表す。

```json
{
  "encoding": "runtime-handle-v1",
  "proteinSetHash": "sha256:...",
  "recordInstanceKey": "record-1",
  "runtimeIds": {
    "f_...": "h_..."
  }
}
```

`recordAnalysisId`、record source label、display alias、feature display metadataは含めない。
これにより、配列・feature membership・record instanceが同じなら、表示名だけの変更でraw
LOSATを再実行しない。

### 3.4 Display binding hash

display binding hashは、export/UI/derived payloadに埋め込まれる表示情報の変更を検出する。

最低限、次をcanonicalizeする。

- `recordAnalysisId`
- `recordSourceId`
- record instance key
- featureごとのdisplay alias
- export時のduplicate ordinal
- derived payloadが保持するfeature metadata

alias、product、gene、note等の表示情報が変わった場合:

- raw cache: hit
- derived cache: 必要に応じてmiss/rebuild
- export ID: 新しい表示へ更新

## 4. ユーザー向けTSV contract

### 4.1 Generated protein LOSAT TSV

内部raw TSVではQUERY/SUBJECTに`runtimeHandle`を保存する。

```text
h_abc...\th_def...\t99.1\t...
```

`Save Raw LOSAT TSV`では、download直前にmanifestを使って第1・第2列だけを`exportId`へ
置換する。

通常は、元レコードの次の優先順位で選んだaliasをそのまま使う。

1. `protein_id`
2. `locus_tag`
3. GFF `ID`
4. 完全location fallback

同一record instance内でaliasが一意なら、余分なrecord instanceやhashを付けない。

```text
CAG38695.1
```

同一record instance内でaliasが重複する場合だけ、`featureAnalysisId`順に決めた短い
1-based ordinalを付ける。

```text
duplicate~1
duplicate~2
```

aliasはUnicode NFC化後、現行のpercent encoding規則を適用する。export IDの一意性scopeは
raw entryのquery側record instance内、またはsubject側record instance内とする。queryとsubjectで
同じaliasを使うことは問題にしない。

hydratorは次を保証する。

- 非comment行は厳密に12列。
- row順、3～12列、数値表現、改行を変更しない。
- QUERY/SUBJECTの全handleをentry metadata＋manifestから解決する。
- 一件でも未解決、別binding、重複解決があればdownload全体を失敗させる。
- 失敗時に内部handleをfallback出力しない。
- 50 MiB確認ダイアログは、内部cache文字数ではなくhydrate後のexport byte数で判定する。

### 4.2 User-uploaded TSV

ユーザーがPairwise Comparisonsへ投入するBLAST/TSVは一切書き換えない。本変更の対象は、
gbdrawがprotein LOSATPから生成した内部raw cacheと、そのdownloadだけである。

### 4.3 Derived/UI reference

runtime handleは内部referenceとしてderived payload、orthogroup member/edge/path、
selection、editor stateで使ってよい。ただしUI表示、tooltip、download、診断メッセージでは
manifest resolverを通し、可能な限りdisplay aliasを表示する。

## 5. Target schema

| Owner | 現行 | Target | 備考 |
|---|---:|---:|---|
| Session envelope | 35 | 36 | Python/Web writerを同時更新 |
| Canonical `renderRequest` | 3 | 3 | 変更しない |
| Protein raw LOSAT cache | 3 | 4 | compact runtime handle |
| Nucleotide raw LOSAT cache | 2 | 2 | 変更しない |
| Derived LOSATP cache | 2 | 3 | protein referenceをruntime handleへ統一 |
| Protein identity manifest | 1 | 2 | runtime/display bindingを分離 |
| v27–34 protein raw candidate | 1 | 1 | 現行lazy migrationを維持 |
| v35 protein raw candidate | — | 1 | schema-3＋manifest-1をlosslessに隔離 |

version 36 writerはprotein raw schema 4、derived schema 3、manifest schema 2だけを書く。
version 35の意味を変更してschema 3を別形式として再利用しない。

### 5.1 Manifest schema 2例

```json
{
  "schema": 2,
  "proteinSets": {
    "sha256:...": {
      "schema": 1,
      "proteins": [
        {
          "featureAnalysisId": "f_...",
          "aaSha256": "...",
          "locationParts": [[7438, 8458, 1]],
          "sameLocationOrdinal": 1
        }
      ]
    }
  },
  "recordAnalyses": {
    "sha256:...": {
      "schema": 1,
      "recordSourceId": "BGC0000708",
      "proteinSetHash": "sha256:...",
      "selector": null,
      "region": null
    }
  },
  "recordInstances": {
    "record-1": {
      "schema": 2,
      "recordAnalysisId": "sha256:...",
      "runtimeBindingHash": "sha256:...",
      "displayBindingHash": "sha256:...",
      "runtimeIds": {
        "f_...": "h_..."
      },
      "featureMetadata": {
        "f_...": {
          "displayAlias": "CAG38695.1",
          "featureSvgId": "..."
        }
      }
    }
  }
}
```

`proteinSets`は従来どおりrecord非依存でdeduplicateする。`runtimeIds`と
`featureMetadata`のkey集合は、参照するprotein setの`featureAnalysisId`集合と完全一致させる。

### 5.2 Protein raw schema 4例

```json
{
  "schema": 4,
  "kind": "raw-losat",
  "identityKind": "protein",
  "idEncoding": "runtime-handle-v1",
  "key": "sha256:...",
  "text": "h_...\\th_...\\t...\\n",
  "program": "blastp",
  "outfmt": "6",
  "args": [],
  "queryProteinSetHash": "sha256:...",
  "subjectProteinSetHash": "sha256:...",
  "queryRuntimeBindingHash": "sha256:...",
  "subjectRuntimeBindingHash": "sha256:...",
  "queryRecordInstanceKey": "record-1",
  "subjectRecordInstanceKey": "record-2"
}
```

raw keyはschema/domain、direction、query/subject protein-set hash、query/subject runtime
binding hash、program、outfmt、意味のあるargsから作る。display binding hashは含めない。

### 5.3 Derived schema 3

derived keyは最低限、次を含む。

- 使用した全raw schema-4 key
- runtime mapping hash
- display mapping hash
- view transform
- filter/converter/orthogroup/collinearity settings

payload内のprotein referenceはすべてruntime handleへ統一する。完全なfeature hashや
readable long transport IDをedge/memberごとに反復しない。

## 6. Migration design

### 6.1 version 27–34 / protein raw schema 2

現行のverified lazy migrationを維持し、最終変換先をschema 4へ変更する。

```text
legacy p_r_ ID
  → full FASTA/AA/program/args/direction検証
  → featureAnalysisId
  → manifest-2 runtimeHandle
  → raw schema 4
```

既存BGC fixtureの`25 hits / 0 misses / 0 jobs`を維持する。empty resultは、現在と同様に
別entryのrecord-token evidenceと完全なmetadataで証明できる場合だけ昇格する。

### 6.2 version 35 / raw schema 3 / manifest schema 1

version 35のschema-3 entryは、旧manifest-1で検証済みのため、LOSATを再実行せず移行できる。
ただしversion 36 current cacheへblind copyしない。

import時:

1. version 35 validatorでmanifest-1、raw schema-3 key、binding、TSV全行を検証する。
2. raw schema-3 entryとmanifest-1を専用の
   `legacyArtifacts.proteinRawV35Candidates`へlosslessに隔離する。
3. current raw mapへは入れない。
4. Save-before-Generateでもcandidateとsource manifestをbyte-equivalentに保持する。

Generate時:

1. inputからcurrent manifest-2を作る。
2. old manifestの`transportIds`を逆引きし、旧QUERY/SUBJECTを`featureAnalysisId`へ解決する。
3. current manifest-2の`runtimeIds`へ置換する。
4. protein-set、AA digest、record instance、direction、program、outfmt、argsを再検証する。
5. 第3～12列とrow順が同一であることを確認する。
6. raw schema 4 keyを再計算し、copy-on-successでcurrent cacheへ昇格する。
7. 成功したcandidateだけlegacy envelopeから除く。

aliasだけが変わっていても、feature/AA/protein setが同じなら昇格を許可する。旧long IDの
readable prefixをcurrent raw identityとして使わない。

### 6.3 version 35 derived schema 2

raw keyとprotein reference semanticsが変わるため、derived schema 2はdirect hitとして再利用しない。

- version 35 import時に、必要ならlossless legacy evidenceとして隔離する。
- migrated raw schema 4が揃った後、Generateでderived schema 3を再構築する。
- LOSAT worker callは0件のままにする。
- 新しいderived schema 3のcommit成功後に旧evidenceを除去する。
- 既存の`orthogroupState`、editor override、selectionはartifact inventory resolverで
  runtime handleへ移し、意味を維持する。

genericな文字列置換だけでderived schemaを昇格してはならない。全reference inventoryと
new derived keyの検証が必要である。

### 6.4 Failure policy

- migration failureはpair-local miss。
- 一つの破損candidateのために他のverified raw entryを破棄しない。
- schema-3 long IDをschema-4 entryとしてrelabelしない。
- 未解決runtime handleをUI、download、derived current cacheへ流さない。
- false hitより明示的な再実行を優先する。

## 7. Cache invalidation matrix

| 変更 | Runtime handle | Raw cache | Derived cache | Export ID |
|---|---|---|---|---|
| filename / mtime / resource rename | same | hit | hit | same |
| save → load | same | hit | hit | same |
| `protein_id` / `locus_tag` / GFF `ID`のみ | same | hit | rebuild if display embedded | change |
| product / gene / noteのみ | same | hit | rebuild if embedded | usually same |
| record source display labelのみ | same | hit | rebuild if embedded | filename/UIのみchange |
| reverse-complement displayのみ | same | hit | rebuild view transform | same |
| AA sequence変更 | feature IDはsameの場合あり | miss via protein set | miss | same/updated |
| protein membership変更 | affected handles/set change | miss | miss | change |
| feature location/strand/ordinal変更 | change | miss | miss | change |
| record instance key変更 | change | miss | miss | duplicate scope再計算 |
| program/outfmt/search args変更 | same | miss | miss | same |
| post-filter/orthogroup option変更 | same | hit | miss | same |
| query/subject reverse | same handles | directional miss unless separately cached | miss | columns reverse |

## 8. Workstreams

### Phase 0: 現状固定と中断prototypeの整理

1. `git diff --ignore-space-at-eol`で実差分だけを確認する。
2. 中断された104 MiB cache pruning prototypeを除去する。
3. pruningとは無関係なvalidator hardening、gallery source同期、layout/reference修正は保持する。
4. oversized Vibrio sessionのraw/manifest/derived reference inventoryをfixture化する。
5. production変更前に、runtime handleをraw＋derivedへread-only変換したサイズ試算を行う。

中断prototypeが入っている可能性があるファイル:

- `gbdraw/session_io.py`
- `gbdraw/web/js/app/losat-cache.js`
- `gbdraw/web/js/services/config.js`
- `tools/refresh_gallery_sessions.py`
- `tests/test_session_io.py`
- `tests/web/losat-cache.test.mjs`
- `tests/test_web_packaging.py`

削除対象は`SESSION_LOSAT_CACHE_BYTE_LIMIT = 109_051_904`、
`pruneSerializedLosatArtifacts` / `prune_serialized_losat_artifacts`とその専用testだけである。
同じファイルにある別修正を巻き戻さない。

### Phase A: Python identity owner

主対象:

- `gbdraw/analysis/protein_colinearity.py`
- `tests/test_protein_colinearity.py`

実装:

1. `build_protein_runtime_handle()`を追加する。
2. manifest schema 2 model/validatorを追加する。
3. runtime binding hashとdisplay binding hashを分ける。
4. extractionのFASTA ID、`CdsProtein.protein_id`、`protein_map` keyをruntime handleへ変更する。
5. readable export ID mapとTSV hydratorをpure functionで実装する。
6. raw schema-4 key/validatorを実装する。
7. old schema-3→schema-4 rewrite helperを実装する。
8. derived reference inventory validatorをschema 3向けに更新する。

Web側に第二のhandle/hash生成実装を作らない。WebはPyodide helperのJSON結果を使う。

### Phase B: Web cache/runtime integration

主対象:

- `gbdraw/web/js/app/python-helpers.js`
- `gbdraw/web/js/app/losat-cache.js`
- `gbdraw/web/js/app/run-analysis.js`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/state.js`

実装:

1. Python extraction結果のmanifest-2/runtime handleをそのまま採用する。
2. JS validatorをraw-4/manifest-2/derived-3へ更新する。
3. LOSAT FASTA、raw result、combined protein map、converter inputをruntime handleへ統一する。
4. derived payloadの全protein referenceをruntime handleで検証する。
5. pair downloadとbulk downloadをhydrator経由にする。
6. download確認サイズをhydrate後UTF-8 bytesで算出する。
7. UI/telemetry/error messageでhandleを直接表示する箇所をresolver経由にする。
8. reverse pair、hidden all-vs-all pair、empty resultのdirection/binding検証を維持する。

### Phase C: Session 36 codecとmigration

主対象:

- `gbdraw/session_io.py`
- `gbdraw/web/js/services/config.js`
- `gbdraw/web/js/services/session-authority.js`
- `gbdraw/api/request_render.py`
- `gbdraw/cli_utils/session.py`
- session codec tests

実装:

1. Python/Webのcurrent session versionを36へ同時更新する。
2. supported versionsへ35を明示的に残す。
3. current writerをraw-4/derived-3/manifest-2限定にする。
4. v27–34 schema-2 candidate pathを維持する。
5. v35 schema-3＋manifest-1専用candidate envelopeを追加する。
6. Load → Save → Load before Generateのlossless round-tripを固定する。
7. successful Generate後だけcandidateをcopy-on-successで除去する。
8. CLIとWebで同じPython migration ownerを使う。

### Phase D: Derived/UI reference migration

1. orthogroup member/edge/path
2. converted pair rows
3. result metadata
4. selected orthogroup alignment feature
5. editor override
6. feature popup/reference metadata

についてscalar key、array、object keyのinventoryを作り、旧long ID、legacy `p_r_`、
runtime handleの許可scopeを明示する。

current sessionでは:

- runtime handleは内部referenceとして許可する。
- full `featureAnalysisId`の大量反復を禁止する。
- old readable transport IDと`p_r_`を禁止する。
- user-visible textではruntime handleを原則表示しない。

### Phase E: Gallery/session refresh

1. unit、codec、browser migration gateが通るまでgalleryを更新しない。
2. transactional toolで全current sessionをversion 36へ更新する。
3. Vibrio sessionが59 raw entryを保持していることを確認する。
4. derived cacheも欠落させず、runtime handleを使うschema 3として保存する。
5. visible SVG geometryがidentity変更だけで変わっていないことを比較する。
6. gallery source SVG/session/resultの同期を確認する。

cache pruningでtestを通さない。全cacheを保持したcompact representationでsize gateを満たす。

### Phase F: 既存計画・release note更新

`docs/LOSAT_CACHE_IDENTITY_AND_LINEAR_SPACING_REGRESSION_IMPLEMENTATION_PLAN.md`の次を更新する。

- identity layer表
- QUERY/SUBJECT contract
- schema matrix
- invalidation matrix
- Workstream A1/A2/A3
- browser acceptance
- risks
- Definition of Done

特に「raw TSV自体がreadable long IDを持つ」を
「内部rawはruntime handle、ユーザーexportはreadable ID」へ置き換える。

## 9. Test plan

### 9.1 Python unit

- 128-bit runtime handle golden vector
- Unicode/NFC input
- alias変更でhandle/runtime binding hashが不変
- product/gene/note変更でraw identityが不変
- record instance変更でhandleが変化
- identical recordを二回配置してもglobal handleが衝突しない
- manifestのruntimeIdsとprotein set feature集合が完全一致
- duplicate handleをfail closed
- duplicate aliasのexport ordinalが決定的
- alias一意時に余分なrecord/hashをexport IDへ付けない
- 12列TSVの第1・第2列だけをhydrate
- comments、empty result、末尾改行を保持
- unknown/wrong-binding handleでdownload変換を拒否
- raw keyがdisplay metadataに依存しない
- derived keyがdisplay mapping変更でmissする

### 9.2 Schema/session

- current version 36 writerがraw-4/derived-3/manifest-2だけを書く
- nucleotide raw schema 2とprotein raw schema 4が混在可能
- version 27–34 readerを維持
- v35 raw schema-3＋manifest-1をlossless quarantine
- v35 Load → Save → Load before Generateでcandidateを維持
- v35 Generateで0 LOSAT jobs、raw-4へ昇格
- migration後の3～12列とrow順がbyte-identical
- empty v35 raw entryを安全に移行
-破損pairだけmissし、他pairはhit
- reverse query/subjectをdirect hitにしない
- gzip round-trip

### 9.3 Web/JS

- manifest-2 structural validation
- runtimeIdsの完全coverageとglobal uniqueness
- raw-4の厳密12列/binding validation
- raw-3をcurrent cacheとして受け入れない
- derived-2をcurrent hitとして受け入れない
- pair downloadにruntime handleが残らない
- bulk downloadにruntime handleが残らない
- export IDsが通常aliasで、duplicate時だけordinal付き
- hydrated byte数で50 MiB promptを判定
- user-uploaded TSVが無変更

### 9.4 Browser acceptance

既存BGC migration fixture:

```text
totalPairs=25
cacheHits=25
cacheMisses=0
uniqueJobs=0
workerCalls=0
```

追加で確認する。

- v34 schema-2→raw-4 migration
- v35 schema-3→raw-4 migration
- Generate → Save → Reload → Generateで0 job
- Load → Save → Load → Generateでもcandidateを失わない
- current raw textに`@...|...~f_<64hex>`が反復されない
- current derived payloadにold long transport IDが残らない
- raw/derivedの全handleがmanifest-2へ解決する
- download TSVのQUERY/SUBJECTは通常alias
- download TSVに`h_...`、`p_r_...`、64桁feature hashが出ない
- preview/export SVG geometryが変更前と一致

Node PlaywrightとPython Playwright fallbackの両adapterで同じexpected fixtureを使う。

### 9.5 Size regression

Vibrio sessionで最低限次をassertする。

- raw entries: 59を維持
- hidden/self/reverse/cross-record rawも欠落しない
- gzip session: `< 100,000,000 bytes`
- expanded JSON: `< 536,870,912 bytes`
- operational target: expanded `<= 500,000,000 bytes`
- operational target: gzip `<= 90,000,000 bytes`
- raw QUERY/SUBJECT bytesが旧実測より大幅に減る
- derived protein-reference bytesも旧実測より大幅に減る

hard gateを超える場合は、まずderived/reference inventoryと重複feature catalogを測定する。
raw/derived entryのsilent pruningを再導入してはならない。

## 10. Verification commands

```bash
node --check gbdraw/web/js/app/run-analysis.js
node --check gbdraw/web/js/app/losat-cache.js
node --check gbdraw/web/js/services/config.js
node tests/web/losat-cache.test.mjs

pytest tests/test_protein_colinearity.py -v
pytest tests/test_session_io.py tests/test_api_session.py -v
pytest tests/test_web_packaging.py -v
pytest tests/test_refresh_gallery_sessions.py -v

python tests/run_losat_cache_browser_acceptance.py

pytest tests/test_output_comparison.py::TestOutputComparison -v
pytest tests/ -v -m "not slow"
ruff check gbdraw/

python -m build
```

browser testでChromium sandbox errorが出た場合は、repository guidanceどおり同じcheckを必要な
sandbox escalation付きで再実行する。Node `@playwright/test`がない場合も、Python Playwright
fallbackを実行し、skipを成功扱いしない。

## 11. 実装順序

1. 中断pruning prototypeだけを除去する。
2. old/new size inventoryとruntime-handle dry-runを固定する。
3. Python runtime handle、manifest-2、export hydratorのpure contractを実装する。
4. raw-4 key/validatorとruntime extractionを実装する。
5. Web FASTA/raw/converter/derived経路をruntime handleへ切り替える。
6. pair/bulk TSV downloadをhydratorへ切り替える。
7. session 36 writer/readerとv35 quarantineを実装する。
8. v27–34とv35のlazy promotionをraw-4へ対応させる。
9. derived/UI reference inventoryとmigrationを完了する。
10. unit/session/browser migration gateを通す。
11. Vibrioだけをtransactional refreshし、sizeと0-job reuseを確認する。
12. 全gallery sessionと必要なreferenceを更新する。
13. full non-slow、ruff、buildを実行する。
14. 既存計画書とrelease noteを最終状態へ更新する。

## 12. Acceptance criteria

次をすべて満たした場合だけ完了とする。

1. `featureAnalysisId`はmanifest authorityとして完全なSHA-256を維持する。
2. new FASTA/raw/derivedはsession-global runtime handleを使う。
3. runtime handleはalias、filename、mtime、display metadataに依存しない。
4. manifest validatorが全handleを一意なrecord instance/featureへ解決する。
5. raw schema 4の全QUERY/SUBJECTが指定query/subject bindingに属する。
6. derived schema 3の全protein referenceがmanifest-2へ解決する。
7. new current artifactにold readable long transport IDと`p_r_...`が残らない。
8. `Save Raw LOSAT TSV`は通常のaliasを出力し、duplicate時だけ短いordinalを付ける。
9. download TSVは12列、row順、3～12列を保持する。
10. 未解決handleをユーザーへ出力しない。
11. user-uploaded TSVを変更しない。
12. alias/display metadata変更だけではraw cacheがmissしない。
13. AA、protein membership、feature identity、instance、program/args変更ではrawがmissする。
14. reverse directionを未検証のdirect hitとして再利用しない。
15. v27–34 BGC fixtureが`25 hits / 0 misses / 0 jobs`になる。
16. v35 schema-3 sessionも`0 jobs`でraw-4へ移行する。
17. Load → Save → Load before Generateでlegacy candidateを失わない。
18. migration failureはpair-localで、他のverified entryを破棄しない。
19. orthogroup override、selection、editor stateが同じfeatureを指し続ける。
20. Vibrio sessionが全59 raw entryを保持してsize hard gateを通る。
21. identity migrationだけではSVG geometryが変わらない。
22. Pythonがhandle/hash/hydrationのsingle ownerで、Web境界とのgolden testが通る。
23. session 36、raw 4、derived 3、manifest 2をPython/Webで同時に書く。
24. full non-slow tests、browser acceptance、ruff、isolated buildが通る。

## 13. Non-goals

- LOSAT/BLASTの検索アルゴリズムやscoreを変更しない。
- nucleotide LOSAT cache schema 2を変更しない。
- `featureAnalysisId`を短縮してauthorityを弱めない。
- runtime handleを外部database IDとして扱わない。
- user-uploaded TSVを正規化・renameしない。
- cache entry pruningをsize対策として導入しない。
- linear spacing、circular layout、SVG X geometryを同時に変更しない。
- 汎用的な全session JSON dictionary compressionを導入しない。

## 14. Risks and mitigations

| Risk | Mitigation |
|---|---|
| truncated digest collision | 128-bitを使い、manifest全体で重複をfail closed |
| handleがUIへ露出 | user-visible境界をinventory化しresolverを必須化 |
| v35 cacheのsilent loss | 専用candidate envelopeとSave-before-Generate test |
| schema-3を誤ってcurrent hit | raw-4 discriminatorとversion-36 validator |
| derived内にlong IDが残る | recursive inventory＋manifest resolution test |
| alias変更で不要なLOSAT再実行 | runtime/display binding hashを分離 |
| downloadが内部tokenのまま | fail-closed hydratorとbrowser download assertion |
| short ID化だけでsize不足 | derived参照も同じhandleへ統一し、Vibrio hard gateで確認 |
| size gateのためcacheを欠落 | 59 raw entry countをacceptanceに固定 |
| Python/Web hash drift | Python single owner＋byte-identical golden fixture |
| broad CRLF diffに実差分が埋もれる | `git diff --ignore-space-at-eol`でreview |

## 15. 次セッション開始時の注意

- 現在のVibrio sessionは、compact handleへ未移行のためgzip約134.5 MBのままである。
- 中断された104 MiB pruning prototypeはsessionをまだ再生成していない。
- untrackedの空staging directory
  `gbdraw/web/gallery/sessions/.gbdraw-gallery-refresh-238xk35_/`
  が残っている可能性がある。空であることを確認してから`rmdir`する。
- worktreeには本件以外の意図した修正もあるため、広範なcheckout/resetを行わない。
- `tests/reference_outputs/`は通常testでは更新せず、geometry changeをreviewした場合だけ
  専用update commandを使う。

次セッションの最初の成果物は、production codeではなく次の三つのfailing testにする。

1. alias変更でもraw keyが同じになるidentity test
2. internal raw handleを通常aliasへhydrateするTSV test
3. v35 schema-3を0 jobでschema-4へ昇格するsession/browser test

この三つが設計境界を固定してからproduction実装へ進む。
